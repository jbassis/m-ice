

dz = model.mesh.dz
bmesh = BoundaryMesh(model.mesh.mesh,'exterior',order=True)
pt = sort_boundary_nodes(bmesh)


# Remove two points that must be determined by manual inspection
pts = vstack((pt[0:436,:],pt[439::,:]))## Identify which points to remove
#pts = pt
#pts = vstack((pt[0:483,:],pt[294:394,:],pt[397::,:]))## Identify which points to remove

#pts = vstack((pts[0:369,:],pts[414::,:]))## Identify which points to remove


# Make a new mesh based on remeshing stuff . . .
# Now remove nodes that are plast the cutoff length and
pt_new = []
pt_flag = None
length_flag = True
xcliff = max_length
xcliff = 12e3

for n in range(len(pts)):
    pt = pts[n]
    # We will stack x points along the calving front if they exceed the distance
    if near(pt[0],0) and pt[1]<model.mesh.bed_fun(0.0):
        #if len(pt_new)==0:
            pt_new.append(pt)
    else:
        if pt[0]<=xcliff:
            if len(pt_new)==0:
                pt_new.append(pt)
            else:
                # If there is at least one point, we calculate the distance
                # between the new and old point
                #dist = np.sqrt((pt[0]-pt_new[-1][0])**2+(pt[1]-pt_new[-1][1])**2)
                pt_new.append(pt)

pt_new = np.array(pt_new)

# The characteristic length is the radius so twice the mesh size
mesh=meshGmsh(pt_new.transpose(),dz*2)


model.mesh.mesh=mesh
model.mesh.mesh.bounding_box_tree().build(model.mesh.mesh)
model.mesh.generate_function_spaces()

model.set_mesh(model.mesh)
model.u_k = interpolate(model.u_k,model.vector2)
model.tempModel.set_mesh(model.mesh.mesh)


Vdg = FunctionSpace(model.mesh.mesh, 'DG',1)
Vcg = FunctionSpace(model.mesh.mesh, 'DG',1)
(xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
  p. return_property(mesh , 1) ,
  p. return_property(mesh , 2),
  p. return_property(mesh , 3))
del p
p = particles(xp, [pstrain,ptemp,pepsII], model.mesh.mesh)
AD = AddDelete(p, p_min, p_max, [interpolate(model.strain,Vdg), interpolate(model.temp,Vdg) , interpolate(model.strain,Vdg)]) # Sweep over mesh to delete/insert particles
AD.do_sweep()


for i in range(i,100000):
   # First need to interpolate tracer quantities to nodes
   #node_vars = particles.tracers_to_nodes()

   # For first time step, we use a small time step.  After that, we use the CFL criterion
   if i==0:
       time_step = time_step_secs/material.time_factor
   else:
       Q0 = FunctionSpace(model.mesh.mesh, "DG", 0)
       quality=np.min(MeshQuality.radius_ratios(model.mesh.mesh).array())
       time_step = CFL*np.min(project(CellDiameter(model.mesh.mesh)/sqrt(inner(u,u)),Q0).compute_vertex_values())
       time_step = np.minimum(time_step_secs/material.time_factor,time_step)
   #time_step = time_step_secs/material.time_factor


   print('Time step',time_step)
   u,pres = model.solve(p,dt=time_step,tolerance=model.tolerance)

   # Do a little bit of accounting to make sure that our time step doesn't violate the CFL criterion
   ux, uz = model.get_velocity();speed = np.sqrt(ux**2+uz**2)
   Q0 = FunctionSpace(model.mesh.mesh, "DG", 0)
   #time_step_CFL = CFL*np.min(project(CellDiameter(model.mesh.mesh)/sqrt(inner(u,u)),Q0).compute_vertex_values())
   #print('Time step CFL',time_step_CFL)
   #if time_step_CFL<time_step:
   #   time_step=time_step_CFL
   #  u,pres = model.solve(node_vars,dt=time_step,tolerance=model.tolerance);#model.u_k = None
   #xm,zm = particles.get_coords()

   #particles.tracers['Strain'][xm<1e3]=particles.tracers['Strain'][xm<1e3]/(1.0+time_step/tau)

   # Decide if we need to remesh
   x,z = model.mesh.get_coords()
   if np.max(x)<max_length:
       remesh_elastic=True
       model.mesh.length = max_length
   else:
       remesh_elastic=False
       model.mesh.length=min_length

   if model.mesh.mesh.hmax()/model.mesh.mesh.hmin() > 10.0:
       remesh_elastic=False
       model.mesh.length=min_length

   quality=np.min(MeshQuality.radius_ratios(model.mesh.mesh).array())
   if quality<0.1:
       remesh_elastic=False
       model.mesh.length=min_length

   #remesh_elastic = False
   print(remesh_elastic)
   if np.mod(i,10)==0:
      #particles.tracers['Strain'][xm<3.5e3]=0.0
      if save_files == True:
          #uxm = particles.tracers['ux']
          #uzm = particles.tracers['uz']
          (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
            p. return_property(mesh , 1) ,
            p. return_property(mesh , 2),
            p. return_property(mesh , 3))

          #eta_visc_m = model.tracers.nodes_to_tracers(model.eta_visc)
          #eta_plas_m = model.tracers.nodes_to_tracers(model.eta_plas)
          #eta_m = model.tracers.nodes_to_tracers(model.eta)
          fname =fname_base + 'glacier_cliff_'+str(i).zfill(3)+'.npz'
          print(fname)
          np.savez(fname, t=t, xm=xp[:,0],zm=xp[:,1],speed=speed,ux=ux,uz=uz,strain=pstrain,
               epsII=pepsII/material.time_factor,temp=ptemp)

          mesh_file_name =fname_base + 'glacier_cliff_'+str(i).zfill(3)+'.xml'
          mesh_file = File(mesh_file_name)
          mesh_file << model.mesh.mesh
          temp_file_name = fname_base + 'temp_'+str(i).zfill(3)+'.hdf'
          temp_file=HDF5File(model.mesh.mesh.mpi_comm(), temp_file_name, 'w')
          temp_file.write(u,'u')
          temp_file.close()


   # Update all quantities (need to update this to simplify it and use RK4 if specified)
   time_step = model.update(u,time_step,p,remesh_elastic=remesh_elastic)
   yr = int(np.mod(t,365))
   #day = round((t - yr)*365,1)
   day = (t - yr)*365
   hr,day = np.modf(day)
   hr = round(hr*24,0)
   title_str = 'Time: '+str(yr).zfill(1)+ 'a '+str(int(day)).zfill(1)+'d '+str(int(hr)).zfill(1)+'hr'



   if remesh_elastic == False:

       Vdg = FunctionSpace(model.mesh.mesh, 'DG',1)
       Vcg = FunctionSpace(model.mesh.mesh, 'DG',1)
       (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
          p. return_property(mesh , 1) ,
          p. return_property(mesh , 2),
          p. return_property(mesh , 3))
       del p
       p = particles(xp, [pstrain,ptemp,pepsII], model.mesh.mesh)

   else:
       p.relocate()

   # Advect particles -- Turn this on to advect particles now that it is removed from Stokes script
   #p.relocate()
   ap = advect_particles(p, model.vector2, model.u_k, "open")
   #ap = advect_particles(p, model.vector2, u,model.facet_marker)
   ap.do_step(time_step)
   AD = AddDelete(p, p_min, p_max, [interpolate(model.strain,Vdg), interpolate(model.temp,Vdg) , interpolate(model.epsII,Vdg)]) # Sweep over mesh to delete/insert particles
   AD.do_sweep()

   # Plotting
   (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
       p. return_property(mesh , 1) ,
       p. return_property(mesh , 2),
       p. return_property(mesh , 3))

   xx=np.linspace(0,length*1.5,101)
   xs=np.linspace(0,length,101)

   plt.figure(1);plt.clf();
   ax1=plt.subplot(2,1,1);
   c=plt.scatter(xp[:,0],xp[:,1],s=0.1,c=np.log10(np.maximum(pstrain,1e-16)),vmin=-4,vmax=1);cbar1=plt.colorbar(c);
   cbar1.set_ticks([-4,1])
   plt.axis('equal')
   plt.title(title_str)
   plt.xlim([0,max_length])
   ax2=plt.subplot(2,1,2)
   c=plt.scatter(xp[:,0],xp[:,1],s=0.1,c=np.log10(pepsII/material.time_factor+1e-16),vmin=-8,vmax=-5);cbar2=plt.colorbar(c);plt.axis('equal')
   #plt.plot(xx,bed_fun_np(xx),'--k',linewidth=2)
   #plt.plot(xs,surf_fun(xs),'--',color='gray')
   cbar2.set_ticks([-8,-5])
   #plt.title(title_str)
   plt.xlim([0,max_length])
   plt.xlabel('Distance (km)')
   for ax in [ax1,ax2]:
       ax.set_xticks([])
       ax.set_yticks([])
       ax.spines['right'].set_visible(False)
       ax.spines['top'].set_visible(False)
       ax.spines['left'].set_visible(False)
       ax.spines['bottom'].set_visible(False)
       ax.plot(xx,bot_fun(xx),color='brown',linewidth=2)
       ax.plot(xx,bed_fun_np(xx),'--k',linewidth=2)
       ax.plot(xs,surf_fun(xs),'--',color='gray')
   ax2.spines['bottom'].set_visible(True)
   ax2.set_xticks([0,3e3,10*ice_thick,max_length])
   ax2.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])

   ax1.spines['bottom'].set_visible(True)
   ax1.set_xticks([0,3e3,10*ice_thick,max_length])
   ax1.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])
   plt.pause(1e-16);
   print('Time step',time_step,'Mesh quality',model.mesh.mesh.hmax()/model.mesh.mesh.hmin(),'quality ratios',quality,'number of negative epsII',sum(pepsII<0))

   # Print some diagnostics to screen for debugging purpose
   t = t+time_step
   print('*******************************************')
   print('Time:  ',t*material.time_factor/material.secpera,'Time step',time_step)
   print(np.max(u.compute_vertex_values()))
   print('*******************************************')
   if t>1.0:
       break
