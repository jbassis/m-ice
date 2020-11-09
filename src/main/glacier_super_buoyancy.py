"""
mpm testing
"""
from importlib import reload


import meshModel
reload(meshModel)
from meshModel import *

import material
reload(material)
from material import *

from scipy.special import erfc


import tempModel
reload(tempModel)
from tempModel import *

import stokes2Dve
reload(stokes2Dve)
from stokes2Dve import *

from leopart import (
    particles,
    RandomRectangle,
    l2projection,
    advect_rk3,
    assign_particle_values,
    AddDelete,
    RandomCell
)

import time

import os

import pylab as plt

import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL').setLevel(logging.ERROR)
parameters['ghost_mode']='shared_facet'

set_log_level(50)
set_log_active(False)

# Turn on plastic failure
plastic = True


# Surface and bottom temperatures
Ts = -20.0
Tb = -20.0

melange = False
buttressing = 1e-32
buttressing_height_max = 25.0 # In meters
buttressing_height_min = 55.0 # In meters


ice_thick = 800.0
Hab = 0.0
fname_dir = 'data/cliff/water_depth_700_super_buoyancy_crevasse3/'
xyield_min = 3e3
water_depth = ice_thick*910.0/1020 - Hab
# Set length of domain and water depth
length= ice_thick*12
#water_depth = ice_thick*910.0/1020 - Hab



# Set mesh resolution and estimate approximate number of points in x/z dir
dz = round(ice_thick/13.333333333/2)
Nx = int(length/dz)
Nz = int(ice_thick/dz)

# Define geometry of domain
surf_slope =  0.02
bed_slope =   0.04
left_vel = 3e3/material.secpera*material.time_factor


L = length+ice_thick*0
bump_width = ice_thick
bump_height =  0.25*ice_thick*0

def bot_fun(x):
   b = -water_depth + bed_slope*(length-x)
   return b

def bed_fun(x):
   b = -water_depth + bed_slope*(length-x) + bump_height*exp(-((L-x)/bump_width)**2)
   return b

def bed_fun_np(x):
   b = -water_depth + bed_slope*(length-x) + bump_height*np.exp(-((L-x)/bump_width)**2)
   return b

def ice_thick_fun(x):
    return ice_thick + (surf_slope-bed_slope)*(L-x)

def surf_fun(x):
   s = -water_depth + ice_thick + (surf_slope)*(L-x)#-  notch_slope*(L+notch_length-x)*(x<(L+notch_length))*(x>(L)) + notch_slope*(L-x-notch_length)*(x>(L-notch_length))*(x<=(L))
   return s

# Initialize mesh
start = time.time()
mesh = MeshModelPoly(surf_fun,bed_fun_np,bed_fun_np,Nx=Nx,Nz=Nz,length=length,dz=dz)
mesh.length = length*1.1
print('Time to make mesh and basis functions',time.time()-start)


# Add crevasse to mesh
dz = mesh.dz
bmesh = BoundaryMesh(mesh.mesh,'exterior',order=True)
pt = sort_boundary_nodes(bmesh)


# Remove two points that must be determined by manual inspection
#pts = vstack((pt[0:380,:],pt[393::,:]))## Identify which points to remove
#pts = pt
#pts = vstack((pt[0:555,:],pt[557:571,:],pt[581::,:]))## Identify which points to remove
idx = 326;
xpt1 = pt[idx,0];
ypt1 = pt[idx,1];
xpt2 = pt[idx+1,0];

xc = 0.5*(xpt1+xpt2)
yc = -200.0
pts = np.vstack((pt[0:idx,:],[xc,yc],pt[idx+1::,:]))## Identify which points to remove


#pts = vstack((pt[0:326,:],pt[326::,:]))## Identify which points to remove
#pts = pt


# Make a new mesh based on remeshing stuff . . .
# Now remove nodes that are plast the cutoff length and
pt_new = []
pt_flag = None
length_flag = True
xcliff = 1e6

for n in range(len(pts)):
    pt = pts[n]
    # We will stack x points along the calving front if they exceed the distance
    if near(pt[0],0) and pt[1]<mesh.bed_fun(0.0):
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
mesh.mesh=meshGmsh(pt_new.transpose(),dz*2)
mesh.mesh.bounding_box_tree().build(mesh.mesh)
mesh.generate_function_spaces()

#_____________________________________________
# Create particles defined on mesh
#xm,zm = mesh.get_coords();xmin = np.min(xm); xmax=np.max(xm);ymin = np.min(zm); ymax=np.max(zm)
#xp = RandomRectangle(Point(xmin, ymin), Point(xmax, ymax)).generate([Nx*10, Nz*10])
p_min = 8
p_max = 16

#p_min = 16
#p_max = 32
print('Generating particles')
gen = RandomCell(mesh.mesh)
xp = gen.generate(p_min)
print('Done generating particles')
#_____________________________________________
# Define function space for strain and temperature function
Vdg = FunctionSpace(mesh.mesh, 'DG',1)

print('Initializing function spaces')
strain_mesh, temp_mesh, epsII_mesh = Function(Vdg), Function(Vdg), Function(Vdg)

# Initial conditions for strain and temp
strain_fun = Expression("0.0", degree=1)
temp_fun = temp_init(Ts,Tb,surf_fun, bed_fun,degree=1)

# Create functions defined on mesh with initial values
strain_mesh.assign(strain_init)
temp_mesh.assign(interpolate(temp_fun,Vdg))
epsII_mesh.assign(strain_init)

# Particle values at nodes
pstrain = assign_particle_values(xp, strain_fun)
pepsII = assign_particle_values(xp, strain_fun)
ptemp = assign_particle_values(xp, temp_fun)
print('Done initializing function spaces')


# Now we initialize the particle class
print('Creating particles')
p = particles(xp, [pstrain,ptemp,pepsII], mesh.mesh)
print('Done particles')


# Make sure we have enough particles per cell
#AD = AddDelete(p, p_min, p_max, [strain_mesh, temp_mesh,epsII_mesh]) # Sweep over mesh to delete/insert particles
#AD.do_sweep()

(xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
    p. return_property(mesh , 1) ,
    p. return_property(mesh , 2),
    p. return_property(mesh , 3))
#_____________________________________________
# Initialize temperature model
print('Initializing temperature')
Tmodel = tempModel(mesh.mesh,Tb=Tb,Ts=Ts)
x = SpatialCoordinate(mesh.mesh)
zb = bot_fun(x[0])
zs = surf_fun(x[0])
Tmodel.set_mesh(mesh.mesh)
Tmodel.set_temp(x,surf_fun,bot_fun,bed_fun)
print('Done initializing temperature')

#_____________________________________________
# Viscosity and material properties
# Define the rheology and make sure we set some of the properties
glenVisc = glenFlow(grain_size=1e-3)
glenVisc.yield_strength=0.75e6#-12.5e3
glenVisc.visc_min = 5e9
glenVisc.visc_max = 1e18
glenVisc.plastic = plastic
glenVisc.yield_min = 20e3
glenVisc.crit_strain = 0.1
glenVisc.mu = 0.0

#_____________________________________________
# Viscosity and material properties
# Set inflow velocity of the domain
right_vel = 0.0 # Outflow velocity is not used

#_____________________________________________
# Define our model and some parameters
model = Stokes2D(mesh,glenVisc,left_vel=left_vel,right_vel=right_vel)
model.tempModel=Tmodel
model.tolerance = 1e-5
model.u_k = None
model.maxit = 25
model.maxit_local = 25
model.local_err_min = 1.0
model.tracers = particles

model.calving_front = True # This applies a stress boundary condition to the right edge
model.alpha = 1.0  # Not used anymore, legacy from when we used overrelaxation
model.water_drag = 1e5 # Quadratic drag term associated with water drag
# Add lateral drag to model
B = 0.75e8
width = 10e3
model.lateral_drag = 2*(4)**(1.0/3)*B/(width**(4.0/3))*0
model.bed_yield_stregth = 400e3
model.friction = 4e6#/(surf_fun(0)-bed_fun(0))/model.rho_i/model.g # Friction coefficient
model.m = 1.0/3.0 # Friction exponent

model.left_wall = 0.0
if melange == True:
    model.right_wall = 5.5e3
else:
    model.right_wall = 1e6
#_____________________________________________
# Maximum time step
time_step_secs = 86400.0/16*2# Time step in seconds

time_step = time_step_secs/material.time_factor # Convert time step to unit we are using


#_____________________________________________
# Start loop for simulation
# Initialize time vector to zero
t = 0.0
it_type = 'Picard'

max_length = 1.375*length# Regrid if length exceeds this value
min_length = max_length-ice_thick # Set new length after regridding to this value
model.mesh.length = max_length # Set this as the max length of the mesh--doesn't actually do anything
save_files = True # Set to True if we want to save output files
#fname_base = fname_dir + 'glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(left_vel/1e3)+'_high_res_T_'+str(Tb)+'_buttressing'+str(buttressing/1e3)+'kPa_CFL/'
fname_base = fname_dir + 'glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(left_vel/1e3)+'_high_res_T_'+str(Tb)+'_CFL/'

if not os.path.exists(fname_base):
    os.makedirs(fname_base)


if save_files==True:
    import shutil
    shutil.copy2('glacier_test_buoyancy.py', fname_base+'glacier_test_buoyancy.py')
    shutil.copy2('stokes2Dve.py', fname_base+'stokes2Dve.py')


input_flux = left_vel*(surf_fun(0.0)-bot_fun(0.0)) # Define input flux at left edge of the domain
CFL = 1.0
model.u_k = None
tau = 0.1*60*60/(60*60*24*365.24) # Relaxation time for upstream plastic strain
i =0
model.strain = Function(model.mesh.Q)
L2_strain = []
tlist = []
model.method = 1
model.buttressing = buttressing # Apply buttressing force (in Pa) in x-direction to portion of ice under water???
model.buttressing_height_max = buttressing_height_max # Apply buttressing force (in Pa) in x-direction to portion of ice under water???
model.buttressing_height_min = buttressing_height_min # Apply buttressing force (in Pa) in x-direction to portion of ice under water???

model.buttressing = 1e-32
for i in range(i,100000):
   #L2_strain.append(assemble(model.strain*dx(model.mesh.mesh))/assemble(Constant(1.0)*dx(model.mesh.mesh)))
   #L2_strain.append(assemble(model.strain*dx(model.mesh.mesh)))
   #tlist.append(t)
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


   print('Time step',time_step)
   u,pres = model.solve(p,dt=time_step,tolerance=model.tolerance)

   # Do a little bit of accounting to make sure that our time step doesn't violate the CFL criterion
   ux, uz = model.get_velocity();speed = np.sqrt(ux**2+uz**2)
   Q0 = FunctionSpace(model.mesh.mesh, "DG", 0)
   time_step_CFL = CFL*np.min(project(CellDiameter(model.mesh.mesh)/sqrt(inner(u,u)),Q0).compute_vertex_values())
   print('Time step CFL',time_step_CFL)
   #if time_step_CFL<0.5*time_step:
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

   #if model.mesh.mesh.hmax()/model.mesh.mesh.hmin() > 10.0:
    #   remesh_elastic=False
     #  model.mesh.length=min_length

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
          temp_file.write(model.strain,'strain')
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
   Vdg = FunctionSpace(model.mesh.mesh, 'DG',1)
   ap = advect_rk3(p, model.vector2, model.u_k, "open")
   ap.do_step(time_step)
   AD = AddDelete(p, p_min, p_max, [interpolate(model.strain,Vdg), interpolate(model.temp,Vdg) , interpolate(model.epsII,Vdg)]) # Sweep over mesh to delete/insert particles
   AD.do_sweep()

   # Plotting
   (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
       p. return_property(mesh , 1) ,
       p. return_property(mesh , 2),
       p. return_property(mesh , 3))
   pstrain[xp[:,0]<xyield_min]=0.0
   p.change_property(pstrain,1)


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
       if melange == True:
           plt.gca()
           plt.axvline(model.right_wall,linestyle='--',color='k')
   ax2.spines['bottom'].set_visible(True)
   ax2.set_xticks([0,3e3,10*ice_thick,max_length])
   ax2.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])

   ax1.spines['bottom'].set_visible(True)
   ax1.set_xticks([0,3e3,10*ice_thick,max_length])
   ax1.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])
   plt.pause(1e-16);
   print('Time step',time_step,'Mesh quality',model.mesh.mesh.hmax()/model.mesh.mesh.hmin(),'quality ratios',quality,'number of negative epsII',sum(pepsII<0),'Percent yielded',np.sum(pstrain>0)/len(pstrain),'Maximum strain',np.max(pstrain))

   # Print some diagnostics to screen for debugging purpose
   t = t+time_step
   print('*******************************************')
   print('Time:  ',t*material.time_factor/material.secpera,'Time step',time_step)
   print(np.max(u.compute_vertex_values()))
   print('*******************************************')
   if t>1.01:
       break


#if model.method==1:
#    plt.figure(4);plt.plot(tlist,L2_strain,label=time_step_secs)
#else:
#    plt.figure(4);plt.plot(tlist,L2_strain,'--',label=time_step_secs)
#plt.legend()
