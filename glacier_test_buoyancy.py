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

from tracers_cython import *

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

import pylab as plt

import logging
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL').setLevel(logging.ERROR)


set_log_level(50)
set_log_active(False)

# Turn on plastic failure
plastic = True

# Thick cliff
cliff_type = 3

# Surface and bottom temperatures
Ts = -20.0
Tb = -20.0

# Geometric variables
# Case 1: Medium cliff
if cliff_type == 1:
    ice_thick = 400.0
    Hab = 45.0
# Case 2: Dry cliff
elif cliff_type == 2:
    ice_thick = 135.0
    Hab = ice_thick
# Case 3: Thick cliff
elif cliff_type == 3:
    ice_thick = 800.0
    Hab = 25.0

# Set length of domain and water depth
length= ice_thick*12
water_depth = ice_thick*910.0/1020 - Hab

# Set mesh resolution and estimate approximate number of points in x/z dir
dz = round(ice_thick/13.333333333)
Nx = int(length/dz)
Nz = int(ice_thick/dz)

# Define geometry of domain
surf_slope =  0.02
bed_slope =   0.01

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

#_____________________________________________
# Create particles defined on mesh
#xm,zm = mesh.get_coords();xmin = np.min(xm); xmax=np.max(xm);ymin = np.min(zm); ymax=np.max(zm)
#xp = RandomRectangle(Point(xmin, ymin), Point(xmax, ymax)).generate([Nx*10, Nz*10])
tracers_per_cell = 4
gen = RandomCell(mesh.mesh)
xp = gen.generate(tracers_per_cell)
#_____________________________________________
# Define function space for strain and temperature function
Vdg = FunctionSpace(mesh.mesh, 'DG',1)

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

# Now we initialize the particle class
p = particles(xp, [pstrain,ptemp,pepsII], mesh.mesh)

# Make sure we have enough particles per cell
p_min = tracers_per_cell
p_max = 12
AD = AddDelete(p, p_min, p_max, [strain_mesh, temp_mesh,epsII_mesh]) # Sweep over mesh to delete/insert particles
AD.do_sweep()

(xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
    p. return_property(mesh , 1) ,
    p. return_property(mesh , 2),
    p. return_property(mesh , 3))
#_____________________________________________
# Initialize temperature model
Tmodel = tempModel(mesh.mesh,Tb=Tb,Ts=Ts)
x = SpatialCoordinate(mesh.mesh)
zb = bot_fun(x[0])
zs = surf_fun(x[0])
Tmodel.set_mesh(mesh.mesh)
Tmodel.set_temp(x,surf_fun,bot_fun,bed_fun)

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
left_vel = 1e3/material.secpera*material.time_factor
right_vel = None # Outflow velocity is not used

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

#_____________________________________________
# Maximum time step
time_step_secs = 86400.0  # Time step in seconds
time_step = time_step_secs/material.time_factor # Convert time step to unit we are using


#_____________________________________________
# Start loop for simulation
# Initialize time vector to zero
t = 0.0
it_type = 'Picard'

max_length = 1.375*length# Regrid if length exceeds this value
min_length = max_length-ice_thick # Set new length after regridding to this value
model.mesh.length = max_length # Set this as the max length of the mesh--doesn't actually do anything
save_files = False# Set to True if we want to save output files
fname_base = 'data/cliff/water_depth_700/glacier_min_yield_20_surf_slope_0.02_med_visc_bed_slope_0.01_flux_1.0_CFL/'
if save_files==True:
    import shutil
    shutil.copy2('glacier_test_buoyancy.py', fname_base+'glacier_test_buoyancy.py')
    shutil.copy2('stokes2Dve.py', fname_base+'stokes2Dve.py')


input_flux = left_vel*(surf_fun(0.0)-bot_fun(0.0)) # Define input flux at left edge of the domain
CFL = 0.2
model.u_k = None
tau = 0.1*60*60/(60*60*24*365.24) # Relaxation time for upstream plastic strain
i =0
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

   if model.mesh.mesh.hmax()/model.mesh.mesh.hmin() > 10.0:
       remesh_elastic=False
       model.mesh.length=min_length

   quality=np.min(MeshQuality.radius_ratios(model.mesh.mesh).array())
   if quality<0.1:
       remesh_elastic=False
       model.mesh.length=min_length

   #remesh_elastic = False
   print(remesh_elastic)



   # Now plot some stuff
   #xm,zm = particles.get_coords()
   """
   plt.figure(1);plt.clf()
   title_str = 'Time: '+str(round(t*material.time_factor/material.secpera,2)).zfill(4)+ ' a'


   # First plot is total plastic strain
   ax1=plt.subplot(3,1,1)
   if glenVisc.plastic==True:
       #c=plt.scatter(xm,zm,s=1,c=np.log10(model.deps_m+1e-16),vmin=-10,vmax=0)
       #cbar1=plt.colorbar(c)
       c=plt.scatter(xm,zm,s=1,c=np.log10(np.maximum(particles.tracers['Strain'],1e-16)),vmin=-4,vmax=1);cbar1=plt.colorbar(c);plt.axis('equal');cbar1.set_ticks([-4,1])
       cbar1.set_label('$\epsilon_p$', fontsize=12)
       plt.axis('equal')
       cbar1.set_ticks([-4,1])
       cbar1.set_label('$\Delta \epsilon_p$', fontsize=12)
       plt.title(title_str)
       plt.xlim([0,max_length])


   # Second plot is effective strain rate
   ax3=plt.subplot(3,1,2)
   epsII_m = model.epsII_m#model.tracers.nodes_to_tracers(model.epsII)
   c=plt.scatter(xm,zm, s=1,c=np.log10(epsII_m/material.time_factor+1e-16),vmin=-13,vmax=-7);cbar3=plt.colorbar(c);plt.axis('equal');cbar3.set_ticks([-13,-7])
   cbar3.set_label('$\epsilon_{II}$ (s$^{-1}$)',fontsize=12)
   plt.xlim([0,max_length])


   # Second plot is effective viscosity
   ax2=plt.subplot(3,1,3)
   #speed_m = np.sqrt(particles.tracers['ux']**2+particles.tracers['uz']**2)
   #c=plt.scatter(xm,zm, s=1,c=speed_m,vmin=0,vmax=2*10**4);cbar2=plt.colorbar(c);plt.axis('equal');cbar2.set_ticks([0,2*10**4])
   #cbar2.set_label('speed', fontsize=12)
   #c=plot(Tmodel.Temp-273.15)
   #cbar2 = plt.colorbar(c)
   #cbar2.set_label('Temp', fontsize=12)

   plt.xlabel('Distance (km)')
   plt.xlim([0,max_length])
   xx=np.linspace(0,length*1.5,101)
   xs=np.linspace(0,length,101)
   b = bot_fun(xx)
   Habx = bot_fun(xx)-bot_fun(xx)*1020.0/920.0
   for ax in [ax1,ax2,ax3]:
       ax.set_xticks([])
       ax.set_yticks([])
       ax.spines['right'].set_visible(False)
       ax.spines['top'].set_visible(False)
       ax.spines['left'].set_visible(False)
       ax.spines['bottom'].set_visible(False)
       ax.plot(xx,bot_fun(xx),color='brown',linewidth=2)
       ax.plot(xx,bed_fun_np(xx),'--k',linewidth=2)
       ax.plot(xx,Habx,'--b',linewidth=2)
       ax.plot(xs,surf_fun(xs),'--',color='gray')
   ax2.spines['bottom'].set_visible(True)
   ax2.set_xticks([0,3e3,10*ice_thick,max_length])
   ax2.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])

   ax1.spines['bottom'].set_visible(True)
   ax1.set_xticks([0,3e3,10*ice_thick,max_length])
   ax1.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])

   ax3.spines['bottom'].set_visible(True)
   ax3.set_xticks([0,3e3,10*ice_thick,max_length])
   ax3.set_xticklabels([0,3,10*ice_thick/1e3,max_length/1e3])
   plt.xlim([0,max_length])

   plt.draw();plt.pause(1e-16);plt.show()
   if np.mod(i,10)==0:
       #particles.tracers['Strain'][xm<3.5e3]=0.0
       if save_files == True:
           #uxm = particles.tracers['ux']
           #uzm = particles.tracers['uz']
           eta_visc_m = model.eta_visc_m
           eta_plas_m = model.eta_plas_m
           eta_m = model.eta_m
           #eta_visc_m = model.tracers.nodes_to_tracers(model.eta_visc)
           #eta_plas_m = model.tracers.nodes_to_tracers(model.eta_plas)
           #eta_m = model.tracers.nodes_to_tracers(model.eta)
           fname =fname_base + 'glacier_cliff_'+str(i).zfill(3)+'.npz'
           print(fname)
           np.savez(fname, t=t, xm=xm,zm=zm,speed=speed,ux=ux,uz=uz,strain=particles.tracers['Strain'],
                epsII=epsII_m/material.time_factor,eta=eta_m*material.time_factor,
                eta_visc_m=eta_visc_m*material.time_factor,
                eta_plas_m=eta_plas_m)

           mesh_file_name =fname_base + 'glacier_cliff_'+str(i).zfill(3)+'.xml'
           mesh_file = File(mesh_file_name)
           mesh_file << model.mesh.mesh
           temp_file_name = fname_base + 'temp_'+str(i).zfill(3)+'.hdf'
           temp_file=HDF5File(mpi_comm_world(), temp_file_name, 'w')
           temp_file.write(Tmodel.Temp,'Temp')
           temp_file.write(u,'u')
           temp_file.close()


   """
   # Update all quantities (need to update this to simplify it and use RK4 if specified)
   time_step = model.update(u,time_step,p,remesh_elastic=remesh_elastic)
   (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
   p. return_property(mesh , 1) ,
   p. return_property(mesh , 2),
   p. return_property(mesh , 3))


   if remesh_elastic == False:

       Vdg = FunctionSpace(model.mesh.mesh, 'DG',1)
       Vcg = FunctionSpace(model.mesh.mesh, 'DG',1)
       (xp , pstrain , ptemp, pepsII) = (p. return_property(mesh , 0) ,
          p. return_property(mesh , 1) ,
          p. return_property(mesh , 2),
          p. return_property(mesh , 3))
       del p
       p = particles(xp, [pstrain,ptemp,pepsII], model.mesh.mesh)


   AD = AddDelete(p, p_min, p_max, [interpolate(model.strain,Vdg), interpolate(model.temp,Vdg) , interpolate(model.epsII,Vdg)]) # Sweep over mesh to delete/insert particles
   AD.do_sweep()

   plt.figure(1);plt.clf();
   plt.subplot(2,1,1);
   c=plt.scatter(xp[:,0],xp[:,1],s=0.1,c=np.log10(np.maximum(pstrain,1e-16)),vmin=-4,vmax=1);plt.colorbar(c);plt.axis('equal')
   plt.subplot(2,1,2)

   c=plt.scatter(xp[:,0],xp[:,1],s=0.1,c=log10(pepsII+1e-16),vmin=-8,vmax=-5);plt.colorbar(c);plt.axis('equal')
   plt.pause(1e-16);
   print('Time step',time_step,'Mesh quality',model.mesh.mesh.hmax()/model.mesh.mesh.hmin(),'quality ratios',quality,'number of negative epsII',sum(pepsII<0))


   print ('coordinates of first particle',xp[0],np.max(pstrain),np.min(pstrain))
   print ('coordinates of first particle',xp[0],np.max(ptemp),np.min(ptemp))

   # Print some diagnostics to screen for debugging purpose
   t = t+time_step
   print('*******************************************')
   print('Time:  ',t*material.time_factor/material.secpera,'Time step',time_step)
   print(np.max(u.compute_vertex_values()))
   print('*******************************************')
