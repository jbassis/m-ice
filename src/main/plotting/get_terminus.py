"""
Try to identify terminus position of glaciers from input files
"""
import pylab as plt
from fenics import *

import numpy as np

from leopart import (
    particles,
    l2projection,
)

from scipy.stats import linregress

# Step 1 read file
surf_slope =  0.02
bed_slope =  -0.06
flux =1.0

step = 0

T = -20.0

fname_dir = '../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_T_'+str(T)+'_CFL/'
if T==-20.0:
    fname_dir = '../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/'
#fname_dir ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/'
#fname_dir ='../data/cliff/water_depth_0/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_CFL/'
fname_dir = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_2.0_high_res_T_-20.0_CFL/'
#fname_dir = '../data/cliff/water_depth_700_high_friction/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(-0.01)+'_flux_'+str(0.0)+'_high_res_T_'+str(T)+'_CFL/'
#fname_dir = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_2.0_high_res_CFL/'
fname_dir  = '../data/buttressing/water_depth_700_buttressing_25.0kPa_removed_50.0day/glacier_surf_slope_0.02_bed_slope_-0.02_flux_2.0_high_res_T_-20.0_CFL/'
fname_base = fname_dir + 'glacier_cliff_'
#fname_base = fname_dir + 'glacier_cliff_'


ice_thick = 800.0
Hab = 25.0
water_depth = ice_thick*910.0/1020 - Hab

ice_thick = 135.0
water_depth = 0.0
length= ice_thick*12

def bed_fun(x):
   b = -water_depth + bed_slope*(length-x)
   return b

#surf_slope =  0.01
#bed_slope =   -0.02
#flux = 4.0
#fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL3/glacier_cliff_'


L = []
Lb = []
t=[]
#u_term = []
h_cliff = []
bed = []
u_term = []
uz_term = []
for step in range(0,100000,10):

    fname_ext = str(step).zfill(3)+'.npz'
    fname = fname_base+fname_ext
    print(fname)


    # Load particle data
    try:
        print('reading file',fname)
        tracer_data = np.load(fname)
    except:
        print('Cannot read file')
        break

    #if tracer_data['t'].item()>0.2:
    #    break

    t.append(tracer_data['t'].item())


    # Load mesh
    temp_file = fname_dir + 'temp_'+str(step).zfill(3)+'.hdf'
    mesh_file = fname_base + str(step).zfill(3)+'.xml'
    mesh = Mesh(mesh_file)

    # Make particle class
    n=len(tracer_data['xm'])
    xp =  np.vstack((tracer_data['xm'],tracer_data['zm'])).transpose().reshape(n,2)
    pstrain = tracer_data['strain']
    pepsII = tracer_data['epsII']
    ptemp = tracer_data['temp']
    p = particles(xp, [pstrain,ptemp,pepsII], mesh)

    # Interpolate particles to mesh
    Vdg = FunctionSpace(mesh, 'DG',1)
    strain = Function(Vdg)
    lstsq_strain = l2projection(p, Vdg, 1) # First variable???
    try:
        lstsq_strain.project_mpm(strain) # Projection is stored in phih0
    except:
        file.read(strain, "strain")
    """
    Vcg = VectorFunctionSpace(mesh, 'CG',2)
    u = Function(Vcg)
    #ux,uz = u.split(deepcopy=True)
    file = HDF5File(MPI.comm_world,temp_file,'r')
    file.read(u, "u")
    """
    # Boundary mesh with portion above sea level marked
    bmesh = BoundaryMesh(mesh,'exterior',order=True)
    x = bmesh.coordinates()
    #ymax = np.max(x[x[:,0]==0,1])
    filter = (x[:,0]>1e-4) & (x[:,1]>0)
    xm = np.min(x[filter,0])
    id = np.argwhere(x[:,0]==xm).item()
    # Check if nodes increasing or decreasing
    if (x[id-1,0]>x[id,0]):
        # Decrease nodes
        inc = -1
        stop = 0
    else:
        # Increase nodes
        inc = 1
        stop = len(x)

    iold= id
    for i in range(id,stop,inc):
        #print(i,iold)


        if x[i,1]>-5.0:
            slope = (x[i,1]-x[iold,1])/(x[i,0]-x[iold,0])
            #print(i,iold,x[i,0],strain(x[i]),strain(x[i])>0.99,slope)
            #print(-slope*180/pi)
            if strain(x[i])>0.1:
                L.append(x[iold,0])
                #u_term.append(u(x[iold]))
                h_cliff.append(x[iold,1])
                #ux,uz = u(x[iold])
                #u_term.append(np.sqrt(ux**2))
                #uz_term.append(uz)
                #f = (x[:,1]<-10) & (strain(x[iold])<0.1) & (x[:,0]<L[-1])
                #idx=np.argmax((x[f,0]-L[-1]))
                #bed.append(x[f,1][idx])
                break
            elif np.abs(slope)>pi/7:
                L.append(x[iold,0])
                #u_term.append(u(x[iold]))
                h_cliff.append(x[iold,1])
                #ux,uz = u(x[iold])
                #u_term.append(np.sqrt(ux**2))
                #uz_term.append(uz)
                #f = (x[:,1]<-10) &  (x[:,0]<L[-1])
                #idx=np.argmax((x[f,0]-L[-1]))
                #bed.append(x[f,1][idx])
                break

            iold = i
    print('Terminus position',L[-1])
    filter = (x[:,0]>1e-4) & (x[:,1]<0)
    if sum(filter)>0:
        xm = np.min(x[filter,0])
        id = np.argwhere(x[:,0]==xm).item()
        # Check if nodes increasing or decreasing
        if (x[id-1,0]>x[id,0]):
            # Decrease nodes
            inc = -1
            stop = 0
        else:
            # Increase nodes
            inc = 1
            stop = len(x)

        iold= id
        for i in range(id,stop,inc):
            #print(i,iold)


            if x[i,1]<-5.0:
                slope = (x[i,1]-x[iold,1])/(x[i,0]-x[iold,0])
                #print(i,iold,x[i,0],strain(x[i]),strain(x[i])>0.99,slope)
                #print(-slope*180/pi)
                if strain(x[i])>0.1:
                    Lb.append(x[iold,0])
                    #u_term.append(u(x[iold]))
                    bed.append(x[iold,1])
                    #f = (x[:,1]<-10) & (strain(x[iold])<0.1) & (x[:,0]<L[-1])
                    #idx=np.argmax((x[f,0]-L[-1]))
                    #bed.append(x[f,1][idx])
                    break
                elif np.abs(slope)>pi/7:
                    Lb.append(x[iold,0])
                    #u_term.append(u(x[iold]))
                    bed.append(x[iold,1])
                    #f = (x[:,1]<-10) &  (x[:,0]<L[-1])
                    #idx=np.argmax((x[f,0]-L[-1]))
                    #bed.append(x[f,1][idx])
                    break

                iold = i
    #if i==stop-inc:
    #    print('No terminus identified')
    #    L.append(np.Nan)

fname = fname_base + 'term_pos.npz'

# Least squares fit to data
if len(t)>len(L):
    print('Not equal')
    t =t[0:-1]
t,L = np.array(t),np.array(L)
dLdt,b,rvalue,pvalue1,err=linregress(t,L)
print('Rate of terminus advance',dLdt,'Errror',err)
np.savez(fname, t=t, L=L,dLdt=dLdt, err=err,pvalue=pvalue1,rvalue=rvalue)

    #break
    # Maybe not needed?
    #pt = sort_boundary_nodes(bmesh)

    # Find top left node??

    # Make part above water?
