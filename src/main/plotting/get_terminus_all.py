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

import os.path
from os import path

def read_files(fname_base):
    L = []
    t=[]
    for step in range(0,1000000,10):

        fname_ext = str(step).zfill(3)+'.npz'
        fname = fname_base+fname_ext
        print(fname)


        # Load particle data
        if path.exists(fname):
            tracer_data = np.load(fname)
        else:
            print('Cannot read file')
            break

        t.append(tracer_data['t'].item())


        # Load mesh
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
        lstsq_strain.project_mpm(strain) # Projection is stored in phih0

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
            if x[i,1]>0.0:
                slope = (x[i,1]-x[iold,1])/(x[i,0]-x[iold,0])
                #print(x[i,0],strain(x[i]),strain(x[i])>0.99,slope,slope<-pi/3)
                #print(-slope*180/pi)
                if strain(x[i])>0.99:
                    L.append(x[iold,0])
                    break
                elif np.abs(slope)>pi/6:
                    L.append(x[iold,0])
                    break

                iold = i
        print('Terminus position',L[-1])
        if i==stop-inc:
            print('No terminus identified')
            L.append(np.Nan)


    fname = fname_base + 'term_pos.npz'
    print(fname)
    # Least squares fit to data
    #if len(t)>len(L):
    #    print('Not equal')
    #    t =t[0:-1]

    print(len(t),len(L))
    t = np.array(t)
    L = np.array(L)
    filter = (t>2/12) #
    t,L = np.array(t),np.array(L)
    dLdt,b,rvalue,pvalue1,err=linregress(t[filter],L[filter])
    print('Rate of terminus advance',dLdt,'Errror',err)
    np.savez(fname, t=t, L=L,dLdt=dLdt, err=err,pvalue=pvalue1,rvalue=rvalue)
    return None


"""
surf_slope = 0.02
bed_slopes = np.array([-0.02,-0.01,0.0,0.01,0.02,0.03,0.04])
fluxes = np.array([0.0,2.0,3.0,4.0,5.0,6.0])
for bed_slope in bed_slopes:
    for flux in fluxes:
        fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
        read_files(fname_base)
"""

surf_slope = 0.01
bed_slopes = np.array([-0.02,0.0,0.01])
fluxes = np.array([4.0,2.0,5.0])
for (bed_slope,flux) in zip(bed_slopes,fluxes):
    fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL3/glacier_cliff_'
    read_files(fname_base)


    #break
    # Maybe not needed?
    #pt = sort_boundary_nodes(bmesh)

    # Find top left node??

    # Make part above water?
