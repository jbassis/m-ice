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

def read_file(fname_base,step):
        fname_ext = str(step).zfill(3)+'.npz'
        fname = fname_base+fname_ext
        print(fname)


        # Load particle data
        if path.exists(fname):
            tracer_data = np.load(fname)
        else:
            print('Cannot read file')


        t=tracer_data['t'].item()


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
                if strain(x[i])>0.1:
                    L=x[iold,0]
                    break
                elif np.abs(slope)>pi/6:
                    L=x[iold,0]
                    break

                iold = i
        print('Terminus position',L)
        # Extract profile centered on terminus
        filter = (x[:,1]>0.0) & (x[:,0]>10) & (x[:,0]<L+5*800)
        xp = x[filter,0]-L
        zp = x[filter,1]
        idx = np.argsort(xp)
        xp = xp[idx]
        zp = zp[idx]

        #fname = fname_base + 'term_pos.npz'
        #print(fname)

        return xp,zp


#"""
surf_slope = 0.02
bed_slope =  0.0
flux = 0.0
fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
fname_base = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.03_flux_2.0_high_res_T_-10.0_CFL/glacier_cliff_'
fname_base = '../data/cliff/water_depth_700_no_failure/glacier_surf_slope_0.02_bed_slope_0.04_flux_3.0_high_res_T_-20.0_CFL/glacier_cliff_'
fname_base = '../data/cliff/water_depth_700_super_buoyancy_crevasse/glacier_surf_slope_0.02_bed_slope_0.02_flux_3.0_high_res_T_-20.0_CFL/glacier_cliff_'
#fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing50.0kPa_CFL/glacier_cliff_'
#for step in range(720,800,20):
for step in range(0,4000,10):

    xp,zp=read_file(fname_base,step)
    plt.clf()

    plt.plot(-xp,zp)
    plt.ylim([0,200])
    plt.pause(1e-16)
