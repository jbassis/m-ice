"""
Pull out files from t=1 a
"""

import pylab as plt
from numpy import *
import numpy as np
from dolfin import *



surf_slope = 0.02
bed_slope =  0.0
flux = 2.0

step = 4050

fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
fname_ext = str(step).zfill(3)+'.npz'

fname = fname_base+fname_ext
print(fname)


# Load data
data = load(fname)
t=data['t']
print('Time',t)

yr = int(mod(t,365))
#day = round((t - yr)*365,1)
day = (t - yr)*365
hr,day = modf(day)
hr = round(hr*24,0)
title_str = 'Time: '+str(yr).zfill(1)+ 'a '+str(int(day)).zfill(1)+'d '+str(int(hr)).zfill(1)+'hr'
print(title_str)

# Load mesh
mesh_file = fname_base + str(step).zfill(3)+'.xml'
mesh_data = Mesh(mesh_file)

plt.clf()
plot(mesh_data)
