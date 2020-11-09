"""
Calculate dLdt from t=0.25-1 for three fluxes
"""

import numpy as np
from scipy.stats import linregress
import pylab as plt
import matplotlib.colors as colors



length = 9600.0
ice_thick = 800.0
surf_slope = 0.02
bed_slopes = np.array([-0.02,-0.01,0.0,0.01,0.02,0.03,0.04])
fluxes = np.array([0.0,2.0,3.0,4.0,5.0,6.0])
dLdt = []
err = []
V = []
S = []
F = []

for i in range(len(fluxes)):
    flux = fluxes[i]

    for j in range(len(bed_slopes)):
        bed_slope = bed_slopes[j]
        S.append(bed_slope-surf_slope)
        V.append(flux)
        flux_in = ((bed_slope-surf_slope)*length + ice_thick)*flux
        F.append(flux_in)
        fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
        fname = fname_base + 'term_pos.npz'
        data=np.load(fname)
        t = data['t']
        L = data['L']
        filter = (t>2/12) #& (t<1.0)
        #filter = (t>2/12 )
        if sum(filter>0):
            t=t[filter]
            L=L[filter]
            dL,b,rvalue,pvalue1,e=linregress(t,L)

            print('Bed slope',bed_slope,'Flux',flux,'Rate of advance (m/a)',dLdt[-1],'err',err[-1],'Maximum time',np.max(t))
        else:
            print('Time series not long enough for flux',flux,'bed slope',bed_slope)
            dL,b,rvalue,pvalue1,e=linregress(t,L)
        dLdt.append(dL)
        err.append(e)

vmin = -4
vmax =  3
v = np.linspace(vmin, vmax, 11, endpoint=True)
norm = colors.DivergingNorm(vmin=vmin, vcenter=0, vmax=vmax)
#clf();c=plt.pcolormesh(slope-0.01/2,flux-1/2,dLdt/1e3,cmap='bwr_r',norm=norm);plt.colorbar(extend='both');plt.clim([vmin,vmax])

def forward(x):
    return x+surf_slope#-surf_slope


def inverse(x):
    return x-surf_slope
#slope = np.array([-0.02,-0.01,0.0,0.01,0.02,0.03,0.04])
#flux = np.array([0,2,4,6])

plt.close('all')
width  = 3.487*2
height = width / 1.618/2
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.subplots_adjust(left=0.073,right=1.05,top=.725,bottom=0.2)

#fig.set_size_inches(width, height)
#plt.draw()
#fig, ax = plt.subplots(constrained_layout=True)
dLdt = np.array(dLdt)
F = np.array(F)
plt.tricontour(S,F/1e3,dLdt/1e3, [0], colors='k',linestyles='--')
c=plt.tricontourf(S,F/1e3,dLdt/1e3,v,cmap='bwr_r',norm=norm,extend='both');
#plt.scatter([-0.03],[4],s=50,c=[351.5094769326309 /1e3],cmap='bwr_r',norm=norm,marker='o',edgecolor='k')
#plt.scatter([-0.01],[2],s=50,c=[-2137.5210508252835/1e3],cmap='bwr_r',norm=norm,marker='o',edgecolor='k')
#plt.scatter([0.0],[5],s=50,c=[487.56420334626046/1e3],cmap='bwr_r',norm=norm,marker='o',edgecolor='k')

cbar=plt.colorbar(extend='both',ticks=[vmin,0,vmax]);#plt.clim([vmin,vmax])
cbar.ax.set_yticklabels([vmin, 0,vmax])
cbar.set_label(r'$dL/dt$ (km/a)')
#plt.text(-0.0+surf_slope+0.005-surf_slope,4.5,'Advance',fontsize=12,weight='bold')
#plt.text(-0.025+surf_slope-surf_slope,2.8,'Retreat',fontsize=12,weight='bold')
#plt.text(-0.0375,1.0,'Collapse',fontsize=12,weight='bold',rotation=90)

plt.xlabel(r'Initial thickness gradient',labelpad = 0)
plt.ylabel(r'Flux ($10^6$ $m^2$/a)')
ax = plt.gca()
secax = ax.secondary_xaxis('top', functions=(forward, inverse))
#secax.set_xlabel('Thickness gradient')
#plt.annotate ('', (-0.04, 7.5), (-0.02, 7.5), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
#plt.text(-0.03-0.0075, 8, r'retrograde bed slope')
#plt.annotate ('', (-0.02, 7.5), (0.02, 7.5), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
#plt.text(0.0-0.006, 8, r'prograde bed slope')

plt.show()


#plt.annotate("", xy=(-0.04, 8), xytext=(0., 8),
#             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"))
#plt.text(-0.02,4,'Neutral',fontsize=12,weight='bold')
"""
plt.errorbar(bed_slopes-surf_slope,dLdt[0,:],yerr=err[0,:],fmt='-o',capsize=5,label='0 km/a')
plt.errorbar(bed_slopes-surf_slope,dLdt[1,:],yerr=err[1,:],fmt='-<',capsize=5,label='2 km/a')
plt.errorbar(bed_slopes-surf_slope,dLdt[2,:],yerr=err[1,:],fmt='-s',capsize=5,label='4 km/a')
plt.errorbar(bed_slopes-surf_slope,dLdt[3,:],yerr=err[1,:],fmt='-d',capsize=5,label='6 km/a')
plt.ylim([-11e3,4e3])
plt.xlabel('Ice thickness gradient')
plt.ylabel('dL/dt (m/a)')
"""



        #break
