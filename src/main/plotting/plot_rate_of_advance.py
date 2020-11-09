"""
Calculate dLdt from t=0.25-1 for three fluxes
"""

import numpy as np
from scipy.stats import linregress
import pylab as plt
import matplotlib.colors as colors


tmin = 2./12
tmax = 10.0

surf_slope = 0.02
bed_slopes = np.array([-0.02,-0.01,0.0,0.01,0.02,0.03,0.04])
fluxes =  np.linspace(0,6,4)
fluxes = np.array([0.0,2.0,3.0,4.0,5.0,6.0])
dLdt = np.empty((len(fluxes),len(bed_slopes)))
err = np.empty((len(fluxes),len(bed_slopes)))

params = {
          'axes.labelsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
plt.rcParams.update(params)


def lsqr_dLdt(fname,tmin,tmax=100.0):
    data=np.load(fname)
    t = data['t']
    L = data['L']
    filter = (t>tmin) & (t<tmax)
    #filter = (t>2/12 )
    if sum(filter>0):
        #t=t[filter]
        #L=L[filter]
        dLdt,b,rvalue,pvalue1,err=linregress(t[filter],L[filter])
    else:
        filter = (t>tmin/2) #& (t<1.0)
        dLdt,b,rvalue,pvalue1,err=linregress(t[filter],L[filter])
    return dLdt,err,t[-1]




for i in range(len(fluxes)):
    flux = fluxes[i]
    for j in range(len(bed_slopes)):
        bed_slope = bed_slopes[j]
        fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
        fname = fname_base + 'term_pos.npz'
        dLdt[i,j],err[i,j],tm = lsqr_dLdt(fname,tmin,tmax)
        print('Bed slope',bed_slope,'Flux',flux,'Rate of advance (m/a)',dLdt[i,j],'err',err[i,j],'Max time',tm)

vmin = -12
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
height = width / 1.618/2*1.25
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.subplots_adjust(left=0.073,right=1.05,top=.75,bottom=0.25)

#fig.set_size_inches(width, height)
#plt.draw()
#fig, ax = plt.subplots(constrained_layout=True)
CS=plt.contour(bed_slopes-surf_slope,fluxes,dLdt/1e3, [-10,0], colors='k',linestyles='--')
c=plt.contourf(bed_slopes-surf_slope,fluxes,dLdt/1e3,v,cmap='bwr_r',norm=norm,extend='both',alpha=0.8)
plt.clabel(CS, inline=1, fontsize=10,fmt='%1.0f')

surf_slope = 0.01
bed_slopes = np.array([-0.02,0.0,0.01])
fluxes = np.array([4.0,2.0,5.0])
for (bed_slope,flux) in zip(bed_slopes,fluxes):
    fname ='../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL3/glacier_cliff_'+ 'term_pos.npz'
    dL1,err1,tm = lsqr_dLdt(fname,tmin,tmax)
    print('Bed slope',bed_slope,'Flux',flux,'Rate of advance (m/a)',dL1,'err',err1,'Max time',tm)
    plt.scatter([bed_slope-surf_slope],[flux],s=50,c=[dL1/1e3],cmap='bwr_r',norm=norm,marker='s',edgecolor='k')
surf_slope = 0.02
#plt.scatter([-0.01],[2],s=50,c=[-2137.5210508252835/1e3],cmap='bwr_r',norm=norm,marker='o',edgecolor='k')
#plt.scatter([0.0],[5],s=50,c=[487.56420334626046/1e3],cmap='bwr_r',norm=norm,marker='o',edgecolor='k')
cticks = [vmin, -9,-6,-3,0, vmax]
cbar=plt.colorbar(extend='both',ticks=cticks);#plt.clim([vmin,vmax])
cbar.ax.set_yticklabels(cticks)
cbar.set_label(r'$dL/dt$ (km/a)')
plt.text(-0.0+surf_slope+0.005-surf_slope,4.5,'Advance',fontsize=12,weight='bold')
plt.text(-0.005+surf_slope-surf_slope,2.8,'Retreat',fontsize=12,weight='bold')
plt.text(-0.0375,1.0,'Collapse',fontsize=12,weight='bold',rotation=90)

plt.xlabel(r'Thickness gradient, $dH/dx$',labelpad = 0)
plt.ylabel('Inflow velocity (km/a)')
ax = plt.gca()
secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.set_xlabel('Bed slope')
#plt.annotate ('', (-0.04, 7.5), (-0.02, 7.5), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
#plt.text(-0.03-0.0075, 8, r'retrograde bed slope')
#plt.annotate ('', (-0.02, 7.5), (0.02, 7.5), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
#plt.text(0.0-0.006, 8, r'prograde bed slope')


plt.annotate ('', (-0.04, -2.15), (-0.0, -2.15), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
plt.text(-0.025-0.0075, -2.95, r'H increases upstream',fontsize=10)
plt.annotate ('', (0.0, -2.15), (0.02, -2.15), arrowprops={'arrowstyle':'<->','linewidth':2},annotation_clip=False)
plt.text(0.001, -2.95, r'H decreases upstream',fontsize=10,wrap=True)

#plt.plot(-0.012,2.3,marker='*',color='k',markersize=10)
#ax.annotate('Thwaites', xy=(-0.01, 2), xytext=(-0.018,1.85),#
#            arrowprops=dict(facecolor='black', shrink=0.015))
plt.text(-0.018,1.85,'Thwaites')
plt.show()
plt.savefig('Figures/Figure_rate_of_advance.pdf',bbox_inches='tight')


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
