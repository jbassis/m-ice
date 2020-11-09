import numpy as np
import pylab as plt

from importlib import reload


import matplotlib
matplotlib.style.use('seaborn-colorblind')

import matplotlib.gridspec as gridspec


import plot_cliff
reload(plot_cliff)
from plot_cliff import *
plt.close('all')

color = '#CC6677'
color2 = '#88CCEE'

width  = 3.5*2.2
height = width / 2.5/2.5
fig=plt.figure(num=2, figsize=(width, height), facecolor='w', edgecolor='k')
fig.set_size_inches(width,height,forward=True)
#ax10 = fig.add_axes([0.125, 0.2, 0.75, 0.75])

#ax10 =  fig.add_axes([snap_start+0.0125, snap_start_vert+0.3, snap_width*3-0.04-0.0175, 0.125])
#ax10.patch.set_alpha(0)

L = 1620
fname1 = '../data/cliff/water_depth_0/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_CFL/glacier_cliff_'+ 'term_pos.npz'
fname2 = '../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_CFL/glacier_cliff_'+ 'term_pos.npz'
fname3 = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_CFL/glacier_cliff_'+ 'term_pos.npz'
fname4 = '../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
fname5 = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
fname6 = '../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_no_buttressing'+ 'term_pos.npz'
fname7 = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_no_buttressing'+ 'term_pos.npz'

gs = fig.add_gridspec(1,6)
plt.subplots_adjust(wspace = 3.0,hspace=0.5)
ax10 = fig.add_subplot(gs[0,0:2])

#ax10 = plt.subplot(3,1,1)
#plot_term_buttressing(ax10,fname1=fname1,fname2=fname2,fname3=fname3,fname4=fname4,fname5=fname5,tmax=250)
plot_term_buttressing(ax10,fname1=fname1,fname2=None,fname3=None,fname4=None,fname5=None,tmax=180)

#data=np.load(fname6)
#ax10.plot(data['t']*365,(data['L']-4800)/1e3,'-.',color='#009E73')
#data=np.load(fname7)
#ax10.plot(data['t']*365,(data['L']-9600)/1e3,'-.',color='#D55E00')
plt.xticks([0,90,180],['day 0', 'day 90', 'day 180'])
#plt.yticks([-2,0.25])
ax10.set_ylim([-2,0.25])
ax10.set_yticks([-2,0.25])
plt.text(30,0.325,'135 m',fontsize=10)
#plt.text(95,-1.25,'buttressing off',fontsize=10,wrap=True,color=color)
#ax10.axvline(90,color='gray',ls='--')


#ax10.axvline(90,color='gray',ls='--')

#plt.text(30,0.1,'buttressing',fontsize=10,color='#D55E00')
#plt.text(100,0.1,'no buttressing',fontsize=10,color='#D55E00')

#ax10.set_xticks([0,100,200])
#ax10.set_xlim([0,200])
#ax10.set_xticklabels(['day 0', 'day 100', 'day 200'])

#plt.ylim([-2,0.25])
#plt.text(200,0.0,u'135 m',color='#0072B2',fontsize=10,fontweight='bold')
#plt.text(24,-1.9,u'400 m',color='#009E73',fontsize=10,fontweight='bold')
#plt.text(120,-0.7,u'800 m',color='#D55E00',fontsize=10,fontweight='bold')

#plt.text(200,-0.5,u'135 m',color='#0072B2',fontsize=10,fontweight='bold')
#plt.text(200,-0.875,u'400 m',color='#009E73',fontsize=10,fontweight='bold')
#plt.text(200,-1.25,u'800 m',color='#D55E00',fontsize=10,fontweight='bold')

plt.text(0,0.325,'A',fontsize=10,fontweight='bold')

#ax11 = plt.subplot(3,1,2)
ax11 = fig.add_subplot(gs[0,2:4])
data5=np.load(fname2)
ax11.plot(data5['t']*365,(data5['L']-4800)/1e3,'-',color=color,linewidth=2)
plot_term_buttressing(ax11,fname1=None,fname2=None,fname3=None,fname4=fname4,fname5=None,tmax=180)
data=np.load(fname6)
ax11.plot(data['t']*365,(data['L']-4800)/1e3,'--',color=color)
plt.xticks([0,90,180],['day 0', 'day 90', 'day 180'])
plt.yticks([-2,0.25])
#plt.text(36,-0.7,'25 kPa',color='#009E73',fontsize=10,fontweight='bold')
#plt.text(24,-1.9,'0 kPa',color='#009E73',fontsize=10,fontweight='bold')
#plt.text(35,-1.9,'25 kPa',fontsize=10)
#plt.text(100,-0.5,'unbuttressed',fontsize=10)

#plt.text(10,-0.12,'buttressed',fontsize=10)
#plt.text(95, 0.0,'buttressing off',fontsize=10,wrap=True,color=color)
#plt.text(95,-0.5,'buttressing on',fontsize=10,wrap=True,color='#0072B2')


ax11.axvline(90,color='gray',ls='--')
plt.text(0,0.325,'B',fontsize=10,fontweight='bold')
plt.text(30,0.325,'400 m',fontsize=10)
ax11.spines['left'].set_visible(False)
ax11.yaxis.set_visible(False)
#plt.text(-5,0.05,'buttressing',fontsize=10,color='#0072B2')

#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#mark_inset(ax10, ax11, loc1=2, loc2=4, fc="none", ec="0.5")
#ax11.yaxis.set_visible(False)
#ax11.spines['left'].set_visible(False)


#ax12 = plt.subplot(3,1,3)
ax12 = fig.add_subplot(gs[0,4:6])
data5=np.load(fname3)
ax12.plot(data5['t']*365,(data5['L']-9600)/1e3,'-',color=color,linewidth=2)
plot_term_buttressing(ax12,fname1=None,fname2=None,fname3=None,fname4=None,fname5=fname5,tmax=180)

data=np.load(fname7)
#ax12.plot(data['t']*365,(data['L']-9600)/1e3,'-.',color='#D55E00')
ax12.plot(data['t']*365,(data['L']-9600)/1e3,'--',color=color)

plt.xticks([0,90,180],['day 0', 'day 90', 'day 180'])
plt.ylim([-2.0,0.25])
plt.yticks([-2,0.25])
ax12.spines['left'].set_visible(False)
ax12.yaxis.set_visible(False)
ax12.axvline(90,color='gray',ls='--')
#plt.text(35,-1.9,'25 kPa',fontsize=10)
#plt.text(95,-1.9,'buttressing off',fontsize=10,color=color)

plt.text(-150,-1.05+0.85,'buttressing on',fontsize=10,color=color2,fontweight='bold')
plt.text(-150,-1.5+0.85,'buttressing off',fontsize=10,color=color,fontweight='bold')

#plt.text(-10,-1.9,'buttressing',fontsize=10,color='#0072B2')

plt.text(0,0.325,'C',fontsize=10,fontweight='bold')
plt.text(30,0.325,'800 m',fontsize=10)
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#mark_inset(ax10, ax12, loc1=2, loc2=4, fc="none", ec="0.5")

#plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(wspace=1.5)


plt.savefig('term_pos.pdf',dpi=800)
plt.show()
