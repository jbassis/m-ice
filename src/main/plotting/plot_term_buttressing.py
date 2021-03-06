"""
Plot effect of temperature on terminus position

"""

import numpy as np
import pylab as plt

import matplotlib
matplotlib.style.use('seaborn-colorblind')


L = 4800.0
fname1 ='../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing50.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
data2=np.load(fname1)

fname2 = '../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
data3=np.load(fname2)

fname3 = '../data/cliff/water_depth_290/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_CFL/glacier_cliff_'+ 'term_pos.npz'
data4=np.load(fname3)

#"""
L = 9600.0
fname1 ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing50.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
data2=np.load(fname1)

fname2 = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_'+ 'term_pos.npz'
data3=np.load(fname2)

fname3 = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_CFL/glacier_cliff_'+ 'term_pos.npz'
data4=np.load(fname3)
#"""



plt.close('all')
width  = 3.5
height = width / 2.5
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.set_size_inches(width,height,forward=True)
ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])
#fig.subplots_adjust(left=0.073,right=1.05,top=.725,bottom=0.2)

#p1=ax.plot(data1['t'],(data1['L']-L)/1e3,'-.')
#p2=ax.plot(data2['t'],(data2['L']-L)/1e3,'-.',linewidth=2,color='crimson')
p3=ax.plot(data3['t'],(data3['L']-L)/1e3,'--',color='gray',linewidth=2)
p4=ax.plot(data4['t'],(data4['L']-L)/1e3,'-',color='slateblue',linewidth=2)
plt.xlabel('Time (a)')
plt.ylabel(r'$\Delta L$ (km)')
plt.ylim([-2.,0.5])

plt.xlim([0.0,0.2])

#plt.text(0.155,-1.25-0.2,u'0 kPa',color=p2[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.155,-1.6-0.35,u'25 kPa',color=p3[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.155,-1.95-0.5,u'50 kPa',color=p4[0].get_color(),fontsize=10,fontweight='bold')

#if L == 9600.0:
#plt.text(0.0025,-1.9,u'50 kPa',color=p2[0].get_color(),fontsize=10,fontweight='bold')
plt.text(0.0025,-1.7,u'25 kPa',color=p3[0].get_color(),fontsize=10,fontweight='bold')
plt.text(0.0025,-1.2,u'0 kPa',color=p4[0].get_color(),fontsize=10,fontweight='bold')
#else:
    #plt.text(0.1575,-0.2,u'50 kPa',color=p2[0].get_color(),fontsize=10,fontweight='bold')
    #plt.text(0.09,-1.35,u'25 kPa',color=p3[0].get_color(),fontsize=10,fontweight='bold')
    #plt.text(0.0225,-1.9,u'0 kPa',color=p4[0].get_color(),fontsize=10,fontweight='bold')




#plt.text(0.001,-0.9,' 0$^\circ$C',color=p1[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.001,-1.4,'-5$^\circ$C',color=p2[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.001,-1.9,'-20$^\circ$C',color=p4[0].get_color(),fontsize=10,fontweight='bold')

#plt.text(0.01,0.8,u'0\u00B0C',color=p1[0].get_color(),fontsize=10)
#plt.text(0.1,.6,u'-5\u00B0C',color=p2[0].get_color(),fontsize=10)
#plt.text(0.14,-2.25,u'-20\u00B0C',color=p4[0].get_color(),fontsize=10)

plt.show()
#ax.legend()
