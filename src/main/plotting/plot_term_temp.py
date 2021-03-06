"""
Plot effect of temperature on terminus position

"""

import numpy as np
import pylab as plt

import matplotlib
matplotlib.style.use('seaborn-colorblind')


L = 9600.0

fname1 ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_3.0_high_res_T_-0_CFL/glacier_cliff_'+ 'term_pos.npz'
data1=np.load(fname1)
t = data1['t']
yr = np.array((np.mod(t,365)),dtype=int)
day1 = (t - yr)*365


fname2 ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_3.0_high_res_T_-5.0_CFL/glacier_cliff_'+ 'term_pos.npz'
data2=np.load(fname2)
t = data2['t']
yr = np.array((np.mod(t,365)),dtype=int)
day2 = (t - yr)*365

fname3 ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_3.0_high_res_T_-15.0_CFL/glacier_cliff_'+ 'term_pos.npz'
data3=np.load(fname3)
t = data3['t']
yr = np.array((np.mod(t,365)),dtype=int)
day3 = (t - yr)*365

fname4 ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_3.0_high_res_CFL/glacier_cliff_'+ 'term_pos.npz'
data4=np.load(fname4)
t = data4['t']
yr = np.array((np.mod(t,365)),dtype=int)
day4 = (t - yr)*365

#fname5 ='../data/cliff/water_depth_700_high_friction/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_CFL/glacier_cliff_'+ 'term_pos.npz'
#data5=np.load(fname5)
#t = data5['t']
#yr = np.array((np.mod(t,365)),dtype=int)
#day5 = (t - yr)*365


plt.close('all')
width  = 3.5
height = width / 2.5
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.set_size_inches(width,height,forward=True)
ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])
#fig.subplots_adjust(left=0.073,right=1.05,top=.725,bottom=0.2)

#p1=ax.plot(data1['t'],(data1['L']-L)/1e3,'-.')
p4=ax.plot(day4,(data4['L']-L)/1e3,'-',color='slateblue',linewidth=2)
p3=ax.plot(day3,(data3['L']-L)/1e3,'--',color='gray',linewidth=2)

p2=ax.plot(day2,(data2['L']-L)/1e3,'-.',linewidth=2,color='crimson')
plt.xlabel('Time (day)')
plt.ylabel(r'$\Delta L$ (km)')
plt.ylim([-2,1])
plt.ylim([-0.5,0.5])

plt.xlim([0.0,0.2*365])

plt.text(0.13*365,0.3,u'-5\u00B0C',color=p2[0].get_color(),fontsize=10,fontweight='bold')
plt.text((0.05+0.11)*365,0.1,u'-15\u00B0C',color=p3[0].get_color(),fontsize=10,fontweight='bold')
plt.text(0.125*365,-0.42,u'-20\u00B0C',color=p4[0].get_color(),fontsize=10,fontweight='bold')

#plt.text(0.001,-0.9,' 0$^\circ$C',color=p1[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.001,-1.4,'-5$^\circ$C',color=p2[0].get_color(),fontsize=10,fontweight='bold')
#plt.text(0.001,-1.9,'-20$^\circ$C',color=p4[0].get_color(),fontsize=10,fontweight='bold')

#plt.text(0.01,0.8,u'0\u00B0C',color=p1[0].get_color(),fontsize=10)
#plt.text(0.1,.6,u'-5\u00B0C',color=p2[0].get_color(),fontsize=10)
#plt.text(0.14,-2.25,u'-20\u00B0C',color=p4[0].get_color(),fontsize=10)

plt.show()
plt.savefig('Figures/term_temp.pdf')
#ax.legend()
