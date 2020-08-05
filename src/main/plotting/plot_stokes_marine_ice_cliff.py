"""
Plot 2D rift profiles
"""
#%matplotlib inline
import pylab as plt
from numpy import *
import numpy as np

time_factor = 86400.0*365.24

# Create figure
width  = 3.487*3.0/1.25
height = width / 1.618
height = width/2.25

# Base directory for file
#fname_base= 'data/tmp13/glacier_cliff_'
#fname_base= 'data/prograde20_mid_res/glacier_cliff_'
#fname_base= 'data/retrograde6_mid_res/glacier_cliff_'
#fname_base= 'data/retrograde_test_old2/glacier_cliff_'
fname_base= 'data/prograde/prograde13_low_res/glacier_cliff_'
fname_base = 'data/prograde/cold_ice/fixed_flux/inflow_6kma_low_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_low_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_high_res_no_slope/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_low_res_old/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_high_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_300_high_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_high_res_retro/glacier_cliff'
fname_out = fname_base#+'junk'

#fname_out = 'data/tmp14/glacier_cliff_'

#12: 791
# Specify file to load
step = 1
nmax = 4650

# Geometric variables
Hab = 25.0
ice_thick = 800.0

ice_thick = 400.0
Hab = 65.0


length= ice_thick*12
water_depth = ice_thick*910.0/1020 - Hab

# Define surface, bottom, etc
if fname_base[5:10] == 'retro':
    surf_slope = 0.02
    bed_slope = -0.01
else:
    surf_slope = 0.02
    bed_slope =  0.01

surf_slope= 0.02
bed_slope = 0.0

notch_height = 200.0*0
notch_length = 40.0
notch_slope = notch_height/notch_length/2
L = length - 60.0
bump_width = ice_thick
bump_height = 0.1*0

def bot_fun(x):
   b = -water_depth + bed_slope*(length-x)
   return b

def bed_fun(x):
   b = -water_depth + bed_slope*(length-x)
   return b

def ice_thick_fun(x):
    return np.maximum(-bot_fun(x)*1020.0/910.0,ice_thick)

def surf_fun(x):
   s = -water_depth + ice_thick + surf_slope*(L-x)#-  notch_slope*(L+notch_length-x)*(x<(L+notch_length))*(x>(L)) + notch_slope*(L-x-notch_length)*(x>(L-notch_length))*(x<=(L))
   return s

xs=linspace(1e3,length,11)
hs=surf_fun(xs)
bed = -793.7594315851419
max_length = 12e3
xx=linspace(1e3,max_length,1001)
plt.close('all')
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.subplots_adjust(left=-0.05, bottom=.07, right=0.96, top=.95)
#plt.tick_params(axis='both', which='major', labelsize=20)

#for step in xrange(8350,9040,10):
for step in xrange(0,2770,10):

   #fname_base= 'data/tmp8/restart_glacier_cliff_'

   fig.clf()
   # File extension
   fname_ext = str(step).zfill(3)+'.npz'


   # Filename
   fname = fname_base+fname_ext
   print fname


   # Load data
   data = load(fname)




   # Triangulation and point coordinates
   x=data['xm'];z=data['zm'];eta=data['eta'];eta_visc=data['eta_visc_m']
   filter =(x< max_length) & (eta>0.0) & (eta_visc>0.0) & (x>1e3)
   x = x[filter];z=z[filter]
   strain = maximum(data['strain'],1e-16)[filter]
   epsII = data['epsII'][filter]
   t=data['t']
   eta = data['eta'][filter]
   eta_visc = data['eta_visc_m'][filter]
   D = np.minimum(np.maximum(1-eta/eta_visc,0.0),1.0)
   #xx=data['xx']
   #hsurf=data['hsurf']
   #hbot=data['hbot']
   #width=data['width']


   # Make plot
   pts_size = 1.0
   ax2=plt.subplot(3,1,1)
   c=ax2.scatter(x,z,s=pts_size,c=np.log10(strain),vmin=-4,vmax=1)
   ax2.fill_between([1e3,max_length],0.0,-800,color='dodgerblue',alpha=0.33,zorder=0)
   ax2.fill_between(xx,bed_fun(xx),-800,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   ax2.plot(xs,hs,'--',color='Gray',linewidth=2)
   ax2.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   #ax2.plot([xs[0],xs[-1]],[hs[-1],hs[-1]],'-.',color='Orange',linewidth=2)
   #plt.ylabel('Distance (m)',fontsize=20)
   yr = int(mod(t,365))
   day = round((t - yr)*365,1)
   #title_str = 'Time: '+str(yr).zfill(0)+ 'a ' + str(day).zfill(0) ,'days'
   #title_str = 'Time: '+str(round(t*365.24,2)).zfill(4)+ ' days'
   title_str = 'Time: '+str(yr).zfill(1)+ 'a '+str(day).zfill(1)+'d'
   ax2.text(4e3,300,title_str,fontsize=14)
   #plt.title(title_str,fontsize=20)
   plt.xlim([1e3,max_length])
   #plt.ylim([-850,150])
   plt.axis('equal')

   # Add colorbar
   cbaxes = fig.add_axes([0.9, 0.695, 0.015, 0.255])
   cbar2 = plt.colorbar(c, cax = cbaxes,extend='both')
   cbar2.set_ticks([-4,1])
   cbar2.set_label(r'$\epsilon_p$', fontsize=14)
   cbar2.set_ticklabels([r'$10^{-4}$',r'$10^{-1}$'])

   # Add plot
   ax3=plt.subplot(3,1,2)
   c=ax3.scatter(x,z,s=pts_size,c=np.log10(epsII),vmin=-8,vmax=-5)
   ax3.fill_between([1e3,max_length],0.0,-800,color='dodgerblue',alpha=0.33,zorder=0)
   ax3.fill_between(xx,bed_fun(xx),-800,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   ax3.plot(xs,hs,'--',color='Gray',linewidth=2)
   ax3.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   #ax3.plot([xs[0],xs[-1]],[hs[-1],hs[-1]],'-.',color='Orange',linewidth=2)
   #plt.ylabel('Distance (m)',fontsize=20)
   plt.xlim([1e3,max_length])
   #plt.ylim([-850,150])
   plt.axis('equal')

   # Add colorbar
   cbaxes = fig.add_axes([0.9, 0.695-0.315, 0.015, 0.255])
   cbar3 = plt.colorbar(c, cax = cbaxes,extend='both')
   cbar3.set_ticks([-8,-5])
   cbar3.set_ticklabels([r'$10^{-8}$',r'$10^{-5}$'])
   cbar3.set_label('$\dot \epsilon_{II}$ (s$^{-1}$)',fontsize=14)

   # Add plot
   ax1=plt.subplot(3,1,3)
   c3=ax1.scatter(x,z,s=pts_size,c=D,vmin=0,vmax=1)
   ax1.fill_between([1e3,max_length],0.0,-800,color='dodgerblue',alpha=0.33,zorder=0)
   ax1.fill_between(xx,bed_fun(xx),-800,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   ax1.plot(xs,hs,'--',color='Gray',linewidth=2)
   ax1.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   #ax1.plot([xs[0],xs[-1]],[hs[-1],hs[-1]],'-.',color='Orange',linewidth=2)
   plt.xlim([1e3,max_length])
   #plt.ylim([-850,150])
   plt.axis('equal')

   # Add colorbar
   cbaxes = fig.add_axes([0.9, 0.695-0.315*2, 0.015, 0.255])
   cbar3 = plt.colorbar(c3, cax = cbaxes)
   cbar3.set_ticks([0,1])
   cbar3.set_label('Enhancement', fontsize=14)
   """
   # Subplot 1: Plastic strain
   ax1=plt.subplot(3,1,3)

   plt.plot(xx,hsurf,'k',linewidth=2)
   plt.plot(xx,hbot,'k',linewidth=2)
   plt.fill_between(xx,hsurf,hbot,color='gray')
   plt.ylabel('Elevation (m)',fontsize=20)
   plt.xlim([0,max_length])
   plt.ylim([-400,50])
   ax1.set_xticks([0,max_length])
   ax1.set_yticks([-400,50])
   ax1.tick_params(labelsize=14)
   ax1.spines['top'].set_visible(False)
   #ax.spines['left'].set_visible(False)
   ax1.spines['bottom'].set_visible(True)
   """


   for ax in [ax1,ax2,ax3]:
       ax.set_xticks([])
       ax.set_yticks([])
       ax.spines['right'].set_visible(False)
       ax.spines['top'].set_visible(False)
       #ax.spines['left'].set_visible(False)
       #ax.spines['left'].set_linewidth(2)
       ax.spines['bottom'].set_visible(False)
       ax.spines['left'].set_visible(False)
       tick_xx = [1e3-150,1e3-150];tick_yy = [-800,200]
       ax.plot(tick_xx,tick_yy,'_',color='k',clip_on=False,linewidth=2)
       ax.vlines(1e3-200, -800,200, colors='k',linewidth=2,clip_on=False)
       ax.text(1e3-500,-200,'1 km',fontsize=14,rotation=90)
       ax.set_xlim([1e3,max_length])


       #ax.set_yticklabels([0,400])
       #ax.tick_params(axis='both', which='major', labelsize=18)
       #ax.vlines(x=-200,ymin=bed,ymax=bed+900,linewidth=1,color='k',zorder=100,clip_on=False)
       #ax.xaxis.set_tick_params(width=1)
       #ax.yaxis.set_tick_params(width=1)
       #ax.set_ylim([-850,150])
       #ax.set_xticks([0,max_length])
       #ax.set_yticks([-380,80])
       #ax.tick_params(labelsize=14)
   #ax1.spines['bottom'].set_visible(True)
   #ax1.tick_params(labelsize=14)
   #ax1.set_xticks([0,max_length])
   #ax1.spines['bottom'].set_linewidth(2)

   ax1.hlines(-1200, 1e3,max_length, colors='k',linewidth=2,clip_on=False)
   ax2.hlines(-1200, 1e3,max_length, colors='w',linewidth=2,clip_on=False)
   ax3.hlines(-1200, 1e3,max_length, colors='w',linewidth=2,clip_on=False)
   tick_x = [1e3,max_length];tick_y = [-1200+60,-1200+60]
   ax1.plot(tick_x,tick_y,'|',color='k',clip_on=False,linewidth=4)
   #ax1.spines['bottom'].set_visible(True)
   #ax1.spines['bottom'].set_lw(1)
   #ax1.tick_params(axis="x",direction="in")
   #ax4.hlines(y=bed-200,xmin=0,xmax=max_length,linewidth=2,color='k',zorder=100,clip_on=False)
   #ax1.set_xticks([0,max_length])

   #tick_xx = [-45,-45];tick_yy = [-360,60]
   #ax1.plot(tick_xx,tick_yy,'_',color='k',clip_on=False,linewidth=2)
   #plt.xlabel('Distance (m)',fontsize=20)
   max_length_txt = str(int((max_length - 1e3)/1e3))+' km'
   ax1.text(max_length/2+1e3,-1600,max_length_txt,fontsize=14)
   #ax1.text(max_length/2,-550,'10 km',fontsize=14)
   plt.pause(1e-16)
   plt.show()
   fname_ext = str(step).zfill(3)+'.png'
   #fname_base= 'data/tmp8/glacier_cliff_'
   fname = fname_out+fname_ext
   plt.savefig(fname)
