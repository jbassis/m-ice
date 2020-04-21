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
height = width/7*2

# Base directory for file
#fname_base= 'data/tmp13/glacier_cliff_'
#fname_base= 'data/prograde20_mid_res/glacier_cliff_'
#fname_base= 'data/retrograde6_mid_res/glacier_cliff_'
#fname_base= 'data/retrograde_test_old2/glacier_cliff_'

fname_base = 'data/cliff/water_depth_688_low_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_0_high_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_0_high_res_slope/glacier_cliff_'
fname_base = 'data/water_depth_0_high_res_low_visc/glacier_cliff_'
#fname_base = 'data/cliff/water_depth_300_high_res/glacier_cliff_'
fname_base = 'data/cliff/water_depth_0_high_res_slope_new/glacier_cliff_'
#fname_base = 'data/cliff/water_depth_688_high_res_new6/glacier_cliff_'


fname_base = 'data/cliff/water_depth_300_high_res_retro_steep_bed/glacier_cliff_'
fname_base = 'data/cliff/water_depth_688_high_res_retro/glacier_cliff_'
fname_base = 'data/cliff/15C/water_depth_306_high_res_flux_1.5_slope-0.01_min_visc/glacier_cliff_'
fname_base = 'data/cliff/15C/water_depth_306_high_res_flux_1.5_slope-0.01/glacier_cliff_'
fname_base = 'data/cliff/15C/water_depth_306_high_res_flux_2.0_slope-0.00_min_visc/glacier_cliff_'
fname_base ='data/cliff/15C/water_depth_306_low_res_flux_0.0_slope-0.01/glacier_cliff_'
fname_base ='data/cliff/15C/water_depth_306_low_res_flux_0.0_slope-0.0/glacier_cliff_'
fname_base = 'data/cliff/15C/water_depth_306_high_res_flux_0.0_slope-0.0/glacier_cliff_'
fname_base = 'data/cliff/flat/water_depth_306_high_res_flux_1.0_slope--0.01/glacier_cliff_'
fname_base ='data/cliff/flat/water_depth_306_high_res_flux_2.0_slope--0.01/glacier_cliff_'
#fname_base = 'data/cliff/water_depth_700/glacier_min_yield_20_surf_slope_0.02_med_visc_bed_slope_0.01_flux_2.0/glacier_cliff_'
fname_base = 'data/cliff/water_depth_300/glacier_min_yield_20_surf_slope_0.02_med_visc_bed_slope_0.0_bump8/glacier_cliff_'

fname_out = fname_base#+'junk'
#fname_out = 'data/cliff/flat/water_depth_306_high_res_flux_0.0_slope-0.0_colder/glacier_cliff_'

#fname_out = 'data/tmp14/glacier_cliff_'

#12: 791
# Specify file to load
step = 1
nmax = 4650

# Geometric variables
if fname_base[23]=='0':
    ice_thick = 135.0
    Hab = ice_thick
elif fname_base[23]=='3':
    ice_thick = 400.0
    Hab = 65.0
    Hab = 60.0
else:
    ice_thick = 800.0
    Hab = 25.0


length= ice_thick*12
water_depth = ice_thick*910.0/1020 - Hab

# Define surface, bottom, etc
if fname_base[5:10] == 'retro':
    surf_slope = 0.02
    bed_slope = -0.01
else:
    surf_slope = 0.02
    bed_slope =  0.01

#surf_slope= 0.02



ice_thick = 400.0
Hab = 60.0
Hab = 45.0
#ice_thick = 135.0
#Hab = ice_thick
#ice_thick = 800.0
#Hab = 25.0
length= ice_thick*12
#length= ice_thick*5
water_depth = ice_thick*910.0/1020 - Hab

#dz = 60.0
dz = round(ice_thick/13.333333333/2)
Nx = int(length/dz)
Nz = int(ice_thick/dz)

# Define surface, bottom, etc
surf_slope =  0.02
bed_slope =   0.01*0

notch_height = 200.0*0
notch_length = 40.0
notch_slope = notch_height/notch_length/2
L = length+ice_thick*0
bump_width = ice_thick
bump_height =  0.25*ice_thick

def bot_fun(x):
   b = -water_depth + bed_slope*(length-x)
   return b

#def bed_fun(x):
 #  b = -water_depth + bed_slope*(length-x) + bump_height*exp(-((L-x)/bump_width)**2)
  # return b

def bed_fun(x):
   b = -water_depth + bed_slope*(length-x) + bump_height*np.exp(-((L-x)/bump_width)**2)
   return b

def ice_thick_fun(x):
    return ice_thick + (surf_slope-bed_slope)*(L-x)

def surf_fun(x):
   s = -water_depth + ice_thick + (surf_slope)*(L-x)#-  notch_slope*(L+notch_length-x)*(x<(L+notch_length))*(x>(L)) + notch_slope*(L-x-notch_length)*(x>(L-notch_length))*(x<=(L))
   return s

xs=linspace(0,length,11)
hs=surf_fun(xs)
bed = -793.7594315851419
max_length = 2.1e3/135.0*ice_thick
xx=linspace(0.0,max_length,1001)
plt.close('all')
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.subplots_adjust(left=-0.0, bottom=.07, right=0.92, top=.95)
#plt.tick_params(axis='both', which='major', labelsize=20)

#for step in xrange(8350,9040,10):
for step in xrange(1640,1890,10):
#for step in xrange(3070,3120,10):


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
   filter =(x< max_length) & (eta>0.0) & (eta_visc>0.0)
   x = x[filter];z=z[filter]
   strain = maximum(data['strain'],1e-16)[filter]
   epsII = data['epsII'][filter]
   t=data['t']
   eta = data['eta'][filter]
   eta_visc = data['eta_visc_m'][filter]
   #xx=data['xx']
   #hsurf=data['hsurf']
   #hbot=data['hbot']
   #width=data['width']


   # Make plot
   pts_size = 0.1
   plt.subplot(2,1,1)
   #c=plt.scatter(x,z,s=pts_size,c=np.log10(epsII),vmin=-8,vmax=-5)
   c=plt.scatter(x,z,s=pts_size,c=np.log10(strain),vmin=-4,vmax=1)
   plt.fill_between([1e3,max_length],0.0,bed_fun(0)-ice_thick/5,color='dodgerblue',alpha=0.33,zorder=0)
   plt.fill_between(xx,bed_fun(xx),bed_fun(0)-ice_thick/5,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   plt.plot(xs,hs,'--',color='Gray',linewidth=2)
   plt.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   ax=plt.gca()
   ax.set_xticks([])
   ax.set_yticks([])
   ax.spines['right'].set_visible(False)
   ax.spines['top'].set_visible(False)
   #ax.spines['left'].set_visible(False)
   #ax.spines['left'].set_linewidth(2)
   ax.spines['bottom'].set_visible(False)
   ax.spines['left'].set_visible(False)
   plt.axis('equal')

   cbaxes = fig.add_axes([0.9, 0.63, 0.015, 0.3])
   cbar2 = plt.colorbar(c, cax = cbaxes,extend='both')
   cbar2.set_ticks([-4,1])
   cbar2.set_ticklabels([r'$10^{-4}$',r'$10^{1}$'])
   cbar2.set_label('$\epsilon_{p}$',fontsize=12)



   text_str = str(int(ice_thick+surf_slope*length+4))+' m'
   text_str2 = str(int(max_length))+' m'

   ax.annotate ('', (-ice_thick/5, bed_fun(0.0)), (-ice_thick/5, bed_fun(0.0)+ice_thick+surf_slope*length+4), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
       text_str, xy=(-ice_thick/5, 120), xycoords='data',
       xytext=(-15, 0), textcoords='offset points',rotation=90)


   yr = int(mod(t,365))
   #day = round((t - yr)*365,1)
   day = (t - yr)*365
   hr,day = modf(day)
   hr = round(hr*24,0)
   title_str = 'Time: '+str(yr).zfill(1)+ 'a '+str(int(day)).zfill(1)+'d '+str(int(hr)).zfill(1)+'hr'
   ax.text(ice_thick/2*10,bed_fun(length/2)+ice_thick+surf_slope*length+20,title_str,fontsize=14)

   plt.subplot(2,1,2)

   c=plt.scatter(x,z,s=pts_size,c=np.log10(epsII),vmin=-8,vmax=-5)
   plt.fill_between([1e3,max_length],0.0,bed_fun(0)-ice_thick/5,color='dodgerblue',alpha=0.33,zorder=0)
   plt.fill_between(xx,bed_fun(xx),bed_fun(0)-ice_thick/5,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   plt.plot(xs,hs,'--',color='Gray',linewidth=2)
   plt.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   ax=plt.gca()
   ax.set_xticks([])
   ax.set_yticks([])
   ax.spines['right'].set_visible(False)
   ax.spines['top'].set_visible(False)
   ax.spines['bottom'].set_visible(False)
   ax.spines['left'].set_visible(False)
   plt.axis('equal')

   cbaxes = fig.add_axes([0.9, 0.18, 0.015, 0.3])
   cbar3 = plt.colorbar(c, cax = cbaxes,extend='both')
   cbar3.set_ticks([-8,-5])
   cbar3.set_ticklabels([r'$10^{-8}$',r'$10^{-5}$'])
   cbar3.set_label('$\dot \epsilon_{II}$ (s$^{-1}$)',fontsize=12)

   ax.annotate ('', (-ice_thick/5, bed_fun(0.0)), (-ice_thick/5, bed_fun(0.0)+ice_thick+surf_slope*length+4), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
      text_str, xy=(-ice_thick/5, 120), xycoords='data',
      xytext=(-15, 0), textcoords='offset points',rotation=90)

   ax.annotate ('', (0.0, bed_fun(0)-ice_thick/3), (int(max_length), bed_fun(0)-ice_thick/3), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
    text_str2, xy=(max_length/2, bed_fun(0)-ice_thick/3), xycoords='data',
    xytext=(0, -10), textcoords='offset points')

   plt.subplots_adjust(wspace=0, hspace=0)
   plt.pause(1e-16)
   plt.show()
   fname_ext = str(step).zfill(3)+'.png'
   #fname_base= 'data/tmp8/glacier_cliff_'
   fname = fname_out+fname_ext
   plt.savefig(fname)
