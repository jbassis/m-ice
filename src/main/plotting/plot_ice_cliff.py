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


water_depth = 700
surf_slope =  0.02
bed_slope =  -0.02
flux = 3.0

# Base directory for file
fname_base = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.01_flux_2.0_high_res/glacier_cliff_'
fname_base = '../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.01_flux_3.5_high_res_CFL/glacier_cliff_'
#fname_base ='../data/cliff/water_depth_'+str(water_depth)+'/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_T_-15.0melange_CFL/glacier_cliff_'
fname_base ='../data/cliff/water_depth_'+str(water_depth)+'/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
#fname_base ='../data/cliff/water_depth_700/glacier_surf_slope_0.02_bed_slope_-0.02_flux_4.0_high_res_T_-15.0melange_CFL/glacier_cliff_'
fname_base = '../data/cliff/water_depth_'+str(water_depth)+'/glacier_surf_slope_0.02_bed_slope_0.0_flux_0.0_high_res_T_-20.0_buttressing25.0kPa_CFL/glacier_cliff_'
fname_base = '../data/buttressing/water_depth_700_buttressing_25.0kPa_removed_50.0day/glacier_surf_slope_0.02_bed_slope_-0.02_flux_2.0_high_res_T_-20.0_CFL/glacier_cliff_'

fname_out = fname_base#+'junk'

#fname_out = '../data/cliff/water_depth_700_super_buoyancy_crevasse/glacier_surf_slope_0.02_bed_slope_0.02_flux_3.0_high_res_T_-20.0_CFL/'
# Specify file to load


#surf_slope= 0.02



ice_thick = 800.0
Hab = 25.0
#Hab = 45.0
#ice_thick = 135.0
#Hab = ice_thick
#ice_thick = 800.0
#ice_thick = 135
#Hab = 25.0
length= ice_thick*12
water_depth = ice_thick*910.0/1020 - Hab
#water_depth = 20
#dz = 60.0
dz = round(ice_thick/13.333333333/2)
Nx = int(length/dz)
Nz = int(ice_thick/dz)

# Define surface, bottom, etc


notch_height = 200.0*0
notch_length = 40.0
notch_slope = notch_height/notch_length/2
L = length+ice_thick*0
bump_width = ice_thick
bump_height =  0.25*ice_thick*0

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
#for step in range(7220,8230,10):
import glob
import os
list_of_files = glob.glob(fname_out+'*.npz') # * means all if need specific format then *.csv
latest_file = max(list_of_files, key=os.path.getctime)
#for step in range(0,680,10):
for latest_file in list_of_files:



   #fname_base= 'data/tmp8/restart_glacier_cliff_'

   fig.clf()


   # File extension
   #fname_ext = 'no_buttressing'+str(step).zfill(3)+'.npz'
   #fname_ext = str(step).zfill(3)+'.npz'


   # Filename
   #fname = fname_base+fname_ext
   fname = latest_file
   print(fname)


   # Load data
   data = load(fname)




   # Triangulation and point coordinates
   x=data['xm'];z=data['zm'];
   filter =(x< max_length)
   x = x[filter];z=z[filter]
   strain = maximum(data['strain'],1e-16)[filter]
   epsII = data['epsII'][filter]
   t=data['t']



   # Make plot
   pts_size = 0.1
   plt.subplot(2,1,1)
   #c=plt.scatter(x,z,s=pts_size,c=np.log10(epsII),vmin=-8,vmax=-5)
   c=plt.scatter(x,z,s=pts_size,c=np.log10(strain),vmin=-4,vmax=1)
   plt.fill_between([1e3,max_length],0.0,np.minimum(bed_fun(0),bed_fun(length))-ice_thick/5,color='dodgerblue',alpha=0.33,zorder=0)
   plt.fill_between(xx,bed_fun(xx),np.minimum(bed_fun(0),bed_fun(length))-ice_thick/5,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
   plt.plot(xs,hs,'--',color='Gray',linewidth=2)
   plt.plot([xs[-1],xs[-1]],[hs[-1],bed_fun(xs[-1])],'--',color='Gray',linewidth=2)
   ax=plt.gca()
   #ax.text(4500,bed_fun(length/2)+ice_thick+surf_slope*length+20,'buttressing',color='red',fontweight='bold',fontsize=14)
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



   text_str = str(int(ice_thick+(surf_slope-bed_slope)*length))+' m'
   text_str2 = str(int(max_length))+' m'

   ax.annotate ('', (-ice_thick/5, bed_fun(0.0)), (-ice_thick/5, bed_fun(0.0)+ice_thick+(surf_slope-bed_slope)*length), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
      text_str, xy=(-ice_thick/5,  0.5*(2*bed_fun(0)+ice_thick+(surf_slope-bed_slope)*length)), xycoords='data',
      xytext=(-15, 0), textcoords='offset points',rotation=90,va='center')


   yr = int(mod(t,365))
   #day = round((t - yr)*365,1)
   day = (t - yr)*365
   hr,day = modf(day)
   hr = round(hr*24,0)
   title_str = 'Time: '+str(yr).zfill(1)+ 'a '+str(int(day)).zfill(1)+'d '+str(int(hr)).zfill(1)+'hr'
   ax.text(ice_thick/2*10,bed_fun(length/2)+ice_thick+surf_slope*length+20,title_str,fontsize=14)

   plt.subplot(2,1,2)

   c=plt.scatter(x,z,s=pts_size,c=np.log10(np.maximum(epsII,1e-16)),vmin=-8,vmax=-5)
   plt.fill_between([1e3,max_length],0.0,np.minimum(bed_fun(0),bed_fun(length))-ice_thick/5,color='dodgerblue',alpha=0.33,zorder=0)
   plt.fill_between(xx,bed_fun(xx),np.minimum(bed_fun(0),bed_fun(length))-ice_thick/5,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
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

   ax.annotate ('', (-ice_thick/5, bed_fun(0.0)), (-ice_thick/5, bed_fun(0.0)+ice_thick+(surf_slope-bed_slope)*length), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
      text_str, xy=(-ice_thick/5,  0.5*(2*bed_fun(0)+ice_thick+(surf_slope-bed_slope)*length)), xycoords='data',
      xytext=(-15, 0), textcoords='offset points',rotation=90,va='center')

   ax.annotate ('', (0.0, np.minimum(bed_fun(0),bed_fun(length))-ice_thick/3), (int(max_length), np.minimum(bed_fun(0),bed_fun(length))-ice_thick/3), arrowprops={'arrowstyle':'<->','linewidth':2})
   ax.annotate(
    text_str2, xy=(max_length/2, np.minimum(bed_fun(0),bed_fun(length))-ice_thick/3), xycoords='data',
    xytext=(0, -10), textcoords='offset points')



   plt.subplots_adjust(wspace=0, hspace=0)
   plt.pause(1e-16)
   plt.show()

   #fname_ext = str(step).zfill(3)+'.png'
   #fname_base= 'data/tmp8/glacier_cliff_'
   #fname = fname_out+fname_ext

   fname = fname[0:-3]+'.png'
   plt.savefig(fname)
