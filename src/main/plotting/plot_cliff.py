"""
Plot 2D rift profiles
"""
#%matplotlib inline
import pylab as plt
from numpy import *
import numpy as np

def plot_snap(fig,water_depth,surf_slope,bed_slope,flux,step,annotate_left=True):
    fname_base ='../data/cliff/water_depth_'+str(water_depth)+'/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_T_-20.0_CFL/glacier_cliff_'


    if water_depth == 0:
        ice_thick = 135
    elif water_depth == 700:
        ice_thick = 800
        fname_base = '../data/cliff/water_depth_700/glacier_surf_slope_'+str(surf_slope)+'_bed_slope_'+str(bed_slope)+'_flux_'+str(flux)+'_high_res_CFL/glacier_cliff_'
    else:
        ice_thick = 400.0
    length= ice_thick*12
    L = length

    def bed_fun(x):
       b = -water_depth + bed_slope*(length-x)
       return b


    def surf_fun(x):
       s = -water_depth + ice_thick + (surf_slope)*(L-x)
       return s


    max_length = 2.1e3/135.0*ice_thick
    max_length = length+ice_thick
    min_length = length-3.5*ice_thick
    xs=linspace(min_length,length,11)
    hs=surf_fun(xs)
    xx=linspace(min_length,max_length,1001)



    # File extension
    fname_ext = str(step).zfill(3)+'.npz'


    # Filename
    fname = fname_base+fname_ext
    print(fname)


    # Load data
    data = load(fname)


    # Triangulation and point coordinates
    x=data['xm'];z=data['zm'];
    filter =(x< max_length) & (x>=min_length)
    x = x[filter];z=z[filter]
    strain = maximum(data['strain'],1e-16)[filter]
    epsII = data['epsII'][filter]
    t=data['t']



    # Make plot
    pts_size = 0.1
    c=plt.scatter(x,z,s=pts_size,c=np.log10(np.maximum(epsII,1e-16)),vmin=-8,vmax=-5,rasterized=True)
    plt.fill_between([min_length,max_length],0.0,np.minimum(bed_fun(min_length),bed_fun(length))-ice_thick/5,color='dodgerblue',alpha=0.33,zorder=0)
    plt.fill_between(xx,bed_fun(xx),np.minimum(bed_fun(min_length),bed_fun(length))-ice_thick/5,color=[139/255.,115/255.,85/255.],alpha=1.0,zorder=1)
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

    """
    cbaxes = fig.add_axes([0.9, 0.18, 0.015, 0.3])
    cbar3 = plt.colorbar(c, cax = cbaxes,extend='both')
    cbar3.set_ticks([-8,-5])
    cbar3.set_ticklabels([r'$10^{-8}$',r'$10^{-5}$'])
    cbar3.set_label('$\dot \epsilon_{II}$ (s$^{-1}$)',fontsize=12)
    """
    yr = int(mod(t,365))
    #day = round((t - yr)*365,1)
    day = (t - yr)*365
    #hr,day = modf(day)
    #hr = round(hr*24,0)
    title_str = 'Day '+str(int(day)).zfill(1)
    ax.text(0.5*(min_length+max_length),bed_fun(length/2)+ice_thick+surf_slope*(length-min_length)+30,title_str,fontsize=10,ha='center')


    if annotate_left == True:
        #height =ice_thick+(surf_slope-bed_slope)*(length-min_length)
        height = ice_thick
        height_str = str(int(height))+ ' m'

        # Annotate ice thickness with arrow
        ax.annotate("", xy=(min_length-ice_thick/10, bed_fun(min_length)), xytext=(min_length-ice_thick/10, bed_fun(min_length)+height),
                     arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",linewidth=2))
        # Label arrow with height
        mid_point = 0.5*(bed_fun(min_length) + bed_fun(min_length)+height)
        ax.text(min_length-ice_thick/2,mid_point,height_str,va='center',rotation=90,clip_on=False)



    # Annotate length with arrow
    length =max_length-min_length
    length_str = str(int(length))+ ' m'
    offset = ice_thick/3
    ax.annotate("", xy=(min_length, bed_fun(min_length)-offset), xytext=(max_length, bed_fun(min_length)-offset),
                 arrowprops=dict(arrowstyle="<->", connectionstyle="arc3",linewidth=2))
    # Label arrow with length
    mid_point = 0.5*(min_length+max_length)
    #ax.text(mid_point,bed_fun(min_length)-2*offset,length_str,ha='center',clip_on=False,zorder=100)
    ax.annotate(length_str,
            xy=(mid_point,bed_fun(min_length)), xytext =(mid_point,bed_fun(min_length)-2*offset),
            xycoords='data',ha='center',annotation_clip=False,clip_on=False)

    return c


    #plt.subplots_adjust(wspace=0, hspace=0)
    #plt.pause(1e-16)
    #plt.show()
    #fname_ext = str(step).zfill(3)+'.png'
    #fname = fname_base+fname_ext
    #print(fname)
color = '#CC6677'
color2 = '#88CCEE'


def plot_term_buttressing(ax,fname1=None,fname2=None,fname3=None,fname4=None,fname5=None,tmax=200):
    if fname1 !=None:
        data3=np.load(fname1)
        t = data3['t']
        yr = np.array((np.mod(t,365)),dtype=int)
        day = (t - yr)*365
        filter = day<=365
        L = data3['L'][0]
        p3=ax.plot(day[filter],(data3['L'][filter]-L)/1e3,'-',color=color,linewidth=2)
        print(p3[0].get_color())
    if fname2 !=None:
        data4=np.load(fname2)
        t = data4['t']
        yr = np.array((np.mod(t,365)),dtype=int)
        day = (t - yr)*365
        filter = day<=365
        L = data4['L'][0]
        p4=ax.plot(day[filter],(data4['L'][filter]-L)/1e3,'-',color='#009E73',linewidth=2)
        print(p4[0].get_color())
    if fname3 !=None:
        data5=np.load(fname3)
        L = data5['L'][0]
        t = data5['t']
        yr = np.array((np.mod(t,365)),dtype=int)
        day = (t - yr)*365
        filter = day<=365
        p5=ax.plot(day[filter],(data5['L'][filter]-L)/1e3,'-',color='#D55E00',linewidth=2)
        print(p5[0].get_color())
    if fname4 !=None:
        data6=np.load(fname4)
        L = data6['L'][0]
        t = data6['t']
        yr = np.array((np.mod(t,365)),dtype=int)
        day = (t - yr)*365
        filter = day<=365
        p6=ax.plot(day[filter],(data6['L'][filter]-L)/1e3,'--',color=color2,linewidth=2)
        print('Max day',day[-1])
    if fname5 !=None:
        data7=np.load(fname5)
        L = data7['L'][0]
        t = data7['t']
        yr = np.array((np.mod(t,365)),dtype=int)
        day = (t - yr)*365
        filter = day<=365
        p7=ax.plot(day[filter],(data7['L'][filter]-L)/1e3,'--',color=color2,linewidth=2)
        #plt.text(65,-0.8,u'800 m',color=p5[0].get_color(),fontsize=8,fontweight='bold')
    #plt.xlabel('Day',fontsize=10)
    plt.ylabel(r'$\Delta L$ (km)',fontsize=10,labelpad=-15)
    plt.ylim([-2,0.25])
    plt.xlim([0.0,tmax])


    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    #ax.xaxis.tick_top()
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    ax.spines['left'].set_position(('outward', 10))
    #ax.spines['left'].set_position('center')
    #ax.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    #ax.spines['right'].set_color('none')
    #ax.spines['top'].set_color('none')

    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


"""
water_depth = 290
surf_slope = 0.02
bed_slope = 0.0
flux = 0.0
step = 50
step = 130
step = 250

plt.clf()
ax = plt.gca()
fig = plt.figure(num=1,facecolor='w', edgecolor='k')

plot_snap(fig,water_depth,surf_slope,bed_slope,flux,step)

ax1.annotate("", xy=(min_length, 0.0), xytext=(min_length, 144),
             arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"))
"""
