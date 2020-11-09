"""
Plot mega-figure showing all simulations
"""

import numpy as np
import pylab as plt

from importlib import reload


import matplotlib
matplotlib.style.use('seaborn-colorblind')

import matplotlib.gridspec as gridspec


import plot_cliff
reload(plot_cliff)
from plot_cliff import *

# Define plot dimensions
fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 7.2#fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height*4.5/5] # 3 rows of figures 1 column width

params = {
          'axes.labelsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': False,
          'figure.figsize': fig_size}


plt.rcParams.update(params)

# Create new figure with correct dimensions
plt.close('all')
fig = plt.figure(1)
fig.patch.set_alpha(0)


snap_height = 0.3
snap_start = 0.025
snap_end = 1
snap_start_vert = 1.0-snap_height
snap_gap_v = 0.0
snap_gap = -0.02
snap_width = (snap_end-snap_start-2*snap_gap)/3

water_depth = 0
surf_slope = 0.02
bed_slope = 0.0
flux = 0.0
ax1 = fig.add_axes([snap_start, snap_start_vert, snap_width, snap_height])
ax1.patch.set_alpha(0)
ax1.text(0.05, 0.825, 'A', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes,fontweight='bold')
plot_snap(ax1,water_depth,surf_slope,bed_slope,flux,130,annotate_left=True)

ax2 = fig.add_axes([snap_start+snap_width+snap_gap, snap_start_vert, snap_width, snap_height])
ax2.patch.set_alpha(0)
ax2.text(0.05, 0.825, 'B', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes,fontweight='bold')
plot_snap(ax2,water_depth,surf_slope,bed_slope,flux,300,annotate_left=False)

ax3 = fig.add_axes([snap_start+2*snap_width+2*snap_gap, snap_start_vert, snap_width, snap_height])
ax3.patch.set_alpha(0)
ax3.text(0.05, 0.825, 'C', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes,fontweight='bold')
plot_snap(ax3,water_depth,surf_slope,bed_slope,flux,680,annotate_left=False)


water_depth = 290
ax4 = fig.add_axes([snap_start, snap_start_vert-snap_height-snap_gap_v, snap_width, snap_height])
ax4.patch.set_alpha(0)
ax4.text(0.05, 0.825, 'D', horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes,fontweight='bold')
plot_snap(ax1,water_depth,surf_slope,bed_slope,flux,50,annotate_left=True)

ax5 = fig.add_axes([snap_start+snap_width+snap_gap, snap_start_vert-snap_height-snap_gap_v, snap_width, snap_height])
ax5.patch.set_alpha(0)
ax5.text(0.05, 0.825, 'E', horizontalalignment='center',verticalalignment='center', transform=ax5.transAxes,fontweight='bold')
plot_snap(ax2,water_depth,surf_slope,bed_slope,flux,110,annotate_left=False)


ax6 = fig.add_axes([snap_start+2*snap_width+2*snap_gap, snap_start_vert-snap_height-snap_gap_v, snap_width, snap_height])
ax6.patch.set_alpha(0)
ax6.text(0.05, 0.825, 'F', horizontalalignment='center',verticalalignment='center', transform=ax6.transAxes,fontweight='bold')
plot_snap(ax3,water_depth,surf_slope,bed_slope,flux,130,annotate_left=False)


water_depth = 700
ax7 = fig.add_axes([snap_start, snap_start_vert-2*snap_height-2*snap_gap_v, snap_width, snap_height])
ax7.patch.set_alpha(0)
ax7.text(0.05, 0.825, 'G', horizontalalignment='center',verticalalignment='center', transform=ax7.transAxes,fontweight='bold')
plot_snap(ax7,water_depth,surf_slope,bed_slope,flux,70,annotate_left=True)

ax8 = fig.add_axes([snap_start+snap_width+snap_gap, snap_start_vert-2*snap_height-2*snap_gap_v, snap_width, snap_height])
ax8.patch.set_alpha(0)
ax8.text(0.05, 0.825, 'H', horizontalalignment='center',verticalalignment='center', transform=ax8.transAxes,fontweight='bold')
plot_snap(ax8,water_depth,surf_slope,bed_slope,flux,120,annotate_left=False)

ax9 = fig.add_axes([snap_start+2*snap_width+2*snap_gap, snap_start_vert-2*snap_height-2*snap_gap_v, snap_width, snap_height])
ax9.patch.set_alpha(0)
ax9.text(0.05, 0.825, 'I', horizontalalignment='center',verticalalignment='center', transform=ax9.transAxes,fontweight='bold')
c = plot_snap(ax9,water_depth,surf_slope,bed_slope,flux,730,annotate_left=False)

cbaxes = fig.add_axes([0.1, 0.07, 0.8, 0.025])
cbar3 = plt.colorbar(c, cax = cbaxes,extend='both',orientation="horizontal")
cbar3.set_ticks([-8,-5])
cbar3.set_ticklabels([r'$10^{-8}$',r'$10^{-5}$'])
cbar3.set_label('$\dot \epsilon_{II}$ (s$^{-1}$)',fontsize=10,labelpad=-17)



plt.show()
plt.savefig('snapshots.pdf')
