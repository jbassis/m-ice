"""
Plot critical thickness gradient as a function of ice temperature

"""

import numpy as np
import pylab as plt

import matplotlib
matplotlib.style.use('seaborn-colorblind')


plt.close('all')
width  = 3.5
height = width / 2.5
fig=plt.figure(num=1, figsize=(width, height), facecolor='w', edgecolor='k')
fig.set_size_inches(width,height,forward=True)
ax = fig.add_axes([0.2, 0.3, 0.75, 0.6])

gradH = [-0.04,-0.08,-0.1,-0.19]
T =[-20,-15,-10,-5]
plt.errorbar(T,gradH,xerr=2.5,yerr=0.01,color='k',fmt='s',ecolor='lightgray',capsize=5)
plt.xlabel(u'Ice temperature (\u00B0C)')
plt.ylabel('Critical thickness gradient')
plt.show()
#ax.legend()
