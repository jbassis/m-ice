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

gradH = [-0.04,-0.06,-0.18]
T =[-20,-10,-5]
plt.errorbar(T,gradH,xerr=2.5,yerr=0.01,color='k',fmt='s',ecolor='gray',capsize=5)
plt.xlabel(u'Ice temperature (\u00B0C)')
plt.ylabel(r'Critical Gradient')
p=np.polyfit(T,gradH,3)
Tp=np.linspace(-25,0,11)
f=np.polyval(p,Tp)
#plt.plot(Tp,f,'--k')
#plt.fill_between(Tp,f,-0.2,color='lightgray')
#plt.text(-22,-0.15,'Collapse')
#plt.text(-10,-0.05,'Collapse')
plt.ylim([-0.2,0.0])
plt.xlim([-25,-0.0])

plt.show()
#ax.legend()
