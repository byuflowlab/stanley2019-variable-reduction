import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from windRoses import *

"""colors"""
robin = '#55dfe6'
nectar = '#e3bf54'
tuscan = '#e37854'
royal = '#a054e3'

robin = 'C0'
nectar = 'C1'
tuscan = 'C1'
royal = 'C3'

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10')

fig = plt.figure(figsize=[6.5,4.5])
ax1 = plt.subplot2grid((2, 3), (0, 0),projection='polar')
ax2 = plt.subplot2grid((2, 3), (0, 1),projection='polar')
ax3 = plt.subplot2grid((2, 3), (0, 2),projection='polar')
ax4 = plt.subplot2grid((2, 3), (1, 0),projection='polar')
ax5 = plt.subplot2grid((2, 3), (1, 1),projection='polar')
ax6 = plt.subplot2grid((2, 3), (1, 2),projection='polar')




nDirections = 24
bottom = 0
width = (2.*np.pi) / float(nDirections)

wd, wf, ws = northIslandRose(nDirections)
wd -= wd[np.argmax(wf*ws**3)]
wd += 270.
# wd +=180./float(nDirections)
for i in range(nDirections):
   wd[i] = radians(wd[i])*-1.
"""wind rose"""
max_height = max(wf)
bars = ax1.bar(wd, wf, width=width, bottom=bottom, color='C2',alpha=0.5, edgecolor='white', linewidth=0.5)
thetaticks = np.arange(0,360,45)
ax1.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax1.set_thetagrids(thetaticks, frac=1.2)
ax1.set_rgrids([0.05,0.1,0.13], angle=-35.)
ax1.set_yticklabels(['0.05','0.1',''],horizontalalignment='center',fontsize=8)

max_height = max(ws)
bars = ax4.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5, edgecolor='white', linewidth=0.5)
thetaticks = np.arange(0,360,45)
ax4.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax4.set_thetagrids(thetaticks, frac=1.2)
ax4.set_rgrids([5,10,13], angle=-35.)
ax4.set_yticklabels(['5 m/s','10 m/s',''],horizontalalignment='center',fontsize=8)





wd, wf, ws = ukiahRose(nDirections)
wd -= wd[np.argmax(wf*ws**3)]
wd += 270.
# wd +=180./float(nDirections)
for i in range(nDirections):
   wd[i] = radians(wd[i])*-1.
"""wind rose"""
max_height = max(wf)
bars = ax2.bar(wd, wf, width=width, bottom=bottom, color='C2',alpha=0.5, edgecolor='white', linewidth=0.5)
thetaticks = np.arange(0,360,45)
ax2.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax2.set_thetagrids(thetaticks, frac=1.2)
ax2.set_rgrids([0.05,0.1,0.13], angle=-35.)
ax2.set_yticklabels(['0.05','0.1',''],horizontalalignment='center',fontsize=8)

max_height = max(ws)
bars = ax5.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5, edgecolor='white', linewidth=0.5)
thetaticks = np.arange(0,360,45)
ax5.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax5.set_thetagrids(thetaticks, frac=1.2)
ax5.set_rgrids([5,10,13], angle=-35.)
ax5.set_yticklabels(['5 m/s','10 m/s',''],horizontalalignment='center',fontsize=8)





wd, wf, ws = victorvilleRose(nDirections)
wd -= wd[np.argmax(wf*ws**3)]
wd += 270.
# wd +=180./float(nDirections)
for i in range(nDirections):
   wd[i] = radians(wd[i])*-1.
"""wind rose"""
max_height = max(wf)
bars = ax3.bar(wd, wf, width=width, bottom=bottom, color='C2',alpha=0.5, edgecolor='white', linewidth=0.5)
# ax9.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=10)

# tick locations
thetaticks = np.arange(0,360,45)
ax3.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax3.set_thetagrids(thetaticks, frac=1.2)
ax3.set_rgrids([0.05,0.1,0.13], angle=-35.)
ax3.set_yticklabels(['0.05','0.1',''],horizontalalignment='center',fontsize=8)

max_height = max(ws)
bars = ax6.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5, edgecolor='white', linewidth=0.5)
thetaticks = np.arange(0,360,45)
ax6.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax6.set_thetagrids(thetaticks, frac=1.2)
ax6.set_rgrids([5,10,13], angle=-35.)
ax6.set_yticklabels(['5 m/s','10 m/s',''],horizontalalignment='center',fontsize=8)




ax4.set_xlabel('North Island',fontsize=8)
ax5.set_xlabel('Ukiah',fontsize=8)
ax6.set_xlabel('Victorville',fontsize=8)

ax1.set_ylabel('wind roses',fontsize=8,labelpad=20)
ax4.set_ylabel('wind speeds',fontsize=8,labelpad=20)

plt.subplots_adjust(top = 0.96, bottom = 0.1, right = 0.99, left = 0.1,
            hspace = 0.2, wspace = 0.35)

ax1.text(0.76*np.pi,0.2,'a',fontsize=11,color='black')
ax2.text(0.76*np.pi,0.2,'b',fontsize=11,color='black')
ax3.text(0.76*np.pi,0.2,'c',fontsize=11,color='black')

# plt.savefig('roses.pdf',transparent=True)
plt.show()
