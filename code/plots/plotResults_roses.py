import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
from windRoses import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes


"""colors"""
robin = '#55dfe6'
nectar = '#e3bf54'
tuscan = '#e37854'
royal = '#a054e3'

robin = 'C0'
nectar = 'C1'
tuscan = 'C1'
royal = 'C3'

# robin = '#ACBD78'

nTurbs = 100
rotorDiameter = np.ones(nTurbs)*130.

spacing = 4.
side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
a = side_length**2
circle_radius = np.sqrt(a/np.pi)

folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_4.0_northIslandRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_4.0_northIslandRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_4.0_northIslandRose'

with open('%s/AEP.txt'%folder_grid) as my_file:
    AEP_grid4 = my_file.readlines()
for p in range(len(AEP_grid4)):
    AEP_grid4[p] = float(AEP_grid4[p])
with open('%s/funcCalls.txt'%folder_grid) as my_file:
    f_grid4 = my_file.readlines()
for p in range(len(f_grid4)):
    f_grid4[p] = float(f_grid4[p])

with open('%s/AEP.txt'%folder_direct) as my_file:
    AEP_direct4 = my_file.readlines()
for p in range(len(AEP_direct4)):
    AEP_direct4[p] = float(AEP_direct4[p])
with open('%s/funcCalls.txt'%folder_direct) as my_file:
    f_direct4 = my_file.readlines()
for p in range(len(f_direct4)):
    f_direct4[p] = float(f_direct4[p])

with open('%s/AEP.txt'%folder_param) as my_file:
    AEP_param4 = my_file.readlines()
for p in range(len(AEP_param4)):
    AEP_param4[p] = float(AEP_param4[p])
with open('%s/funcCalls.txt'%folder_param) as my_file:
    f_param4 = my_file.readlines()
for p in range(len(f_param4)):
    f_param4[p] = float(f_param4[p])


folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_4.0_ukiahRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_4.0_ukiahRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_4.0_ukiahRose'

with open('%s/AEP.txt'%folder_grid) as my_file:
    AEP_grid6 = my_file.readlines()
for p in range(len(AEP_grid6)):
    AEP_grid6[p] = float(AEP_grid6[p])
with open('%s/funcCalls.txt'%folder_grid) as my_file:
    f_grid6 = my_file.readlines()
for p in range(len(f_grid6)):
    f_grid6[p] = float(f_grid6[p])

with open('%s/AEP.txt'%folder_direct) as my_file:
    AEP_direct6 = my_file.readlines()
for p in range(len(AEP_direct6)):
    AEP_direct6[p] = float(AEP_direct6[p])
with open('%s/funcCalls.txt'%folder_direct) as my_file:
    f_direct6 = my_file.readlines()
for p in range(len(f_direct6)):
    f_direct6[p] = float(f_direct6[p])

with open('%s/AEP.txt'%folder_param) as my_file:
    AEP_param6 = my_file.readlines()
for p in range(len(AEP_param6)):
    AEP_param6[p] = float(AEP_param6[p])
with open('%s/funcCalls.txt'%folder_param) as my_file:
    f_param6 = my_file.readlines()
for p in range(len(f_param6)):
    f_param6[p] = float(f_param6[p])


folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_4.0_victorvilleRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_4.0_victorvilleRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_4.0_victorvilleRose'

with open('%s/AEP.txt'%folder_grid) as my_file:
    AEP_grid8 = my_file.readlines()
for p in range(len(AEP_grid8)):
    AEP_grid8[p] = float(AEP_grid8[p])
with open('%s/funcCalls.txt'%folder_grid) as my_file:
    f_grid8 = my_file.readlines()
for p in range(len(f_grid8)):
    f_grid8[p] = float(f_grid8[p])

with open('%s/AEP.txt'%folder_direct) as my_file:
    AEP_direct8 = my_file.readlines()
for p in range(len(AEP_direct8)):
    AEP_direct8[p] = float(AEP_direct8[p])
with open('%s/funcCalls.txt'%folder_direct) as my_file:
    f_direct8 = my_file.readlines()
for p in range(len(f_direct8)):
    f_direct8[p] = float(f_direct8[p])

with open('%s/AEP.txt'%folder_param) as my_file:
    AEP_param8 = my_file.readlines()
for p in range(len(AEP_param8)):
    AEP_param8[p] = float(AEP_param8[p])
with open('%s/funcCalls.txt'%folder_param) as my_file:
    f_param8 = my_file.readlines()
for p in range(len(f_param8)):
    f_param8[p] = float(f_param8[p])

AEP_grid4 = AEP_grid4/np.max(AEP_direct4)
AEP_param4 = AEP_param4/np.max(AEP_direct4)
AEP_direct4 = AEP_direct4/np.max(AEP_direct4)

AEP_grid6 = AEP_grid6/np.max(AEP_direct6)
AEP_param6 = AEP_param6/np.max(AEP_direct6)
AEP_direct6 = AEP_direct6/np.max(AEP_direct6)

AEP_grid8 = AEP_grid8/np.max(AEP_direct8)
AEP_param8 = AEP_param8/np.max(AEP_direct8)
AEP_direct8 = AEP_direct8/np.max(AEP_direct8)

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10')
fig = plt.figure(figsize=[6.5,6.])
ax4 = plt.subplot2grid((3, 3), (1, 0))
ax5 = plt.subplot2grid((3, 3), (1, 1),sharey=ax4,sharex=ax4)
ax6 = plt.subplot2grid((3, 3), (1, 2),sharey=ax4,sharex=ax4)
ax1 = plt.subplot2grid((3, 3), (2, 0))
ax2 = plt.subplot2grid((3, 3), (2, 1),sharey=ax1,sharex=ax1)
ax3 = plt.subplot2grid((3, 3), (2, 2),sharey=ax1,sharex=ax1)


# ax10 = plt.subplot2grid((4, 3), (1, 0), polar=True)
# ax11 = plt.subplot2grid((4, 3), (1, 1), polar=True)
# ax12 = plt.subplot2grid((4, 3), (1, 2), polar=True)

ax7 = plt.subplot2grid((3, 3), (0, 0), polar=True)
ax8 = plt.subplot2grid((3, 3), (0, 1), polar=True)
ax9 = plt.subplot2grid((3, 3), (0, 2), polar=True)


nDirections = 30
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
bars = ax7.bar(wd, wf, width=width, bottom=bottom, color=royal,alpha=0.5)
thetaticks = np.arange(0,360,45)
ax7.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax7.set_thetagrids(thetaticks, frac=1.2)
ax7.set_rgrids([0.04,0.08], angle=-35.)
ax7.set_yticklabels(['0.04','0.08'],horizontalalignment='center',fontsize=8)
ax7.set_ylim(0.,0.1)

# max_height = max(ws)
# bars = ax10.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5)
# ax10.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=9)
# ax10.set_thetagrids(thetaticks, frac=1.2)
# ax10.set_rgrids([5,10], angle=-35.)
# ax10.set_yticklabels(['5','10 m/s'],horizontalalignment='left',fontsize=9)
# ax10.set_ylim(0.,13.)

wd, wf, ws = ukiahRose(nDirections)
wd -= wd[np.argmax(wf*ws**3)]
wd += 270.
# wd +=180./float(nDirections)
for i in range(nDirections):
   wd[i] = radians(wd[i])*-1.
"""wind rose"""
max_height = max(wf)
bars = ax8.bar(wd, wf, width=width, bottom=bottom, color=royal,alpha=0.5)
thetaticks = np.arange(0,360,45)
ax8.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax8.set_thetagrids(thetaticks, frac=1.2)
ax8.set_rgrids([0.04,0.08], angle=-35.)
ax8.set_yticklabels(['0.04','0.08'],horizontalalignment='center',fontsize=8)
ax8.set_ylim(0.,0.1)

# max_height = max(ws)
# bars = ax11.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5)
# ax11.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=9)
# ax11.set_thetagrids(thetaticks, frac=1.2)
# ax11.set_rgrids([5,10], angle=-35.)
# ax11.set_yticklabels(['5','10 m/s'],horizontalalignment='left',fontsize=9)
# ax11.set_ylim(0.,13.)


wd, wf, ws = victorvilleRose(nDirections)
wd -= wd[np.argmax(wf*ws**3)]
wd += 270.
# wd +=360./float(nDirections)
for i in range(nDirections):
   wd[i] = radians(wd[i])*-1.
"""wind rose"""
max_height = max(wf)
bars = ax9.bar(wd, wf, width=width, bottom=bottom, color=royal,alpha=0.5)
# ax9.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=9)

# tick locations
thetaticks = np.arange(0,360,45)
ax9.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
ax9.set_thetagrids(thetaticks, frac=1.2)
ax9.set_rgrids([0.04,0.08], angle=-35.)
ax9.set_yticklabels(['0.04','0.08'],horizontalalignment='center',fontsize=8)
ax9.set_ylim(0.,0.1)

# max_height = max(ws)
# bars = ax12.bar(wd, ws, width=width, bottom=bottom, color=tuscan,alpha=0.5)
# ax12.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=9)
# ax12.set_thetagrids(thetaticks, frac=1.2)
# ax12.set_rgrids([5,10], angle=-35.)
# ax12.set_yticklabels(['5','10 m/s'],horizontalalignment='left',fontsize=9)
# ax12.set_ylim(0.,13.)

ax4.set_title('North Island',fontsize=8,rotation=0.)
ax5.set_title('Ukiah',fontsize=8,rotation=0.)
ax6.set_title('Victorville',fontsize=8,rotation=0.)

ax7.axhline(linewidth=0.5)
ax8.axhline(linewidth=0.5)
ax9.axhline(linewidth=0.5)
# ax10.axhline(linewidth=0.5)
# ax11.axhline(linewidth=0.5)
# ax12.axhline(linewidth=0.5)
ax7.grid(linewidth=0.5)
ax8.grid(linewidth=0.5)
ax9.grid(linewidth=0.5)
# ax10.grid(linewidth=0.5)
# ax11.grid(linewidth=0.5)
# ax12.grid(linewidth=0.5)

ax8.set_xlabel('wind roses',fontsize=8,labelpad=2)
# ax10.set_ylabel('wind speeds',fontsize=9,labelpad=15)




# min4 = np.min(np.array([np.min(f_param4),np.min(f_direct4),np.min(f_grid4),np.min(f_param6),np.min(f_direct6),np.min(f_grid6),np.min(f_param8),np.min(f_direct8),np.min(f_grid8)]))
# min6 = np.min(np.array([np.min(f_param6),np.min(f_direct6),np.min(f_grid6)]))
# min8 = np.min(np.array([np.min(f_param8),np.min(f_direct8),np.min(f_grid8)]))
min4 = 10.
min6 = min4
min8 = min4

# max4 = np.max(np.array([np.max(f_param4),np.max(f_direct4),np.max(f_grid4),np.max(f_param6),np.max(f_direct6),np.max(f_grid6),np.max(f_param8),np.max(f_direct8),np.max(f_grid8)]))
# max6 = np.max(np.array([np.max(f_param6),np.max(f_direct6),np.max(f_grid6)]))
# max8 = np.max(np.array([np.max(f_param8),np.max(f_direct8),np.max(f_grid8)]))
max4 = 2000
max6 = max4
max8 = max4

b4 = np.logspace(np.log10(min4),np.log10(max4),25)
b6 = np.logspace(np.log10(min6),np.log10(max6),25)
b8 = np.logspace(np.log10(min8),np.log10(max8),25)

num = 100

ax1.hist(f_grid4[0:num],bins=b4,color=tuscan,alpha=0.5)
ax2.hist(f_grid6[0:num],bins=b6,color=tuscan,alpha=0.5)
ax3.hist(f_grid8[0:num],bins=b8,color=tuscan,alpha=0.5)

ax1.hist(f_param4[0:num],bins=b4,color=robin,alpha=0.5)
ax2.hist(f_param6[0:num],bins=b6,color=robin,alpha=0.5)
ax3.hist(f_param8[0:num],bins=b8,color=robin,alpha=0.5)

ax1.hist(f_direct4[0:num],bins=b4,color=royal,alpha=0.5)
ax2.hist(f_direct6[0:num],bins=b6,color=royal,alpha=0.5)
ax3.hist(f_direct8[0:num],bins=b8,color=royal,alpha=0.5)

ax1.hist(f_param4[0:num],bins=b4,color=robin,alpha=0.5)
ax2.hist(f_param6[0:num],bins=b6,color=robin,alpha=0.5)
ax3.hist(f_param8[0:num],bins=b8,color=robin,alpha=0.5)

min4 = np.min(np.array([np.min(AEP_param4),np.min(AEP_direct4),np.min(AEP_grid4)]))
min6 = np.min(np.array([np.min(AEP_param6),np.min(AEP_direct6),np.min(AEP_grid6)]))
min8 = np.min(np.array([np.min(AEP_param8),np.min(AEP_direct8),np.min(AEP_grid8)]))

max4 = np.max(np.array([np.max(AEP_param4),np.max(AEP_direct4),np.max(AEP_grid4)]))
max6 = np.max(np.array([np.max(AEP_param6),np.max(AEP_direct6),np.max(AEP_grid6)]))
max8 = np.max(np.array([np.max(AEP_param8),np.max(AEP_direct8),np.max(AEP_grid8)]))

min4 = 0.9
min6 = min4
min8 = min4

max4 = 1.025
max6 = max4
max8 = max4

b4 = np.linspace(min4,max4,25)
b6 = np.linspace(min6,max6,25)
b8 = np.linspace(min8,max8,25)

# b6 = np.linspace(min4,max4,25)
# b8 = np.linspace(min4,max4,25)

# b4 = np.linspace(0.9,1.05,25)
# b6 = np.linspace(0.9,1.05,25)
# b8 = np.linspace(0.9,1.05,25)


ax4.hist(AEP_grid4[0:num],bins=b4,color=tuscan,alpha=0.5,label='grid')
ax5.hist(AEP_grid6[0:num],bins=b6,color=tuscan,alpha=0.5)
ax6.hist(AEP_grid8[0:num],bins=b8,color=tuscan,alpha=0.5)
ax4.hist(AEP_param4[0:num],bins=b4,color=robin,alpha=0.5,label='parameterized')
ax5.hist(AEP_param6[0:num],bins=b6,color=robin,alpha=0.5)
ax6.hist(AEP_param8[0:num],bins=b8,color=robin,alpha=0.5)
ax4.hist(AEP_direct4[0:num],bins=b4,color=royal,alpha=0.5,label='direct')
ax5.hist(AEP_direct6[0:num],bins=b6,color=royal,alpha=0.5)
ax6.hist(AEP_direct8[0:num],bins=b8,color=royal,alpha=0.5)
ax4.hist(AEP_param4[0:num],bins=b4,color=robin,alpha=0.5)
ax5.hist(AEP_param6[0:num],bins=b6,color=robin,alpha=0.5)
ax6.hist(AEP_param8[0:num],bins=b8,color=robin,alpha=0.5)

ax4.legend(loc=2,prop={'size': 8})

ax5.set_xlabel('normalized optimal AEP',fontsize=8)
ax2.set_xlabel('function calls', fontsize=8)

ax4.set_ylabel('# optimization runs',fontsize=8)
ax1.set_ylabel('# optimization runs', fontsize=8)


windDirections, windFrequencies, windSpeeds = ukiahRose(30)
print windDirections
print windFrequencies

windDirections = 90.-windDirections
for i in range(len(windDirections)):
    windDirections[i] = np.deg2rad(windDirections[i])

ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')

ticks = np.array([10.,100.,1000.,2000.])
ax1.set_xticks(ticks)
ax1.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))
ax2.set_xticks(ticks)
ax2.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))
ax3.set_xticks(ticks)
ax3.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))


ticks = np.array([0.9,0.95,1.0])
ax4.set_xticks(ticks)
ax4.set_xticklabels(('0.9','0.95','1.0'))
ax5.set_xticks(ticks)
ax5.set_xticklabels(('0.9','0.95','1.0'))
ax6.set_xticks(ticks)
ax6.set_xticklabels(('0.9','0.95','1.0'))
ax1.set_xlim(10.,5000.)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('bottom')
ax5.spines['right'].set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.yaxis.set_ticks_position('left')
ax5.xaxis.set_ticks_position('bottom')
ax6.spines['right'].set_visible(False)
ax6.spines['top'].set_visible(False)
ax6.yaxis.set_ticks_position('left')
ax6.xaxis.set_ticks_position('bottom')

ax4.set_xlim(0.9,1.01)
ax5.set_xlim(0.9,1.01)
ax6.set_xlim(0.9,1.01)


plt.subplots_adjust(top = 0.94, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.6, wspace = 0.2)


ax4.text(0.907,55.,'a',fontsize=11,color='black')
ax5.text(0.907,55.,'b',fontsize=11,color='black')
ax6.text(0.907,55.,'c',fontsize=11,color='black')
ax1.text(12.5,27.,'d',fontsize=11,color='black')
ax2.text(12.5,27.,'e',fontsize=11,color='black')
ax3.text(12.5,27.,'f',fontsize=11,color='black')

ax1.tick_params(axis='both', which='major', labelsize=8)
ax1.tick_params(axis='both', which='minor', labelsize=8)
ax2.tick_params(axis='both', which='major', labelsize=8)
ax2.tick_params(axis='both', which='minor', labelsize=8)
ax3.tick_params(axis='both', which='major', labelsize=8)
ax3.tick_params(axis='both', which='minor', labelsize=8)
ax4.tick_params(axis='both', which='major', labelsize=8)
ax4.tick_params(axis='both', which='minor', labelsize=8)
ax5.tick_params(axis='both', which='major', labelsize=8)
ax5.tick_params(axis='both', which='minor', labelsize=8)
ax6.tick_params(axis='both', which='major', labelsize=8)
ax6.tick_params(axis='both', which='minor', labelsize=8)

plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/results_roses.pdf',transparent=True)
plt.show()
