import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from windRoses import *

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


folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/circle_4.0_northIslandRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/circle_4.0_northIslandRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/circle_4.0_northIslandRose'

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


folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/square_4.0_northIslandRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/square_4.0_northIslandRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/square_4.0_northIslandRose'

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


ax7 = plt.subplot2grid((3, 3), (0, 0))
ax8 = plt.subplot2grid((3, 3), (0, 1))
ax9 = plt.subplot2grid((3, 3), (0, 2))

ax7.axis('equal')
ax8.axis('equal')
ax9.axis('equal')

windDirections, windFrequencies, windSpeeds = northIslandRose(30)
wind_angle = windDirections[np.argmax(windFrequencies)]




"""amalia boundary"""
locations = np.loadtxt('layout_amalia.txt')
xBounds = locations[:, 0]
yBounds = locations[:, 1]
xBounds = xBounds - min(xBounds) - (max(xBounds)-min(xBounds))/2.
yBounds = yBounds - min(yBounds) - (max(yBounds)-min(yBounds))/2.
locations[:, 0] = xBounds
locations[:, 1] = yBounds
boundaryVertices, boundaryNormals = calculate_boundary(locations)
xBounds = boundaryVertices[:, 0]
yBounds = boundaryVertices[:, 1]
nBounds = len(xBounds)
points = np.zeros((nBounds,2))
points[:, 0] = xBounds
points[:, 1] = yBounds
hull = sp.spatial.ConvexHull(points)
area = hull.volume
area_ratio = area/(np.pi*circle_radius**2)
xBounds = xBounds/np.sqrt(area_ratio)
yBounds = yBounds/np.sqrt(area_ratio)

xBounds = np.append(xBounds,xBounds[0])
yBounds = np.append(yBounds,yBounds[0])

ax7.plot(xBounds,yBounds,'-',color='C0',alpha=0.5)


"""circle boundary"""
nBounds = 200
xBounds = np.zeros(nBounds)
yBounds = np.zeros(nBounds)
theta = np.linspace(0.,2.*np.pi-2.*np.pi/float(nBounds),nBounds)
for i in range(nBounds):
    xBounds[i] = circle_radius*np.cos(theta[i])
    yBounds[i] = circle_radius*np.sin(theta[i])

x = np.zeros_like(xBounds)
x[:] = xBounds[:]
y = np.zeros_like(yBounds)
y[:] = yBounds[:]
xBounds = x*np.cos(np.deg2rad(wind_angle)) - y*np.sin(np.deg2rad(wind_angle))
yBounds = x*np.sin(np.deg2rad(wind_angle)) + y*np.cos(np.deg2rad(wind_angle))

xBounds = np.append(xBounds,xBounds[0])
yBounds = np.append(yBounds,yBounds[0])

ax8.plot(xBounds,yBounds,'-',color='C0',alpha=0.5)


"""square boundary rotated 30 deg from dominant wind direction"""
nBounds = 4
x = np.array([-side_length/2.,side_length/2.,side_length/2.,-side_length/2.])
y = np.array([-side_length/2.,-side_length/2.,side_length/2.,side_length/2.])
xBounds = x*np.cos(np.deg2rad(30.)) - y*np.sin(np.deg2rad(30.))
yBounds = x*np.sin(np.deg2rad(30.)) + y*np.cos(np.deg2rad(30.))

xBounds = np.append(xBounds,xBounds[0])
yBounds = np.append(yBounds,yBounds[0])

ax9.plot(xBounds,yBounds,'-',color='C0',alpha=0.5)

ax7.set_ylim(-3300.,3300.)
ax8.set_ylim(-3300.,3300.)
ax9.set_ylim(-3300.,3300.)

ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
ax7.spines['left'].set_visible(False)
ax8.spines['right'].set_visible(False)
ax8.spines['top'].set_visible(False)
ax8.spines['bottom'].set_visible(False)
ax8.spines['left'].set_visible(False)
ax9.spines['right'].set_visible(False)
ax9.spines['top'].set_visible(False)
ax9.spines['bottom'].set_visible(False)
ax9.spines['left'].set_visible(False)

ax7.set_xticks(())
ax8.set_xticks(())
ax9.set_xticks(())
ax7.set_yticks(())
ax8.set_yticks(())
ax9.set_yticks(())


ax4.set_title('Amalia',fontsize=8,rotation=0.)
ax5.set_title('circle',fontsize=8,rotation=0.)
ax6.set_title('square',fontsize=8,rotation=0.)

ax8.set_xlabel('boundary shape',fontsize=8)




min4 = np.min(np.array([np.min(f_param4),np.min(f_direct4),np.min(f_grid4)]))
min6 = np.min(np.array([np.min(f_param6),np.min(f_direct6),np.min(f_grid6)]))
min8 = np.min(np.array([np.min(f_param8),np.min(f_direct8),np.min(f_grid8)]))

max4 = np.max(np.array([np.max(f_param4),np.max(f_direct4),np.max(f_grid4)]))
max6 = np.max(np.array([np.max(f_param6),np.max(f_direct6),np.max(f_grid6)]))
max8 = np.max(np.array([np.max(f_param8),np.max(f_direct8),np.max(f_grid8)]))

b4 = np.logspace(np.log10(min4),np.log10(max4),25)
b6 = np.logspace(np.log10(min6),np.log10(max6),25)
b8 = np.logspace(np.log10(min8),np.log10(max8),25)

num = 100

ax1.hist(f_grid4[0:num],bins=b4,color='C1',alpha=0.5)
ax2.hist(f_grid6[0:num],bins=b6,color='C1',alpha=0.5)
ax3.hist(f_grid8[0:num],bins=b8,color='C1',alpha=0.5)

ax1.hist(f_param4[0:num],bins=b4,color='C0',alpha=0.5)
ax2.hist(f_param6[0:num],bins=b6,color='C0',alpha=0.5)
ax3.hist(f_param8[0:num],bins=b8,color='C0',alpha=0.5)

ax1.hist(f_direct4[0:num],bins=b4,color='C2',alpha=0.5)
ax2.hist(f_direct6[0:num],bins=b6,color='C2',alpha=0.5)
ax3.hist(f_direct8[0:num],bins=b8,color='C2',alpha=0.5)

ax1.hist(f_param4[0:num],bins=b4,color='C0',alpha=0.5)
ax2.hist(f_param6[0:num],bins=b6,color='C0',alpha=0.5)
ax3.hist(f_param8[0:num],bins=b8,color='C0',alpha=0.5)

min4 = np.min(np.array([np.min(AEP_param4),np.min(AEP_direct4),np.min(AEP_grid4)]))
min6 = np.min(np.array([np.min(AEP_param6),np.min(AEP_direct6),np.min(AEP_grid6)]))
min8 = np.min(np.array([np.min(AEP_param8),np.min(AEP_direct8),np.min(AEP_grid8)]))

max4 = np.max(np.array([np.max(AEP_param4),np.max(AEP_direct4),np.max(AEP_grid4)]))
max6 = np.max(np.array([np.max(AEP_param6),np.max(AEP_direct6),np.max(AEP_grid6)]))
max8 = np.max(np.array([np.max(AEP_param8),np.max(AEP_direct8),np.max(AEP_grid8)]))

b4 = np.linspace(min4,max4,25)
b6 = np.linspace(min6,max6,25)
b8 = np.linspace(min8,max8,25)

b6 = np.linspace(min4,max4,25)
b8 = np.linspace(min4,max4,25)


ax4.hist(AEP_grid4[0:num],bins=b4,color='C1',alpha=0.5,label='grid')
ax5.hist(AEP_grid6[0:num],bins=b6,color='C1',alpha=0.5)
ax6.hist(AEP_grid8[0:num],bins=b8,color='C1',alpha=0.5)
ax4.hist(AEP_param4[0:num],bins=b4,color='C0',alpha=0.5,label='parameterized')
ax5.hist(AEP_param6[0:num],bins=b6,color='C0',alpha=0.5)
ax6.hist(AEP_param8[0:num],bins=b8,color='C0',alpha=0.5)
ax4.hist(AEP_direct4[0:num],bins=b4,color='C2',alpha=0.5,label='direct')
ax5.hist(AEP_direct6[0:num],bins=b6,color='C2',alpha=0.5)
ax6.hist(AEP_direct8[0:num],bins=b8,color='C2',alpha=0.5)
ax4.hist(AEP_param4[0:num],bins=b4,color='C0',alpha=0.5)
ax5.hist(AEP_param6[0:num],bins=b6,color='C0',alpha=0.5)
ax6.hist(AEP_param8[0:num],bins=b8,color='C0',alpha=0.5)

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


# ax4.plot(windDirections, windFrequencies, '-', color=tuscan, linewidth=3,zorder=0)
# ax4.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])

# ax4.set_rgrids([2,4,6], angle=-35.)
# ax4.set_yticklabels(['2 m/s','4 m/s','6 m/s'],horizontalalignment='center')

ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')

ticks = np.array([10.,100.,1000.,5000.])
ax1.set_xticks(ticks)
ax1.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))
ax2.set_xticks(ticks)
ax2.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))
ax3.set_xticks(ticks)
ax3.set_xticklabels((r'10$^1$',r'10$^2$',r'10$^3$',''))
ax1.set_xlim(10.,5000.)

ticks = np.array([0.9,0.95,1.0])
ax4.set_xticks(ticks)
ax4.set_xticklabels(('0.9','0.95','1.0'))
ax5.set_xticks(ticks)
ax5.set_xticklabels(('0.9','0.95','1.0'))
ax6.set_xticks(ticks)
ax6.set_xticklabels(('0.9','0.95','1.0'))

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


plt.subplots_adjust(top = 0.94, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.6, wspace = 0.2)


ax4.text(0.907,41.,'a',fontsize=11,color='black')
ax5.text(0.907,41.,'b',fontsize=11,color='black')
ax6.text(0.907,41.,'c',fontsize=11,color='black')
ax1.text(12.5,35.,'d',fontsize=11,color='black')
ax2.text(12.5,35.,'e',fontsize=11,color='black')
ax3.text(12.5,35.,'f',fontsize=11,color='black')

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

# plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/results_bounds.pdf',transparent=True)
plt.show()
