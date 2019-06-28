import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# from windRoses import *

"""colors"""
robin = '#B2DBD5'
nectar = '#EAB364'
tuscan = '#B2473E'
olive = '#ACBD78'

robin = 'C0'
nectar = 'C1'
tuscan = 'C1'
royal = 'C3'

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10',size=8)


# folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_4.0_victorvilleRose'
# folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_4.0_victorvilleRose'
# folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_4.0_victorvilleRose'

# with open('%s/AEP.txt'%folder_grid) as my_file:
#     AEP_grid = my_file.readlines()
# for p in range(len(AEP_grid)):
#     AEP_grid[p] = float(AEP_grid[p])
#
# with open('%s/AEP.txt'%folder_direct) as my_file:
#     AEP_direct = my_file.readlines()
# for p in range(len(AEP_direct)):
#     AEP_direct[p] = float(AEP_direct[p])
#
# with open('%s/AEP.txt'%folder_param) as my_file:
#     AEP_param = my_file.readlines()
# for p in range(len(AEP_param)):
#     AEP_param[p] = float(AEP_param[p])
#
#
# print 'direct: ',  np.max(AEP_direct[0:100])
# print 'param: ',  np.max(AEP_param[0:100])
# print 'grid: ',  np.max(AEP_grid[0:100])


d_base = 889885.834065
p_base = 886639.398754
g_base = 873169.69392
b_base = 833453.470716884
base = np.array([b_base,g_base,p_base,d_base])

print (d_base-b_base)/b_base
print (p_base-b_base)/b_base
print (g_base-b_base)/b_base

d_6 = 960932.660336
p_6 = 960476.8279
g_6 = 954136.187431
b_6 = 934698.2333983908
six = np.array([b_6,g_6,p_6,d_6])

print (d_6-b_6)/b_6
print (p_6-b_6)/b_6
print (g_6-b_6)/b_6

d_8 = 994438.529361
p_8 = 994902.113881
g_8 = 992642.870223
b_8 = 973337.4197978062
eight = np.array([b_8,g_8,p_8,d_8])

print (d_8-b_8)/b_8
print (p_8-b_8)/b_8
print (g_8-b_8)/b_8

d_circle = 887909.788146
p_circle = 886315.754861
g_circle = 882782.998293
b_circle = 825614.6393008
circle = np.array([b_circle,g_circle,p_circle,d_circle])

print (d_circle-b_circle)/b_circle
print (p_circle-b_circle)/b_circle
print (g_circle-b_circle)/b_circle

d_square = 892719.886631
p_square = 892883.370424
g_square = 886954.157891
b_square = 830776.6255677297
square = np.array([b_square,g_square,p_square,d_square])

print (d_square-b_square)/b_square
print (p_square-b_square)/b_square
print (g_square-b_square)/b_square

d_ukiah = 963872.640356
p_ukiah = 959998.056736
g_ukiah = 927694.932959
b_ukiah = 913115.0771472208
ukiah = np.array([b_ukiah,g_ukiah,p_ukiah,d_ukiah])

print (d_ukiah-b_ukiah)/b_ukiah
print (p_ukiah-b_ukiah)/b_ukiah
print (g_ukiah-b_ukiah)/b_ukiah

d_victorville = 1107248.95024
p_victorville = 1105525.8565
g_victorville = 1089947.64585
b_victorville = 1047172.3399274771
victorville = np.array([b_victorville,g_victorville,p_victorville,d_victorville])

print (d_victorville-b_victorville)/b_victorville
print (p_victorville-b_victorville)/b_victorville
print (g_victorville-b_victorville)/b_victorville

fig = plt.figure(figsize=[6.5,3.])
ax1 = plt.subplot2grid((1, 3), (0, 0))
ax2 = plt.subplot2grid((1, 3), (0, 1))
ax3 = plt.subplot2grid((1, 3), (0, 2))

size = 5

ax1.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color=royal,markersize=size,label='direct')
ax1.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color=robin,markersize=size,label='parameterized')
ax1.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color=tuscan,markersize=size,label='grid')

# ax1.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color=tuscan,markersize=size)
# ax1.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color=royal,markersize=size)
# ax1.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color=robin,markersize=size)

ax1.plot(2,100.*(six[1]-six[0])/six[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax1.plot(2,100.*(six[3]-six[0])/six[0],'o',alpha=0.5,color=royal,markersize=size)
ax1.plot(2,100.*(six[2]-six[0])/six[0],'o',alpha=0.5,color=robin,markersize=size)

ax1.plot(3,100.*(eight[1]-eight[0])/eight[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax1.plot(3,100.*(eight[3]-eight[0])/eight[0],'o',alpha=0.5,color=royal,markersize=size)
ax1.plot(3,100.*(eight[2]-eight[0])/eight[0],'o',alpha=0.5,color=robin,markersize=size)


ax2.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax2.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color=royal,markersize=size)
ax2.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color=robin,markersize=size)

ax2.plot(2,100.*(ukiah[1]-ukiah[0])/ukiah[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax2.plot(2,100.*(ukiah[3]-ukiah[0])/ukiah[0],'o',alpha=0.5,color=royal,markersize=size)
ax2.plot(2,100.*(ukiah[2]-ukiah[0])/ukiah[0],'o',alpha=0.5,color=robin,markersize=size)

ax2.plot(3,100.*(victorville[1]-victorville[0])/victorville[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax2.plot(3,100.*(victorville[3]-victorville[0])/victorville[0],'o',alpha=0.5,color=royal,markersize=size)
ax2.plot(3,100.*(victorville[2]-victorville[0])/victorville[0],'o',alpha=0.5,color=robin,markersize=size)

ax1.legend(loc=1,prop={'size': 8})


ax3.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax3.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color=royal,markersize=size)
ax3.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color=robin,markersize=size)

ax3.plot(2,100.*(circle[1]-circle[0])/circle[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax3.plot(2,100.*(circle[3]-circle[0])/circle[0],'o',alpha=0.5,color=royal,markersize=size)
ax3.plot(2,100.*(circle[2]-circle[0])/circle[0],'o',alpha=0.5,color=robin,markersize=size)

ax3.plot(3,100.*(square[1]-square[0])/square[0],'o',alpha=0.5,color=tuscan,markersize=size)
ax3.plot(3,100.*(square[3]-square[0])/square[0],'o',alpha=0.5,color=royal,markersize=size)
ax3.plot(3,100.*(square[2]-square[0])/square[0],'o',alpha=0.5,color=robin,markersize=size)

ax1.set_yticks((0.,2.,4.,6.,8.))
ax2.set_yticks((0.,2.,4.,6.,8.))
ax3.set_yticks((0.,2.,4.,6.,8.))

ax1.set_yticklabels(('0','2','4','6','8'))
ax2.set_yticklabels(('','','','',''))
ax3.set_yticklabels(('','','','',''))

ax1.set_ylim(0.,8.)
ax2.set_ylim(0.,8.)
ax3.set_ylim(0.,8.)

ax1.set_xlim(0.8,3.2)
ax2.set_xlim(0.8,3.2)
ax3.set_xlim(0.8,3.2)

ax1.set_xticks((1,2,3))
ax2.set_xticks((1,2,3))
ax3.set_xticks((1,2,3))

ax1.set_xticklabels(('4 D','6 D','8 D'),fontsize=8)
ax2.set_xticklabels(('North\nIsland','Ukiah','Victorville'),fontsize=8,rotation=45.)
ax3.set_xticklabels(('Amalia','circle','square'),fontsize=8,rotation=45.)

ax1.set_ylabel('% increase over baseline',fontsize=8)

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

plt.subplots_adjust(top = 0.93, bottom = 0.4, right = 0.98, left = 0.1,
            hspace = 0.6, wspace = 0.2)

ax1.text(1.,8.4,'a',fontsize=11,color='black')
ax2.text(1.,8.4,'b',fontsize=11,color='black')
ax3.text(1.,8.4,'c',fontsize=11,color='black')

# ax1.set_xlabel('turbine spacing',fontsize=8)
# ax2.set_xlabel('wind rose',fontsize=8)
# ax3.set_xlabel('boundary shape',fontsize=8)

ax1.text(2.,-3.,'turbine spacing',fontsize=8,horizontalalignment='center')
ax2.text(2.,-3.,'wind rose',fontsize=8,horizontalalignment='center')
ax3.text(2.,-3.,'boundary shape',fontsize=8,horizontalalignment='center')

#
plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/results_max.pdf',transparent=True)
plt.show()
