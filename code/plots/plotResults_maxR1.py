import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/')
from position_constraints import *
import scipy as sp
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
# from windRoses import *

import matplotlib as mpl
mpl.rc('font', family = 'serif', serif = 'cmr10',size=8)

# bound = 'amalia'
# spac = '4.0'
# rose = 'victorvilleRose'
# folder_direct = '/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/code/plots/revised_data/opt-results-revision/direct/%s_%s_%s'%(bound,spac,rose)
# folder_param = '/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/code/plots/revised_data/opt-results-revision/param/%s_%s_%s'%(bound,spac,rose)
# folder_grid = '/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/code/plots/revised_data/opt-results-revision/grid/%s_%s_%s_my1'%(bound,spac,rose)
# #
# with open('%s/AEPopt.txt'%folder_grid) as my_file:
#     AEP_grid = my_file.readlines()
# for p in range(len(AEP_grid)):
#     AEP_grid[p] = float(AEP_grid[p])
#
# with open('%s/AEPopt.txt'%folder_direct) as my_file:
#     AEP_direct = my_file.readlines()
# for p in range(len(AEP_direct)):
#     AEP_direct[p] = float(AEP_direct[p])
#
# with open('%s/AEPopt.txt'%folder_param) as my_file:
#     AEP_param = my_file.readlines()
# for p in range(len(AEP_param)):
#     AEP_param[p] = float(AEP_param[p])


# print 'direct: ',  np.max(AEP_direct[0:100])
# print 'param: ',  np.max(AEP_param[0:100])
# print 'grid: ',  np.max(AEP_grid[0:100])

# #opt base
# # base 745490.812401
# # U 609826.62169
# # V 1026981.87294
# # circle 829796.857156
# # square 845120.195615
# # 6 944478.512624
# # 8 1059055.59708
#
# d_base = 1319337.4624830007
# p_base = 1321848.6760216402
# g_base = 1285013.9936945664
# b_base = 1272384.14711
# base = np.array([b_base,g_base,p_base,d_base])
#
# # print (d_base-b_base)/b_base
# # print (p_base-b_base)/b_base
# # print (g_base-b_base)/b_base
#
# d_6 = 1444130.5355461645
# p_6 = 1447939.702717834
# g_6 = 1428351.7311277087
# b_6 = 1409997.66677
# six = np.array([b_6,g_6,p_6,d_6])
#
# # print (d_6-b_6)/b_6
# # print (p_6-b_6)/b_6
# # print (g_6-b_6)/b_6
#
# d_8 = 1504627.0471007032
# p_8 = 1507076.2508215704
# g_8 = 1495577.354103449
# b_8 = 1478346.00892
# eight = np.array([b_8,g_8,p_8,d_8])
#
# # print (d_8-b_8)/b_8
# # print (p_8-b_8)/b_8
# # print (g_8-b_8)/b_8
#
# d_circle = 1315815.5893947128
# p_circle = 1321429.8945051455
# # g_circle = 1259051.6321162586
# g_circle = 1288343.196518258
# b_circle = 1262718.1614
# circle = np.array([b_circle,g_circle,p_circle,d_circle])
#
# # print (d_circle-b_circle)/b_circle
# # print (p_circle-b_circle)/b_circle
# # print (g_circle-b_circle)/b_circle
#
# d_square = 1329733.8046032637
# p_square = 1331340.662503046
# g_square = 1302801.6310573649
# b_square = 1269726.70336
# square = np.array([b_square,g_square,p_square,d_square])
#
# # print (d_square-b_square)/b_square
# # print (p_square-b_square)/b_square
# # print (g_square-b_square)/b_square
#
# d_ukiah = 1415278.9473485365
# p_ukiah = 1417439.4876242485
# g_ukiah = 1381832.8806432036
# b_ukiah = 1373083.98385
# ukiah = np.array([b_ukiah,g_ukiah,p_ukiah,d_ukiah])
#
# # print (d_ukiah-b_ukiah)/b_ukiah
# # print (p_ukiah-b_ukiah)/b_ukiah
# # print (g_ukiah-b_ukiah)/b_ukiah
#
# d_victorville = 1497711.1393903312
# p_victorville = 1498899.2947109619
# g_victorville = 1462698.469289606
# b_victorville = 1451109.14502
# victorville = np.array([b_victorville,g_victorville,p_victorville,d_victorville])
#
# # print (d_victorville-b_victorville)/b_victorville
# # print (p_victorville-b_victorville)/b_victorville
# # print (g_victorville-b_victorville)/b_victorville





#true base
# U 1425429.79046
# V 1520007.41548
# 6 1446206.6559
# 8 1506032.52603
# c 1301035.62243
# s 1307208.73706

d_base = 1357709.4021
p_base = 1360629.14391
g_base = 1326187.41196
b_base = 1311195.09809
base = np.array([b_base,g_base,p_base,d_base])

print (d_base-g_base)/d_base

# print (d_base-b_base)/b_base
# print (p_base-b_base)/b_base
# print (g_base-b_base)/b_base

d_6 = 1470505.76598
p_6 = 1475519.11985
g_6 = 1456771.55278
b_6 = 1446206.6559
six = np.array([b_6,g_6,p_6,d_6])

print (d_6-g_6)/d_6

# print (d_6-b_6)/b_6
# print (p_6-b_6)/b_6
# print (g_6-b_6)/b_6

d_8 = 1519337.74429
p_8 = 1524358.24114
g_8 = 1513847.77186
b_8 = 1506032.52603
eight = np.array([b_8,g_8,p_8,d_8])

print (d_8-g_8)/d_8

# print (d_8-b_8)/b_8
# print (p_8-b_8)/b_8
# print (g_8-b_8)/b_8

d_circle = 1355977.26239
p_circle = 1360955.2908
# g_circle = 1300385.34741
g_circle = 1329373.22315
b_circle = 1301035.62243
circle = np.array([b_circle,g_circle,p_circle,d_circle])

print (d_circle-g_circle)/d_circle

# print (d_circle-b_circle)/b_circle
# print (p_circle-b_circle)/b_circle
# print (g_circle-b_circle)/b_circle

d_square = 1366805.77434
p_square = 1370147.21542
g_square = 1346112.53376
b_square = 1307208.73706
square = np.array([b_square,g_square,p_square,d_square])

print (d_square-g_square)/d_square

# print (d_square-b_square)/b_square
# print (p_square-b_square)/b_square
# print (g_square-b_square)/b_square

d_ukiah = 1470807.90913
p_ukiah = 1474659.63653
g_ukiah = 1442684.54864
b_ukiah = 1425429.79046
ukiah = np.array([b_ukiah,g_ukiah,p_ukiah,d_ukiah])

print (d_ukiah-g_ukiah)/d_ukiah

# print (d_ukiah-b_ukiah)/b_ukiah
# print (p_ukiah-b_ukiah)/b_ukiah
# print (g_ukiah-b_ukiah)/b_ukiah

d_victorville = 1561992.89748
p_victorville = 1564138.02158
# g_victorville = 1527258.98838
g_victorville = 1529512.52545
b_victorville = 1520007.41548
victorville = np.array([b_victorville,g_victorville,p_victorville,d_victorville])

print (d_victorville-g_victorville)/d_victorville

# print (d_victorville-b_victorville)/b_victorville
# print (p_victorville-b_victorville)/b_victorville
# print (g_victorville-b_victorville)/b_victorville


fig = plt.figure(figsize=[6.5,2.5])
ax1 = plt.subplot2grid((1, 3), (0, 0))
ax2 = plt.subplot2grid((1, 3), (0, 1))
ax3 = plt.subplot2grid((1, 3), (0, 2))

size = 5

ax1.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color='C2',markersize=size,label='direct')
ax1.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color='C0',markersize=size,label='parameterized')
ax1.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color='C1',markersize=size,label='grid')

# ax1.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color='C1',markersize=size)
# ax1.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color='C2',markersize=size)
# ax1.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color='C0',markersize=size)

ax1.plot(2,100.*(six[1]-six[0])/six[0],'o',alpha=0.5,color='C1',markersize=size)
ax1.plot(2,100.*(six[3]-six[0])/six[0],'o',alpha=0.5,color='C2',markersize=size)
ax1.plot(2,100.*(six[2]-six[0])/six[0],'o',alpha=0.5,color='C0',markersize=size)

ax1.plot(3,100.*(eight[1]-eight[0])/eight[0],'o',alpha=0.5,color='C1',markersize=size)
ax1.plot(3,100.*(eight[3]-eight[0])/eight[0],'o',alpha=0.5,color='C2',markersize=size)
ax1.plot(3,100.*(eight[2]-eight[0])/eight[0],'o',alpha=0.5,color='C0',markersize=size)


ax2.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color='C1',markersize=size)
ax2.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color='C2',markersize=size)
ax2.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color='C0',markersize=size)

ax2.plot(2,100.*(ukiah[1]-ukiah[0])/ukiah[0],'o',alpha=0.5,color='C1',markersize=size)
ax2.plot(2,100.*(ukiah[3]-ukiah[0])/ukiah[0],'o',alpha=0.5,color='C2',markersize=size)
ax2.plot(2,100.*(ukiah[2]-ukiah[0])/ukiah[0],'o',alpha=0.5,color='C0',markersize=size)

ax2.plot(3,100.*(victorville[1]-victorville[0])/victorville[0],'o',alpha=0.5,color='C1',markersize=size)
ax2.plot(3,100.*(victorville[3]-victorville[0])/victorville[0],'o',alpha=0.5,color='C2',markersize=size)
ax2.plot(3,100.*(victorville[2]-victorville[0])/victorville[0],'o',alpha=0.5,color='C0',markersize=size)

ax1.legend(loc=1,prop={'size': 8})


ax3.plot(1,100.*(base[1]-base[0])/base[0],'o',alpha=0.5,color='C1',markersize=size)
ax3.plot(1,100.*(base[3]-base[0])/base[0],'o',alpha=0.5,color='C2',markersize=size)
ax3.plot(1,100.*(base[2]-base[0])/base[0],'o',alpha=0.5,color='C0',markersize=size)

ax3.plot(2,100.*(circle[1]-circle[0])/circle[0],'o',alpha=0.5,color='C1',markersize=size)
ax3.plot(2,100.*(circle[3]-circle[0])/circle[0],'o',alpha=0.5,color='C2',markersize=size)
ax3.plot(2,100.*(circle[2]-circle[0])/circle[0],'o',alpha=0.5,color='C0',markersize=size)

ax3.plot(3,100.*(square[1]-square[0])/square[0],'o',alpha=0.5,color='C1',markersize=size)
ax3.plot(3,100.*(square[3]-square[0])/square[0],'o',alpha=0.5,color='C2',markersize=size)
ax3.plot(3,100.*(square[2]-square[0])/square[0],'o',alpha=0.5,color='C0',markersize=size)

ax1.set_yticks((0.,2.,4.,6.,8.))
ax2.set_yticks((0.,2.,4.,6.,8.))
ax3.set_yticks((0.,2.,4.,6.,8.))

ax1.set_yticklabels(('0','2','4','6'))
ax2.set_yticklabels(('','','',''))
ax3.set_yticklabels(('','','',''))

ax1.set_ylim(0.,6.)
ax2.set_ylim(0.,6.)
ax3.set_ylim(0.,6.)
#
# ax1.set_xlim(0.8,3.2)
# ax2.set_xlim(0.8,3.2)
# ax3.set_xlim(0.8,3.2)

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

plt.subplots_adjust(top = 0.9, bottom = 0.2, right = 0.98, left = 0.1,
            hspace = 0.6, wspace = 0.2)

ax1.text(1.,6.4,'a',fontsize=11,color='black')
ax2.text(1.,6.4,'b',fontsize=11,color='black')
ax3.text(1.,6.4,'c',fontsize=11,color='black')

# ax1.set_xlabel('turbine spacing',fontsize=8)
# ax2.set_xlabel('wind rose',fontsize=8)
# ax3.set_xlabel('boundary shape',fontsize=8)

ax1.text(2.,-3.,'turbine spacing',fontsize=8,horizontalalignment='center')
ax2.text(2.,-3.,'wind rose',fontsize=8,horizontalalignment='center')
ax3.text(2.,-3.,'boundary shape',fontsize=8,horizontalalignment='center')

#
# plt.savefig('/Users/ningrsrch/Dropbox/Projects/stanley2019-variable-reduction/paper/paper-figures/results_maxR1.pdf',transparent=True)
plt.show()
