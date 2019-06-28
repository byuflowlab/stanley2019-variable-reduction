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





folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_6.0_northIslandRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_6.0_northIslandRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_6.0_northIslandRose'

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


folder_direct = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/direct/amalia_8.0_northIslandRose'
folder_param = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/param/amalia_8.0_northIslandRose'
folder_grid = '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients/opt-results/grid/amalia_8.0_northIslandRose'


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



ax9 = plt.subplot2grid((3, 3), (0, 2))
ax7 = plt.subplot2grid((3, 3), (0, 0),sharey=ax9,sharex=ax9)
ax8 = plt.subplot2grid((3, 3), (0, 1),sharey=ax9,sharex=ax9)



# ax7.set_xlim(-4800.,4800.)
# ax8.set_xlim(-4800.,4800.)
# ax9.set_xlim(-4800.,4800.)
ax7.axis('equal')
ax8.axis('equal')
ax9.axis('equal')


# s = 4.
# x = np.array([-2.*s,-s,0.,s,2.*s])
# y = np.array([-2.*s,-s,0.,s,2.*s])
# X,Y = np.meshgrid(x,y)
# X = np.ndarray.flatten(X)
# Y = np.ndarray.flatten(Y)
# for i in range(len(X)):
#     circ = plt.Circle((X[i],Y[i]),0.5,facecolor=robin,edgecolor=robin,alpha=0.5)
#     ax7.add_patch(circ)
#
# s = 6.
# x = np.array([-2.*s,-s,0.,s,2.*s])
# y = np.array([-2.*s,-s,0.,s,2.*s])
# X,Y = np.meshgrid(x,y)
# X = np.ndarray.flatten(X)
# Y = np.ndarray.flatten(Y)
# for i in range(len(X)):
#     circ = plt.Circle((X[i],Y[i]),0.5,facecolor=robin,edgecolor=robin,alpha=0.5)
#     ax8.add_patch(circ)
#
# s = 8.
# x = np.array([-2.*s,-s,0.,s,2.*s])
# y = np.array([-2.*s,-s,0.,s,2.*s])
# X,Y = np.meshgrid(x,y)
# X = np.ndarray.flatten(X)
# Y = np.ndarray.flatten(Y)
# for i in range(len(X)):
#     circ = plt.Circle((X[i],Y[i]),0.5,facecolor=robin,edgecolor=robin,alpha=0.5)
#     ax9.add_patch(circ)
#
# ax7.set_ylim(-19,19)
# ax8.set_ylim(-19,19)
# ax9.set_ylim(-19,19)

nTurbs = 100
rotorDiameter = np.ones(nTurbs)*130.

spacing = 4.
side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
a = side_length**2
circle_radius = np.sqrt(a/np.pi)

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

bx4 = np.append(xBounds,xBounds[0])
by4 = np.append(yBounds,yBounds[0])

x4 = np.array([-1055.45573535,  -582.10981314,  -108.76389094,   364.58203126,
       -1356.51783135,  -883.17190915,  -409.82598695,    63.51993525,
         536.86585745,  1010.21177966, -1657.57992736, -1184.23400516,
        -710.88808296,  -237.54216076,   235.80376144,   709.14968365,
        1182.49560585,  1655.84152805, -1958.64202337, -1485.29610117,
       -1011.95017897,  -538.60425677,   -65.25833456,   408.08758764,
         881.43350984,  1354.77943204,  1828.12535424, -1786.35819718,
       -1313.01227498,  -839.66635278,  -366.32043057,   107.02549163,
         580.37141383,  1053.71733603,  1527.06325823,  2000.40918044,
       -2087.42029319, -1614.07437099, -1140.72844879,  -667.38252658,
        -194.03660438,   279.30931782,   752.65524002,  1226.00116223,
        1699.34708443,  2172.69300663, -1915.136467  , -1441.79054479,
        -968.44462259,  -495.09870039,   -21.75277819,   451.59314401,
         924.93906622,  1398.28498842,  1871.63091062, -2216.19856301,
       -1742.8526408 , -1269.5067186 ,  -796.1607964 ,  -322.8148742 ,
         150.531048  ,   623.87697021,  1097.22289241,  1570.56881461,
        2043.91473681, -2043.91473681, -1570.56881461, -1097.22289241,
        -623.87697021,  -150.531048  ,   322.8148742 ,   796.1607964 ,
        1269.5067186 ,  1742.8526408 ,  2216.19856301, -1871.63091062,
       -1398.28498842,  -924.93906622,  -451.59314401,    21.75277819,
         495.09870039,   968.44462259,  1441.79054479,  1915.136467  ,
       -2172.69300663, -1699.34708443, -1226.00116223,  -752.65524002,
        -279.30931782,   194.03660438,   667.38252658,  1140.72844879,
       -2000.40918044, -1527.06325823, -1053.71733603,  -580.37141383,
        -107.02549163,   366.32043057,  -881.43350984,  -408.08758764])
y4 = np.array([-2.84007553e+03, -2.84007553e+03, -2.84007553e+03, -2.84007553e+03,
       -2.36672961e+03, -2.36672961e+03, -2.36672961e+03, -2.36672961e+03,
       -2.36672961e+03, -2.36672961e+03, -1.89338369e+03, -1.89338369e+03,
       -1.89338369e+03, -1.89338369e+03, -1.89338369e+03, -1.89338369e+03,
       -1.89338369e+03, -1.89338369e+03, -1.42003777e+03, -1.42003777e+03,
       -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -1.42003777e+03,
       -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -9.46691844e+02,
       -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
       -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
       -4.73345922e+02, -4.73345922e+02, -4.73345922e+02, -4.73345922e+02,
       -4.73345922e+02, -4.73345922e+02, -4.73345922e+02, -4.73345922e+02,
       -4.73345922e+02, -4.73345922e+02, -2.27373675e-13, -2.27373675e-13,
       -2.27373675e-13, -2.27373675e-13, -2.27373675e-13, -2.27373675e-13,
       -2.27373675e-13, -2.27373675e-13, -2.27373675e-13,  4.73345922e+02,
        4.73345922e+02,  4.73345922e+02,  4.73345922e+02,  4.73345922e+02,
        4.73345922e+02,  4.73345922e+02,  4.73345922e+02,  4.73345922e+02,
        4.73345922e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
        9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
        9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  1.42003777e+03,
        1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
        1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
        1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
        1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
        2.36672961e+03,  2.36672961e+03,  2.36672961e+03,  2.36672961e+03,
        2.36672961e+03,  2.36672961e+03,  2.84007553e+03,  2.84007553e+03])

spacing = 6.
side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
a = side_length**2
circle_radius = np.sqrt(a/np.pi)

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

bx6 = np.append(xBounds,xBounds[0])
by6 = np.append(yBounds,yBounds[0])

x6 = np.array([-1583.18360302,  -873.16471971,  -163.14583641,   546.87304689,
       -2034.77674703, -1324.75786373,  -614.73898043,    95.27990288,
         805.29878618,  1515.31766948, -2486.36989105, -1776.35100774,
       -1066.33212444,  -356.31324114,   353.70564217,  1063.72452547,
        1773.74340877,  2483.76229208, -2937.96303506, -2227.94415176,
       -1517.92526845,  -807.90638515,   -97.88750185,   612.13138146,
        1322.15026476,  2032.16914806,  2742.18803137, -2679.53729577,
       -1969.51841247, -1259.49952916,  -549.48064586,   160.53823744,
         870.55712075,  1580.57600405,  2290.59488735,  3000.61377065,
       -3131.13043978, -2421.11155648, -1711.09267318, -1001.07378987,
        -291.05490657,   418.96397673,  1128.98286003,  1839.00174334,
        2549.02062664,  3259.03950994, -2872.70470049, -2162.68581719,
       -1452.66693389,  -742.64805059,   -32.62916728,   677.38971602,
        1387.40859932,  2097.42748263,  2807.44636593, -3324.29784451,
       -2614.27896121, -1904.2600779 , -1194.2411946 ,  -484.2223113 ,
         225.79657201,   935.81545531,  1645.83433861,  2355.85322192,
        3065.87210522, -3065.87210522, -2355.85322192, -1645.83433861,
        -935.81545531,  -225.79657201,   484.2223113 ,  1194.2411946 ,
        1904.2600779 ,  2614.27896121,  3324.29784451, -2807.44636593,
       -2097.42748263, -1387.40859932,  -677.38971602,    32.62916728,
         742.64805059,  1452.66693389,  2162.68581719,  2872.70470049,
       -3259.03950994, -2549.02062664, -1839.00174334, -1128.98286003,
        -418.96397673,   291.05490657,  1001.07378987,  1711.09267318,
       -3000.61377065, -2290.59488735, -1580.57600405,  -870.55712075,
        -160.53823744,   549.48064586, -1322.15026476,  -612.13138146])
y6 = np.array([-4.26011330e+03, -4.26011330e+03, -4.26011330e+03, -4.26011330e+03,
       -3.55009442e+03, -3.55009442e+03, -3.55009442e+03, -3.55009442e+03,
       -3.55009442e+03, -3.55009442e+03, -2.84007553e+03, -2.84007553e+03,
       -2.84007553e+03, -2.84007553e+03, -2.84007553e+03, -2.84007553e+03,
       -2.84007553e+03, -2.84007553e+03, -2.13005665e+03, -2.13005665e+03,
       -2.13005665e+03, -2.13005665e+03, -2.13005665e+03, -2.13005665e+03,
       -2.13005665e+03, -2.13005665e+03, -2.13005665e+03, -1.42003777e+03,
       -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -1.42003777e+03,
       -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -1.42003777e+03,
       -7.10018883e+02, -7.10018883e+02, -7.10018883e+02, -7.10018883e+02,
       -7.10018883e+02, -7.10018883e+02, -7.10018883e+02, -7.10018883e+02,
       -7.10018883e+02, -7.10018883e+02,  6.82121026e-13,  6.82121026e-13,
        6.82121026e-13,  6.82121026e-13,  6.82121026e-13,  6.82121026e-13,
        6.82121026e-13,  6.82121026e-13,  6.82121026e-13,  7.10018883e+02,
        7.10018883e+02,  7.10018883e+02,  7.10018883e+02,  7.10018883e+02,
        7.10018883e+02,  7.10018883e+02,  7.10018883e+02,  7.10018883e+02,
        7.10018883e+02,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
        1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
        1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  2.13005665e+03,
        2.13005665e+03,  2.13005665e+03,  2.13005665e+03,  2.13005665e+03,
        2.13005665e+03,  2.13005665e+03,  2.13005665e+03,  2.13005665e+03,
        2.84007553e+03,  2.84007553e+03,  2.84007553e+03,  2.84007553e+03,
        2.84007553e+03,  2.84007553e+03,  2.84007553e+03,  2.84007553e+03,
        3.55009442e+03,  3.55009442e+03,  3.55009442e+03,  3.55009442e+03,
        3.55009442e+03,  3.55009442e+03,  4.26011330e+03,  4.26011330e+03])

spacing = 8.
side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
a = side_length**2
circle_radius = np.sqrt(a/np.pi)

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

bx8 = np.append(xBounds,xBounds[0])
by8 = np.append(yBounds,yBounds[0])

x8 = np.array([-2110.91147069, -1164.21962629,  -217.52778188,   729.16406252,
       -2713.03566271, -1766.3438183 ,  -819.6519739 ,   127.0398705 ,
        1073.73171491,  2020.42355931, -3315.15985473, -2368.46801032,
       -1421.77616592,  -475.08432152,   471.60752289,  1418.29936729,
        2364.9912117 ,  3311.6830561 , -3917.28404675, -2970.59220234,
       -2023.90035794, -1077.20851353,  -130.51666913,   816.17517527,
        1762.86701968,  2709.55886408,  3656.25070849, -3572.71639436,
       -2626.02454996, -1679.33270555,  -732.64086115,   214.05098326,
        1160.74282766,  2107.43467206,  3054.12651647,  4000.81836087,
       -4174.84058638, -3228.14874197, -2281.45689757, -1334.76505317,
        -388.07320876,   558.61863564,  1505.31048005,  2452.00232445,
        3398.69416885,  4345.38601326, -3830.27293399, -2883.58108959,
       -1936.88924518,  -990.19740078,   -43.50555638,   903.18628803,
        1849.87813243,  2796.56997684,  3743.26182124, -4432.39712601,
       -3485.70528161, -2539.0134372 , -1592.3215928 ,  -645.62974839,
         301.06209601,  1247.75394041,  2194.44578482,  3141.13762922,
        4087.82947363, -4087.82947363, -3141.13762922, -2194.44578482,
       -1247.75394041,  -301.06209601,   645.62974839,  1592.3215928 ,
        2539.0134372 ,  3485.70528161,  4432.39712601, -3743.26182124,
       -2796.56997684, -1849.87813243,  -903.18628803,    43.50555638,
         990.19740078,  1936.88924518,  2883.58108959,  3830.27293399,
       -4345.38601326, -3398.69416885, -2452.00232445, -1505.31048005,
        -558.61863564,   388.07320876,  1334.76505317,  2281.45689757,
       -4000.81836087, -3054.12651647, -2107.43467206, -1160.74282766,
        -214.05098326,   732.64086115, -1762.86701968,  -816.17517527])
y8 = np.array([-5.68015107e+03, -5.68015107e+03, -5.68015107e+03, -5.68015107e+03,
       -4.73345922e+03, -4.73345922e+03, -4.73345922e+03, -4.73345922e+03,
       -4.73345922e+03, -4.73345922e+03, -3.78676738e+03, -3.78676738e+03,
       -3.78676738e+03, -3.78676738e+03, -3.78676738e+03, -3.78676738e+03,
       -3.78676738e+03, -3.78676738e+03, -2.84007553e+03, -2.84007553e+03,
       -2.84007553e+03, -2.84007553e+03, -2.84007553e+03, -2.84007553e+03,
       -2.84007553e+03, -2.84007553e+03, -2.84007553e+03, -1.89338369e+03,
       -1.89338369e+03, -1.89338369e+03, -1.89338369e+03, -1.89338369e+03,
       -1.89338369e+03, -1.89338369e+03, -1.89338369e+03, -1.89338369e+03,
       -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
       -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
       -9.46691844e+02, -9.46691844e+02, -4.54747351e-13, -4.54747351e-13,
       -4.54747351e-13, -4.54747351e-13, -4.54747351e-13, -4.54747351e-13,
       -4.54747351e-13, -4.54747351e-13, -4.54747351e-13,  9.46691844e+02,
        9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
        9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
        9.46691844e+02,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
        1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
        1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  2.84007553e+03,
        2.84007553e+03,  2.84007553e+03,  2.84007553e+03,  2.84007553e+03,
        2.84007553e+03,  2.84007553e+03,  2.84007553e+03,  2.84007553e+03,
        3.78676738e+03,  3.78676738e+03,  3.78676738e+03,  3.78676738e+03,
        3.78676738e+03,  3.78676738e+03,  3.78676738e+03,  3.78676738e+03,
        4.73345922e+03,  4.73345922e+03,  4.73345922e+03,  4.73345922e+03,
        4.73345922e+03,  4.73345922e+03,  5.68015107e+03,  5.68015107e+03])

ax7.plot(bx4,by4,'--k',linewidth=0.5)
ax8.plot(bx6,by6,'--k',linewidth=0.5)
ax9.plot(bx8,by8,'--k',linewidth=0.5)
for i in range(nTurbs):
    circ = plt.Circle((x4[i],y4[i]),rotorDiameter[0]/2.,facecolor=robin,edgecolor=robin,alpha=0.5)
    ax7.add_patch(circ)
    circ = plt.Circle((x6[i],y6[i]),rotorDiameter[0]/2.,facecolor=robin,edgecolor=robin,alpha=0.5)
    ax8.add_patch(circ)
    circ = plt.Circle((x8[i],y8[i]),rotorDiameter[0]/2.,facecolor=robin,edgecolor=robin,alpha=0.5)
    ax9.add_patch(circ)

# plt.suptitle('average turbine spacing (rotor diameters)',fontsize=10,y=0.998)
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

ax4.set_title('4',fontsize=8,rotation=0.)
ax5.set_title('6',fontsize=8,rotation=0.)
ax6.set_title('8',fontsize=8,rotation=0.)

ax8.set_xlabel('average turbine spacing (rotor diameters)',fontsize=8)




min4 = np.min(np.array([np.min(f_param4),np.min(f_direct4),np.min(f_grid4)]))
min6 = np.min(np.array([np.min(f_param6),np.min(f_direct6),np.min(f_grid6)]))
min8 = np.min(np.array([np.min(f_param8),np.min(f_direct8),np.min(f_grid8)]))

max4 = np.max(np.array([np.max(f_param4),np.max(f_direct4),np.max(f_grid4)]))
max6 = np.max(np.array([np.max(f_param6),np.max(f_direct6),np.max(f_grid6)]))
max8 = np.max(np.array([np.max(f_param8),np.max(f_direct8),np.max(f_grid8)]))

# min4 = np.min(np.array([np.min(f_direct4),np.min(f_grid4)]))
# min6 = np.min(np.array([np.min(f_direct6),np.min(f_grid6)]))
# min8 = np.min(np.array([np.min(f_direct8),np.min(f_grid8)]))
#
# max4 = np.max(np.array([np.max(f_direct4),np.max(f_grid4)]))
# max6 = np.max(np.array([np.max(f_direct6),np.max(f_grid6)]))
# max8 = np.max(np.array([np.max(f_direct8),np.max(f_grid8)]))

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

# min4 = np.min(np.array([np.min(AEP_direct4),np.min(AEP_grid4)]))
# min6 = np.min(np.array([np.min(AEP_direct6),np.min(AEP_grid6)]))
# min8 = np.min(np.array([np.min(AEP_direct8),np.min(AEP_grid8)]))
#
# max4 = np.max(np.array([np.max(AEP_direct4),np.max(AEP_grid4)]))
# max6 = np.max(np.array([np.max(AEP_direct6),np.max(AEP_grid6)]))
# max8 = np.max(np.array([np.max(AEP_direct8),np.max(AEP_grid8)]))

b4 = np.linspace(min4,max4,25)
b6 = np.linspace(min6,max6,25)
b8 = np.linspace(min8,max8,25)

b6 = np.linspace(min4,max4,25)
b8 = np.linspace(min4,max4,25)


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


ax4.text(0.907,62.,'a',fontsize=11,color='black')
ax5.text(0.907,62.,'b',fontsize=11,color='black')
ax6.text(0.907,62.,'c',fontsize=11,color='black')
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

plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/results_spacing.pdf',transparent=True)
plt.show()
