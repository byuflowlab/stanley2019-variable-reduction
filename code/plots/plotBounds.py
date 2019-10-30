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
tuscan = 'C2'
royal = 'C3'

# robin = '#ACBD78'

nTurbs = 100
rotorDiameter = np.ones(nTurbs)*130.

spacing = 4.
side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
a = side_length**2
circle_radius = np.sqrt(a/np.pi)

# fig = plt.figure(figsize=[6.5,3.])
fig = plt.figure(figsize=[4.5,2.])
ax1 = plt.subplot2grid((1, 3), (0, 0))
ax2 = plt.subplot2grid((1, 3), (0, 1),sharex=ax1,sharey=ax1)
ax3 = plt.subplot2grid((1, 3), (0, 2),sharex=ax1,sharey=ax1)

ax1.axis('equal')
ax2.axis('equal')
ax3.axis('equal')

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

ax1.plot(xBounds,yBounds,'-',alpha=0.5,color=robin)


"""circle boundary"""
nBounds = 200
xBounds = np.zeros(nBounds)
yBounds = np.zeros(nBounds)
theta = np.linspace(0.,2.*np.pi,nBounds)
for i in range(nBounds):
    xBounds[i] = circle_radius*np.cos(theta[i])
    yBounds[i] = circle_radius*np.sin(theta[i])

ax2.plot(xBounds,yBounds,'-',alpha=0.5,color=robin)


"""square boundary rotated 30 deg from dominant wind direction"""
nBounds = 4
x = np.array([-side_length/2.,side_length/2.,side_length/2.,-side_length/2.])
y = np.array([-side_length/2.,-side_length/2.,side_length/2.,side_length/2.])
xBounds = x*np.cos(np.deg2rad(30.)) - y*np.sin(np.deg2rad(30.))
yBounds = x*np.sin(np.deg2rad(30.)) + y*np.cos(np.deg2rad(30.))

xBounds = np.append(xBounds,xBounds[0])
yBounds = np.append(yBounds,yBounds[0])

ax3.plot(xBounds,yBounds,'-',alpha=0.5,color=robin)

ax1.set_ylim(-5000.,5000.)
# ax2.set_ylim(-3300.,3300.)
# ax3.set_ylim(-3300.,3300.)

ax1.set_xlim(-5000.,5000.)
# ax2.set_xlim(-3300.,3300.)
# ax3.set_xlim(-3300.,3300.)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)

ax1.set_xticks(())
ax2.set_xticks(())
ax3.set_xticks(())
ax1.set_yticks(())
ax2.set_yticks(())
ax3.set_yticks(())


ax1.set_xlabel('Amalia',fontsize=8,labelpad=-10.)
ax2.set_xlabel('circle',fontsize=8,labelpad=-10.)
ax3.set_xlabel('square',fontsize=8,labelpad=-10.)

plt.subplots_adjust(top = 1.0, bottom = 0.1, right = 0.96, left = 0.04,
            hspace = 0.0, wspace = 0.)

ax1.text(-2800.,2800.,'A',fontsize=14,weight='bold',color='black')
ax2.text(-2800.,2800.,'B',fontsize=14,weight='bold',color='black')
ax3.text(-2800.,2800.,'C',fontsize=14,weight='bold',color='black')

# plt.savefig('/Users/ningrsrch/Dropbox/Projects/reduction/paper/paper-figures/bounds.pdf',transparent=True)
plt.show()
