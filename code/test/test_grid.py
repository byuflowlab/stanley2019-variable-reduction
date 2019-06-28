import numpy as np
import constraints
import sys
sys.dont_write_bytecode = True
import time
import sys
import grid_param_fortran
sys.path.insert(1, '/home/flowlab/PJ/reduction')
from position_constraints import *
sys.path.insert(0, '/home/flowlab/PJ/reduction/exact-gradients')
from var_reduction_exact import *
import matplotlib.pyplot as plt


nTurbines = 100
locations = np.loadtxt('/home/flowlab/PJ/reduction/layout_amalia.txt')
boundaryVertices, boundaryNormals = calculate_boundary(locations)

bx = boundaryVertices[:,0]
by = boundaryVertices[:,1]
bx = bx-(min(bx)+max(bx))/2.
by = by-(min(by)+max(by))/2.
boundaryVertices[:,0] = bx[:]
boundaryVertices[:,1] = by[:]

dx,dy,shear,rotate,turbs_per_row,x_start,y0 = startGrid(nTurbines,
                    boundaryVertices,boundaryNormals,my=1.)

turbineX, turbineY = grid_param_fortran.makegrid_fortran(nTurbines, dx, dy, shear, rotate, turbs_per_row, x_start, y0)

plt.figure(1)
for i in range(nTurbines):
    circ = plt.Circle((turbineX[i],turbineY[i]),50.)
    plt.gca().add_patch(circ)
bx = np.append(bx,bx[0])
by = np.append(by,by[0])
plt.plot(bx,by,'--k')
plt.axis('equal')
plt.show()

# turbineZ = np.ones(nTurbines)*100.
# rotorDiameter = np.ones(nTurbines)*130.
# windDirections = np.array([0.])
# windSpeeds = np.array([10.])
# windFrequencies = np.array([1.])
# shearExp = 0.15
# relaxationFactor = 1.0
# rated_ws = 10.
# rated_power = 3.35
# cut_in_speed = 3.
# cut_out_speed = 25.
# zref = 100.
# z0 = 0.
# AEP = gaus_aep_grid.calcaep_grid(turbineZ, rotorDiameter, windDirections, dx, dy, shear, rotate,
#                     turbs_per_row, x_start, y0, windSpeeds, windFrequencies, shearExp, relaxationFactor,
#                     rated_ws, rated_power, cut_in_speed, cut_out_speed, zref, z0)
#
# print AEP
# #
# # subroutine calcAEP_grid(nTurbines, nDirections, nRows, turbineZ, rotorDiameter, windDirections,&
# #             &dx, dy, shear, rotate, turbs_per_row, x_start, y0,&
# #             &windSpeeds, windFrequencies, shearExp, relaxationFactor, rated_ws, rated_power,&
# #             &cut_in_speed, cut_out_speed, zref, z0, AEP)
