import numpy as np
import constraints
import sys
sys.dont_write_bytecode = True
import time
import sys
sys.path.insert(0, '/fslhome/pjstanle/compute/reduction')
from position_constraints import *

nTurbines = 5
turbineX = np.random.rand(nTurbines)*1200.
turbineY = np.random.rand(nTurbines)*1200.

#circular boundary
boundaryVertices = np.zeros((1,2))
boundaryVertices[0] = np.array([1000.,0.])
boundaryNormals = np.zeros((1,2))
boundaryNormals[0] = np.array([0.,0.])

s_fortran, b_fortran = constraints.constraints_position(turbineX, turbineY, boundaryVertices, boundaryNormals)

# FD
s_fd_x = np.zeros((nTurbines,((nTurbines-1)*nTurbines/2)))
s_fd_y = np.zeros((nTurbines,((nTurbines-1)*nTurbines/2)))
b_fd_x = np.zeros((nTurbines,nTurbines,1))
b_fd_y = np.zeros((nTurbines,nTurbines,1))
d = 1E-5
for i in range(nTurbines):
    x = np.zeros(nTurbines)
    y = np.zeros(nTurbines)
    x[:] = turbineX
    y[:] = turbineY
    x[i] += d
    y[i] += d
    s_x, b_x = constraints.constraints_position(x, turbineY, boundaryVertices, boundaryNormals)
    s_y, b_y = constraints.constraints_position(turbineX, y, boundaryVertices, boundaryNormals)
    s_fd_x[i] = (s_x-s_fortran)/d
    s_fd_y[i] = (s_y-s_fortran)/d
    b_fd_x[i] = (b_x-b_fortran)/d
    b_fd_y[i] = (b_y-b_fortran)/d


dx = np.zeros((nTurbines,nTurbines))
dy = np.eye(nTurbines)
ss,ss_d,bd,bd_d = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                        boundaryVertices,boundaryNormals)

print ss_d-s_fd_y
print bd_d-b_fd_y

dx = np.eye(nTurbines)
dy = np.zeros((nTurbines,nTurbines))
ss,ss_d,bd,bd_d = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                        boundaryVertices,boundaryNormals)

print ss_d-s_fd_x
print bd_d-b_fd_x
