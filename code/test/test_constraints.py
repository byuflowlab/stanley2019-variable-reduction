import numpy as np
import constraints
import sys
sys.dont_write_bytecode = True
import time
import sys
# sys.path.insert(0, '/home/flowlab/PJ/reduction')
sys.path.insert(0, '/home/flowlab/PJ/reduction/exact-gradients')
from position_constraints import *

nTurbines = 15
turbineX = np.random.rand(nTurbines)*1200.
turbineY = np.random.rand(nTurbines)*1200.


#spacing
ss_fortran = constraints.turbine_spacing_squared(turbineX,turbineY)
ss_python = turbineSpacingSquared(turbineX,turbineY)
print ss_fortran - ss_python

#boundary circle
boundaryVertices = np.zeros((1,2))
boundaryVertices[0] = np.array([1000.,0.])
boundaryNormals = np.zeros((1,2))
boundaryNormals[0] = np.array([0.,0.])
bc_fortran = constraints.boundary_distances(turbineX, turbineY, boundaryVertices, boundaryNormals)
bc_python = circularBoundary(turbineX, turbineY, 1000.)
print np.ndarray.flatten(bc_fortran) - bc_python

#boundary polygon
locations = np.loadtxt('/home/flowlab/PJ/reduction/layout_amalia.txt')
boundaryVertices, boundaryNormals = calculate_boundary(locations)
bp_fortran = constraints.boundary_distances(turbineX, turbineY, boundaryVertices, boundaryNormals)
bp_python = arbitraryBoundary(turbineX, turbineY, boundaryVertices, boundaryNormals)
print bp_fortran - bp_python

s_fortran, b_fortran = constraints.constraints_position(turbineX, turbineY, boundaryVertices, boundaryNormals)
print ss_python-s_fortran
print bp_python-b_fortran
