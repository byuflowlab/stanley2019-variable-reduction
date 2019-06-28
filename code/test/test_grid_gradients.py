import numpy as np
import constraints
import sys
sys.dont_write_bytecode = True
import time
import sys
sys.path.insert(1, '/home/flowlab/PJ/reduction')
from position_constraints import *
sys.path.insert(0, '/home/flowlab/PJ/reduction/exact-gradients')
from var_reduction_exact import *

import grid_param_fortran

nTurbines = 15
locations = np.loadtxt('/home/flowlab/PJ/reduction/layout_amalia.txt')
boundaryVertices, boundaryNormals = calculate_boundary(locations)

bx = boundaryVertices[:,0]
by = boundaryVertices[:,1]
bx = bx-(min(bx)+max(bx))/2.
by = by-(min(by)+max(by))/2.
boundaryVertices[:,0] = bx[:]
boundaryVertices[:,1] = by[:]

dx,dy,shear,rotate,turbs_per_row,x_start,y0 = startGrid(nTurbines,
                    boundaryVertices,boundaryNormals,my=1.0)
rotate = float(np.random.rand(1))*10.-5.

dxd = 1.
dyd = 0.
sheard = 0.
rotated = 0.
turbinex,dx_ddx,turbiney,dy_ddx = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

dxd = 0.
dyd = 1.
sheard = 0.
rotated = 0.
turbinex,dx_ddy,turbiney,dy_ddy = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

dxd = 0.
dyd = 0.
sheard = 1.
rotated = 0.
turbinex,dx_dshear,turbiney,dy_dshear = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

dxd = 0.
dyd = 0.
sheard = 0.
rotated = 1.
turbinex,dx_drotate,turbiney,dy_drotate = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

dxd = 0.
dyd = 0.
sheard = 0.
rotated = 0.
d = 1.E-5
delx,_,dely,_ = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx+d,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)
fdx_dx = -(turbinex-delx)/d
fdy_dx = -(turbiney-dely)/d

delx,_,dely,_ = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy+d,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)
fdx_dy = -(turbinex-delx)/d
fdy_dy = -(turbiney-dely)/d

delx,_,dely,_ = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear+d,sheard,rotate,rotated,turbs_per_row,x_start,y0)
fdx_dshear = -(turbinex-delx)/d
fdy_dshear = -(turbiney-dely)/d

delx,_,dely,_ = grid_param_fortran.makegrid_fortran_dv(nTurbines,dx,dxd,dy,dyd,shear,sheard,rotate+d,rotated,turbs_per_row,x_start,y0)
fdx_drotate = -(turbinex-delx)/d
fdy_drotate = -(turbiney-dely)/d


# print 'turbineX'
# print 'fdx_dx: ', fdx_dx
# print 'dx_ddx: ', dx_ddx
print fdx_dx-dx_ddx

# print 'fdx_dy: ', fdx_dy
# print 'dx_ddy: ', dx_ddy
print fdx_dy-dx_ddy

# print 'fdx_dshear: ', fdx_dshear
# print 'dx_ddshear: ', dx_dshear
print fdx_dshear-dx_dshear

# print 'fdx_drotate: ', fdx_drotate
# print 'dx_drotate: ', dx_drotate
print fdx_drotate-dx_drotate


# print 'turbineY'
# print 'fdy_dx: ', fdy_dx
# print 'dy_ddx: ', dy_ddx
print fdy_dx-dy_ddx

# print 'fdy_dy: ', fdy_dy
# print 'dy_ddy: ', dy_ddy
print fdy_dy-dy_ddy

# print 'fdy_dshear: ', fdy_dshear
# print 'dy_dshear: ', dy_dshear
print fdy_dshear-dy_dshear

# print 'fdy_drotate: ', fdy_drotate
# print 'dy_drotate: ', dy_drotate
print fdy_drotate-dy_drotate
