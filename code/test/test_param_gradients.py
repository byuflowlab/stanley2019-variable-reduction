import numpy as np
import constraints
import sys
sys.dont_write_bytecode = True
import time
import sys
import param_fortran
sys.path.insert(1, '/home/flowlab/PJ/reduction')
from position_constraints import *
sys.path.insert(0, '/home/flowlab/PJ/reduction/exact-gradients')
from var_reduction_exact import *
import matplotlib.pyplot as plt


locations = np.loadtxt('/home/flowlab/PJ/reduction/layout_amalia.txt')
boundaryVertices, boundaryNormals = calculate_boundary(locations)

bx = boundaryVertices[:,0]
by = boundaryVertices[:,1]
bx = bx-(min(bx)+max(bx))/2.
by = by-(min(by)+max(by))/2.
bx = np.append(bx,bx[0])
by = np.append(by,by[0])

start = 0.
nTurbines = 20
n = nTurbines/2

dx,dy,shear,rotate,turbs_per_row,x_start,y0 = startGrid(n,
                    boundaryVertices*0.95,boundaryNormals,my=1.)

dxd = 1.
dyd = 0.
sheard = 0.
rotated = 0.
startd = 0.
turbinex,dx_ddx,turbiney,dy_ddx = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)

dxd = 0.
dyd = 1.
sheard = 0.
rotated = 0.
startd = 0.
turbinex,dx_ddy,turbiney,dy_ddy = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)

dxd = 0.
dyd = 0.
sheard = 1.
rotated = 0.
startd = 0.
turbinex,dx_dshear,turbiney,dy_dshear = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)

dxd = 0.
dyd = 0.
sheard = 0.
rotated = 1.
startd = 0.
turbinex,dx_drotate,turbiney,dy_drotate = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)

dxd = 0.
dyd = 0.
sheard = 0.
rotated = 0.
startd = 1.
turbinex,dx_dstart,turbiney,dy_dstart = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)

startd = 0.

d = 1.E-5
DX, _, DY, _ = param_fortran.turbinelocations_dv(nTurbines,n,n,dx+d,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)
fdx_x = (DX-turbinex)/d
fdy_x = (DY-turbiney)/d

DX, _, DY, _ = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy+d,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)
fdx_y = (DX-turbinex)/d
fdy_y = (DY-turbiney)/d

DX, _, DY, _ = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear+d,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)
fdx_shear = (DX-turbinex)/d
fdy_shear = (DY-turbiney)/d

DX, _, DY, _ = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate+d,rotated,turbs_per_row,x_start,y0,
                            start,startd,bx,by)
fdx_rotate = (DX-turbinex)/d
fdy_rotate = (DY-turbiney)/d

DX, _, DY, _ = param_fortran.turbinelocations_dv(nTurbines,n,n,dx,dxd,
                            dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                            start+d,startd,bx,by)
fdx_start = (DX-turbinex)/d
fdy_start = (DY-turbiney)/d


print "turbineX"
print fdx_x-dx_ddx
print fdx_y-dx_ddy
print fdx_shear-dx_dshear
print fdx_rotate-dx_drotate
print fdx_start-dx_dstart

print "turbineY"
print fdy_x-dy_ddx
print fdy_y-dy_ddy
print fdy_shear-dy_dshear
print fdy_rotate-dy_drotate
print fdy_start-dy_dstart
