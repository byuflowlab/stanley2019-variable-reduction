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




# N = 50
start = 0.
# turbineX, turbineY = param_fortran.makeboundary(N,bx,by,start)
nTurbines = 200
n = nTurbines/2
dx,dy,shear,rotate,turbs_per_row,x_start,y0 = startGrid(n,
                    boundaryVertices*0.95,boundaryNormals,my=1.)

turbineX,turbineY = param_fortran.turbinelocations(nTurbines,n,n,dx,dy,shear,rotate,turbs_per_row,
                                x_start,y0,start,bx,by)


plt.figure(1)
for i in range(len(turbineX)):
    circ = plt.Circle((turbineX[i],turbineY[i]),50.)
    plt.gca().add_patch(circ)
plt.plot(bx,by,'--k')
# plt.plot(bx[0],by[0],'or')
# plt.plot(bx[1],by[1],'ob')
plt.axis('equal')
plt.show()
