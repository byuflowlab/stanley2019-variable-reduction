import numpy as np
import matplotlib.pyplot as plt
import math
import random
import warnings
import grid_param_fortran
warnings.simplefilter('ignore', np.RankWarning)
import sys
import grid_param_fortran
sys.path.insert(0, '/home/flowlab/PJ/reduction')
from position_constraints import *
sys.dont_write_bytecode = True


def startLoop(N,d,my,step,plus_minus,rotate,turbs_per_row,x_start,y0,boundaryVertices,boundaryNormals,t):
    d = d+plus_minus*d*step
    shear = np.tan(np.deg2rad(20.))*d
    X,Y = grid_param_fortran.makegrid_fortran(N, d, d*my, shear, rotate, turbs_per_row, x_start, y0)
    bounds = arbitraryBoundary(X, Y, boundaryVertices, boundaryNormals)
    num_inside = 0

    # plt.figure(1)
    # plt.clf()
    # for i in range(len(X)):
    #     circ = plt.Circle((X[i],Y[i]),130./2.)
    #     plt.gca().add_patch(circ)
    # bx = boundaryVertices[:,0]
    # by = boundaryVertices[:,1]
    # bx = np.append(bx,bx[0])
    # by = np.append(by,by[0])
    # plt.plot(bx,by,'--k')
    # plt.axis('equal')
    # plt.pause(0.0001)

    for i in range(t**2):
        if min(bounds[i]) >= 0.:
            num_inside += 1
    return num_inside,d


#fill an arbirary boundary with a baseline grid
def startGrid(nTurbines,boundaryVertices,boundaryNormals,rotate=0.,my=1.):
    xBounds = boundaryVertices[:, 0]
    yBounds = boundaryVertices[:, 1]

    t = 25
    range_x = (max(xBounds)-min(xBounds))
    range_y = (max(yBounds)-min(yBounds))
    d = 1.5*max(np.array([range_x,range_y]))/float(t)
    # rotate = 0.
    nRows = t
    turbs_per_row = np.ones(nRows,dtype=int)*t
    N = np.sum(turbs_per_row)

    # x_start = np.ones(t)*min(xBounds)-range_x*0.5
    x_start = np.zeros(t)
    y0 = min(yBounds)-0.05*range_y

    shear = np.tan(np.deg2rad(20.))*d
    X,Y = grid_param_fortran.makegrid_fortran(N, d, d*my, shear, rotate, turbs_per_row, x_start, y0)

    bounds = arbitraryBoundary(X, Y, boundaryVertices, boundaryNormals)
    num_inside = 0
    for i in range(t**2):
        if min(bounds[i]) >= 0.:
            num_inside += 1

    converged = False
    step = 0.25
    print 'setting up base grid'
    # print nTurbines
    while converged == False:
        while num_inside < nTurbines:
            num_inside,d = startLoop(N,d,my,step,-1.,rotate,turbs_per_row,x_start,y0,boundaryVertices,boundaryNormals,t)
        while num_inside >= nTurbines:
            num_inside,d = startLoop(N,d,my,step,1.,rotate,turbs_per_row,x_start,y0,boundaryVertices,boundaryNormals,t)

        if step < 1.E-8:
            converged = True
        else:
            step = step/2.
        # print step

    print 'final convergence'
    num_inside,d = startLoop(N,d,my,step,-1.,rotate,turbs_per_row,x_start,y0,boundaryVertices,boundaryNormals,t)
    small_step = step/100.
    while num_inside >= nTurbines:
        num_inside,d = startLoop(N,d,my,small_step,1.,rotate,turbs_per_row,x_start,y0,boundaryVertices,boundaryNormals,t)

    while num_inside < nTurbines:
        d = d-d*small_step
        shear = np.tan(np.deg2rad(20.))*d
        X,Y = grid_param_fortran.makegrid_fortran(N, d, d*my, shear, rotate, turbs_per_row, x_start, y0)
        bounds = arbitraryBoundary(X, Y, boundaryVertices, boundaryNormals)
        num_inside = 0
        for i in range(t**2):
            if min(bounds[i]) > 0.:
                num_inside += 1

    x = np.zeros(nTurbines)
    y = np.zeros(nTurbines)
    turb = 0
    index = 0
    while turb < nTurbines:
        if min(bounds[index]) >= 0.:
                x[turb] = X[index]
                y[turb] = Y[index]
                turb += 1
        index += 1

    dx = d
    dy = d*my
    turbineX = np.cos(-np.deg2rad(rotate))*x - np.sin(-np.deg2rad(rotate))*y
    turbineY = np.sin(-np.deg2rad(rotate))*x + np.cos(-np.deg2rad(rotate))*y
    y0 = np.min(turbineY)
    turbs_per_row = np.array([1])
    for i in range(nTurbines-1):
        if turbineY[i+1] < turbineY[i]+1E-4 and turbineY[i+1] > turbineY[i]-1E-4:
            turbs_per_row[-1] += 1
        else:
            turbs_per_row = np.append(turbs_per_row,1)

    shear = dy*np.tan(np.deg2rad(20.))
    nRows = len(turbs_per_row)
    x_start = np.zeros(nRows)
    x_start[0] = 0.0
    index = 1
    # for i in range(nTurbines-1):
    #     if not (turbineY[i+1] < turbineY[i]+1.E-4 and turbineY[i+1] > turbineY[i]-1.E-4):
    #         xs = turbineX[i+1]
    #         x_start[index] = (xs-shear-turbineX[sum(turbs_per_row[0:index-1])])/dx
    #         index += 1
    for i in range(nRows-1):
        x_start[i+1] = x_start[i] + (turbineX[sum(turbs_per_row[0:i+1])] - shear - turbineX[sum(turbs_per_row[0:i])])/dx

    # plt.figure(1)
    # for i in range(nTurbines):
    #     circ = plt.Circle((turbineX[i],turbineY[i]),130./2.,alpha=0.5,color='blue')
    #     plt.gca().add_patch(circ)
    #     circ = plt.Circle((x[i],y[i]),130./2.,alpha=0.5,color='red')
    #     plt.gca().add_patch(circ)
    # bx = boundaryVertices[:,0]
    # by = boundaryVertices[:,1]
    # bx = np.append(bx,bx[0])
    # by = np.append(by,by[0])
    # plt.plot(bx,by,'--k')
    # plt.axis('equal')
    # plt.show()

    return dx,dy,shear,rotate,turbs_per_row,x_start,y0


#makes a grid of points/turbines
def makeGrid(dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start):
    nRows = len(turbs_per_row)
    nTurbs = sum(turbs_per_row)
    x = np.zeros(nTurbs)
    y = np.zeros(nTurbs)
    index = 0
    for i in range(nRows):
        for j in range(turbs_per_row[i]):
            x[index] = x_start[i] + dx*j + float(i)*dy*np.tan(np.deg2rad(shearX))
            y[index] = y0 + dy*i + (float(j)+col_start[i])*dx*np.tan(np.deg2rad(shearY))
            index += 1
    xc = x - (max(x)+min(x))/2.
    yc = y - (max(y)+min(y))/2.
    xr = np.cos(np.deg2rad(rotate))*xc - np.sin(np.deg2rad(rotate))*yc
    yr = np.sin(np.deg2rad(rotate))*xc + np.cos(np.deg2rad(rotate))*yc
    return xr,yr


def find_lengths(x,y,nPoints):
    length = np.zeros(len(x)-1)
    for i in range(nPoints):
        length[i] = np.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)
    return length


#places turbines spaced equally around the boundray of a wind farm
def makeBoundary(xBounds,yBounds,start,diameter,nTurbs,frac=0.5):

    if xBounds[-1] != xBounds[0]:
        xBounds = np.append(xBounds,xBounds[0])
        yBounds = np.append(yBounds,yBounds[0])

    nBounds = len(xBounds)
    lenBound = find_lengths(xBounds,yBounds,len(xBounds)-1)
    circumference = sum(lenBound)
    spacing = 2.*np.sqrt(2.)
    x = np.array([])
    y = np.array([])

    N = int(circumference/(spacing*diameter))
    if N > nTurbs*frac:
        N = int(nTurbs*frac)
    bound_loc = np.linspace(start,start+circumference-circumference/float(N),N)
    for i in range(N):
        if bound_loc[i] > circumference:
            bound_loc[i] = bound_loc[i]%circumference
        while bound_loc[i] < 0.:
            bound_loc[i] += circumference

    for i in range(N):
        done = False
        for j in range(nBounds):
            if done == False:
                if bound_loc[i] < sum(lenBound[0:j+1]):
                    point_x = xBounds[j] + (xBounds[j+1]-xBounds[j])*(bound_loc[i]-sum(lenBound[0:j]))/lenBound[j]
                    point_y = yBounds[j] + (yBounds[j+1]-yBounds[j])*(bound_loc[i]-sum(lenBound[0:j]))/lenBound[j]
                    done = True
                    x = np.append(x,point_x)
                    y = np.append(y,point_y)

    return x,y


#determines the starting grid variables
def startTurbines(xBounds,yBounds,start,diameter,nTurbs,frac=0.5,my=2.):

    boundX, boundY = makeBoundary(xBounds,yBounds,start,diameter,nTurbs,frac=frac)
    turbsRemaining = nTurbs - len(boundX)

    nBounds = len(xBounds)
    locations = np.zeros((nBounds,2))
    locations[:, 0] = xBounds*0.8
    locations[:, 1] = yBounds*0.8
    boundaryVertices, boundaryNormals = calculate_boundary(locations)
    dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start = startGrid(turbsRemaining,boundaryVertices,boundaryNormals,my=my)

    return dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start


#defines the layout of every turbine in a wind farm with the boundary section and the grid section
def turbineLocs(xBounds,yBounds,start,diameter,nTurbs,dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start,frac=0.5):
    boundX, boundY = makeBoundary(xBounds,yBounds,start,diameter,nTurbs,frac=frac)
    gridX,gridY = makeGrid(dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start)
    turbineX = np.append(boundX,gridX)
    turbineY = np.append(boundY,gridY)

    return turbineX,turbineY



if __name__=="__main__":
    # locations = np.loadtxt('layout_amalia.txt')
    # xBounds = locations[:, 0]
    # yBounds = locations[:, 1]
    # xBounds = xBounds - min(xBounds) - (max(xBounds)-min(xBounds))/2.
    # yBounds = yBounds - min(yBounds) - (max(yBounds)-min(yBounds))/2.

    circle_radius = 5000.
    nBounds = 20
    locations = np.zeros((nBounds,2))
    xBounds = np.zeros(nBounds)
    yBounds = np.zeros(nBounds)
    theta = np.linspace(0.,2.*np.pi-2.*np.pi/float(nBounds),nBounds)
    for i in range(nBounds):
        xBounds[i] = circle_radius*np.cos(theta[i])
        yBounds[i] = circle_radius*np.sin(theta[i])

    locations[:, 0] = xBounds
    locations[:, 1] = yBounds

    boundaryVertices, boundaryNormals = calculate_boundary(locations)
    xBounds = boundaryVertices[:, 0]
    yBounds = boundaryVertices[:, 1]
    start = 0.
    diameter = 130.
    nTurbs = 100

    # dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start = startGrid(50,boundaryVertices,boundaryNormals)
    # print 'shearX: ', shearX
    # print 'shearY: ', shearY
    # print 'turbs_per_row: ', turbs_per_row
    # print 'dy: ', dy
    # # x_start = np.zeros(len(turbs_per_row))
    # turbineX, turbineY = makeGrid(dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start)
    # nTurbs = len(turbineX)


    dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start = startTurbines(xBounds,yBounds,start,diameter,nTurbs,frac=0.5,my=1.)
    turbineX, turbineY = turbineLocs(xBounds,yBounds,start,diameter,nTurbs,dx,dy,shearX,shearY,rotate,turbs_per_row,x_start,y0,col_start,frac=0.5)

    xBounds = boundaryVertices[:,0]
    yBounds = boundaryVertices[:,1]
    xBounds = np.append(xBounds,xBounds[0])
    yBounds = np.append(yBounds,yBounds[0])
    plt.plot(xBounds,yBounds,'--k')
    plt.plot(turbineX,turbineY,'o')
    for i in range(nTurbs):
        plt.text(turbineX[i],turbineY[i],'%s'%i)
    plt.axis('equal')
    plt.show()
