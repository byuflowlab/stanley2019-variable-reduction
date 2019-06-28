from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2, CONMIN
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from aep_calc import *
import constraints
import sys
sys.path.insert(0, '/home/flowlab/PJ/reduction/')
from windRoses import *
from position_constraints import calculate_boundary
sys.dont_write_bytecode = True


def random_start(N,D,xmin,xmax,ymin,ymax):
    x = np.zeros(N)
    y = np.zeros(N)
    i = 0
    while i < N:
        good = True
        xtemp = float(np.random.rand(1))*(xmax-xmin)+xmin
        ytemp = float(np.random.rand(1))*(ymax-ymin)+ymin
        for j in range(i):
            dist = np.sqrt((x[j]-xtemp)**2+(y[j]-ytemp)**2)
            if dist < D:
                good = False
        if good == True:
            x[i] = xtemp
            y[i] = ytemp
            i += 1
    return x,y


def objfunc(x):
    global rotorDiameter
    global turbineZ
    global windDirections
    global windSpeeds
    global windFrequencies
    global boundaryVertices
    global boundaryNormals
    global nCalls
    global index

    nCalls += 1

    nTurbs = len(x['x'])/2
    turbineX = x['x'][0:nTurbs]
    turbineY = x['x'][nTurbs:2*nTurbs]

    # show = True
    # if show == True:
    #     boundx = boundaryVertices[:,0]
    #     boundy = boundaryVertices[:,1]
    #     boundx = np.append(boundx,boundx[0])
    #     boundy = np.append(boundy,boundy[0])
    #
    #     plt.figure(1)
    #     plt.clf()
    #     for i in range(nTurbs):
    #         circ = plt.Circle((turbineX[i],turbineY[i]), rotorDiameter[i]/2.,facecolor="blue",edgecolor="blue",alpha=0.2)
    #         plt.gca().add_patch(circ)
    #     plt.plot(boundx,boundy,'--k')
    #     plt.axis('equal')
    #     plt.draw()
    #     plt.pause(0.001)

    AEP = fast_calc_AEP(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies)
    obj = -AEP
    print obj
    ss, bd = constraints.constraints_position(turbineX, turbineY, boundaryVertices, boundaryNormals)
    print '1'
    bounds = np.zeros(nTurbs)
    index = np.zeros(nTurbs)
    for i in range(nTurbs):
        bounds[i] = np.min(bd[i])
        index[i] = np.argmin(bd[i])
    print '2'
    con = np.append((ss-(2.*rotorDiameter[0])**2)/1.E5,bounds/1.E3)
    fail = 0
    print '3'
    return obj, con


def gradfunc(x,obj,con):
    global rotorDiameter
    global turbineZ
    global windDirections
    global windSpeeds
    global windFrequencies
    global boundaryVertices
    global boundaryNormals
    global nCalls
    global index

    nCalls += 1
    print '4'
    nTurbs = len(x['x'])/2
    turbineX = x['x'][0:nTurbs]
    turbineY = x['x'][nTurbs:2*nTurbs]
    print '5'
    """objective"""
    g_obj = np.zeros(2*nTurbs)

    g_obj[0:nTurbs],g_obj[nTurbs:2*nTurbs] = fast_calc_gradients(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies)
    g_obj = -g_obj
    print '6'
    """constraints"""
    dx = np.eye(nTurbs)
    dy = np.zeros((nTurbs,nTurbs))
    _,ss_dx,_,bd_dx = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                            boundaryVertices,boundaryNormals)

    dx = np.zeros((nTurbs,nTurbs))
    dy = np.eye(nTurbs)
    _,ss_dy,_,bd_dy = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                            boundaryVertices,boundaryNormals)
    print '7'

    g_con = np.zeros(((nTurbs-1)*nTurbs/2 + nTurbs,2*nTurbs))
    for i in range((nTurbs-1)*nTurbs/2):
        g_con[i][0:nTurbs] = ss_dx[:,i]
        g_con[i][nTurbs:2*nTurbs] = ss_dy[:,i]
    print '8'
    for i in range(nTurbs):
        for j in range(nTurbs):
            g_con[i+(nTurbs-1)*nTurbs/2][j] = bd_dx[i][i][index[i]]
            g_con[i+(nTurbs-1)*nTurbs/2][j+nTurbs] = bd_dy[i][i][index[i]]
    print '9'
    fail = False

    return g_obj, g_con, fail


if __name__ == "__main__":
    global rotorDiameter
    global turbineZ
    global windDirections
    global windSpeeds
    global windFrequencies
    global boundaryVertices
    global boundaryNormals
    global nCalls

    run = 1
    nTurbs = 3
    rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    spacing = 4.
    boundary = 'amalia'
    # boundary = 'circle'
    # boundary = 'square'

    folder = 'opt-results/exact/direct/%s_%s_%s'%(boundary,spacing,rose)
    if not os.path.exists(folder):
        os.makedirs(folder)

    # nDirections = 30
    # nSpeeds = 8
    nDirections = 23
    nSpeeds = 5
    if rose == 'northIslandRose':
        windDirections, windFrequencies, windSpeeds = northIslandRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = northIslandRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle
    elif rose == 'ukiahRose':
        windDirections, windFrequencies, windSpeeds = ukiahRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = ukiahRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle
    elif rose == 'victorvilleRose':
        windDirections, windFrequencies, windSpeeds = victorvilleRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = victorvilleRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle

    rotorDiameter = np.ones(nTurbs)*130.
    turbineZ = np.ones(nTurbs)*110.

    minSpacing = 2.0

    side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
    a = side_length**2
    circle_radius = np.sqrt(a/np.pi)
    if boundary == 'circle':
        boundaryVertices = np.zeros((1,2))
        boundaryVertices[0] = np.array([circle_radius,0.])
        boundaryNormals = np.zeros((1,2))
        boundaryNormals[0] = np.array([0.,0.])
    elif boundary == 'square':
        nBounds = 4
        x = np.array([-side_length/2.,side_length/2.,side_length/2.,-side_length/2.])
        y = np.array([-side_length/2.,-side_length/2.,side_length/2.,side_length/2.])
        xBounds = x*np.cos(np.deg2rad(30.)) - y*np.sin(np.deg2rad(30.))
        yBounds = x*np.sin(np.deg2rad(30.)) + y*np.cos(np.deg2rad(30.))

        locations = np.zeros((nBounds,2))
        locations[:, 0] = xBounds
        locations[:, 1] = yBounds
        boundaryVertices, boundaryNormals = calculate_boundary(locations)
    elif boundary == 'amalia':
        locations = np.loadtxt('/home/flowlab/PJ/reduction/layout_amalia.txt')
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
        area_ratio = area/a
        xBounds = xBounds/np.sqrt(area_ratio)
        yBounds = yBounds/np.sqrt(area_ratio)

        locations = np.zeros((len(xBounds),2))
        locations[:, 0] = xBounds
        locations[:, 1] = yBounds
        boundaryVertices, boundaryNormals = calculate_boundary(locations)

    if boundary == 'circle':
        xmin = -circle_radius
        xmax = circle_radius
        ymin = -circle_radius
        ymax = circle_radius
    elif boundary == 'square':
        xmax = side_length/2.
        xmin = -side_length/2.
        ymax = side_length/2.
        ymin = -side_length/2.
    elif boundary =='amalia':
        xmax = max(xBounds)
        xmin = min(xBounds)
        ymax = max(yBounds)
        ymin = min(yBounds)

    x,y = random_start(nTurbs,rotorDiameter[0],xmin,xmax,ymin,ymax)
    nCalls = 0
    turbs = np.append(x,y)
    # obj,con = objfunc(turbs)
    # gradfunc(turbs,obj,con)

    opt_prob = Optimization('AEP max',objfunc)
    opt_prob.addVarGroup('x',nTurbs*2,'c',value=turbs)
    opt_prob.addObj('obj')
    opt_prob.addConGroup('con',(nTurbs-1)*nTurbs/2 + nTurbs,lower=0.)
    opt = SNOPT(sens_type=gradfunc)
    opt(opt_prob)

    # num = 1
    # for i in range(num):
    #     nCalls = 0
    #     print 'iteration: ', i
    #     if boundary == 'circle':
    #         x = np.random.rand(nTurbs)*(2.*circle_radius)-circle_radius
    #         y = np.random.rand(nTurbs)*(2.*circle_radius)-circle_radius
    #     elif boundary == 'square':
    #         x = np.random.rand(nTurbs)*(side_length)-side_length/2.
    #         y = np.random.rand(nTurbs)*(side_length)-side_length/2.
    #     elif boundary == 'amalia':
    #         x = np.random.rand(nTurbs)*(max(xBounds)-min(xBounds))+min(xBounds)
    #         y = np.random.rand(nTurbs)*(max(yBounds)-min(yBounds))+min(yBounds)
    #
    #     """Optimization"""
    #     optProb = Optimization('Wind_Farm_AEP', obj_func)
    #     optProb.addObj('obj')
    #
    #     optProb.addVarGroup('x', nTurbs, type='c', lower=xmin, upper=xmax, value=x)
    #     optProb.addVarGroup('y', nTurbs, type='c', lower=ymin, upper=ymax, value=y)
    #
    #     num_cons_sep = (nTurbs-1)*nTurbs/2
    #     optProb.addConGroup('sep', num_cons_sep, lower=(2.*rotorDiameter[0])**2, upper=None)
    #     optProb.addConGroup('bound', nTurbs*len(boundaryVertices), lower=0., upper=None)
    #
    #     opt = SNOPT()
    #     opt.setOption('Scale option',0)
    #     opt.setOption('Derivative level',3)
    #     opt.setOption('Iterations limit',1000000)
    #
    #     # opt.setOption('Summary file','current_summary_square_gaus.out')
    #     opt.setOption('Major optimality tolerance',1.e-5)
    #     opt.setOption('Major feasibility tolerance',1.e-6)
    #
    #     # res = opt(optProb)
    #     opt(optProb,sens_type=grad_func)
    #
    #     x = res.xStar['x']
    #     y = res.xStar['y']
    #
    #     input = {'x':x,'y':y}
    #     funcs,_ = obj_func(input)
    #
    #     separation = min(funcs['sep'])
    #     boundary = min(funcs['bound'])
    #
    #     print 'AEP opt: ', -funcs['obj']
    #     print 'function calls: ', nCalls
    #     print 'boundary constraint: ', boundary
    #     print 'separation constraint: ', separation
    #
    #     # if separation > -1.E-4 and boundary > -1.E-4:
    #     #     print 'AEP opt: ', -funcs['obj']
    #     #
    #     #     file = open('%s/AEP.txt'%folder, 'a')
    #     #     file.write('%s'%(-funcs['obj']) + '\n')
    #     #     file.close()
    #     #
    #     #     file = open('%s/funcCalls.txt'%folder, 'a')
    #     #     file.write('%s'%nCalls + '\n')
    #     #     file.close()
