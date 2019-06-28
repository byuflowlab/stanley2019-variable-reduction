from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
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


class AEP_calc(Component):
    """ Calculates aep with simple gaussian wake model """

    def __init__(self, nTurbines, nDirections):

        # print 'entering windframe __init__ - analytic'

        super(AEP_calc, self).__init__()

        self.nTurbines = nTurbines

        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines))
        self.add_param('turbineY', val=np.zeros(nTurbines))
        self.add_param('turbineZ', val=np.zeros(nTurbines))
        self.add_param('rotorDiameter', val=np.zeros(nTurbines))

        self.add_param('windDirections', val=np.zeros(nDirections))
        self.add_param('windSpeeds', val=np.zeros(nDirections))
        self.add_param('windFrequencies', val=np.zeros(nDirections))


        self.add_output('negAEP', val=0.0, pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        global nCalls_obj
        nCalls_obj += 1

        turbineX = params['turbineX']
        turbineY = params['turbineY']
        turbineZ = params['turbineZ']
        rotorDiameter = params['rotorDiameter']

        windDirections = params['windDirections']
        windSpeeds = params['windSpeeds']
        windFrequencies = params['windFrequencies']

        AEP, daep_dx, daep_dy = fast_calc_gradients(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
                    windSpeeds, windFrequencies)

        self.daep_dx = daep_dx
        self.daep_dy = daep_dy
        unknowns['negAEP'] = -AEP

    def linearize(self, params, unknowns, resids):

        # initialize Jacobian dict
        J = {}

        d_dx = np.zeros((1,self.nTurbines))
        d_dy = np.zeros((1,self.nTurbines))
        for i in range(self.nTurbines):
            d_dx[0][i] = -self.daep_dx[i]
            d_dy[0][i] = -self.daep_dy[i]
        # populate Jacobian dict
        J[('negAEP', 'turbineX')] = d_dx
        J[('negAEP', 'turbineY')] = d_dy

        return J


class constraint_calc(Component):
    """ Calculates spacing and boundary constraints"""

    def __init__(self, nTurbines, nBoundaries):

        super(constraint_calc, self).__init__()

        self.nTurbines = nTurbines
        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines))
        self.add_param('turbineY', val=np.zeros(nTurbines))
        self.add_param('rotorDiameter', val=np.zeros(nTurbines))

        self.add_param('boundaryVertices', val=np.zeros((nBoundaries,2)))
        self.add_param('boundaryNormals', val=np.zeros((nBoundaries,2)))

        self.add_output('spacing_constraint', val=np.zeros((nTurbines-1)*nTurbines/2), pass_by_object=True)
        self.add_output('boundary_constraint', val=np.zeros(nTurbines), pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        global nCalls_con
        nCalls_con += 1

        turbineX = params['turbineX']
        turbineY = params['turbineY']
        rotorDiameter = params['rotorDiameter']

        boundaryVertices = params['boundaryVertices']
        boundaryNormals = params['boundaryNormals']

        dx = np.eye(self.nTurbines)
        dy = np.zeros((self.nTurbines,self.nTurbines))
        _,ss_dx,_,bd_dx = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                                boundaryVertices,boundaryNormals)

        dx = np.zeros((self.nTurbines,self.nTurbines))
        dy = np.eye(self.nTurbines)
        ss,ss_dy,bd,bd_dy = constraints.constraints_position_dv(turbineX,dx,turbineY,dy,
                                boundaryVertices,boundaryNormals)

        bounds = np.zeros(nTurbines)
        index = np.zeros(nTurbines)
        for i in range(nTurbines):
            bounds[i] = np.min(bd[i])
            index[i] = np.argmin(bd[i])

        self.index = index
        self.ss_dx = ss_dx
        self.ss_dy = ss_dy
        self.bd_dx = bd_dx
        self.bd_dy = bd_dy

        unknowns['spacing_constraint'] = ss-(2.*rotorDiameter[0])**2
        unknowns['boundary_constraint'] = bounds

    def linearize(self, params, unknowns, resids):

        # initialize Jacobian dict
        J = {}

        # populate Jacobian dict
        J[('spacing_constraint', 'turbineX')] = self.ss_dx.T
        J[('spacing_constraint', 'turbineY')] = self.ss_dy.T

        db_dx = np.zeros((self.nTurbines,self.nTurbines))
        db_dy = np.zeros((self.nTurbines,self.nTurbines))
        for i in range(nTurbines):
            db_dx[i][i] = self.bd_dx[i][i][self.index[i]]
            db_dy[i][i] = self.bd_dy[i][i][self.index[i]]
        J[('boundary_constraint','turbineX')] = db_dx
        J[('boundary_constraint','turbineY')] = db_dy

        return J


if __name__ == "__main__":
    global nCalls_obj
    global nCalls_con

    run = 20
    nTurbines = 100
    rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    spacing = 4.
    boundary = 'amalia'
    # boundary = 'circle'
    # boundary = 'square'

    folder = 'opt-results/direct/%s_%s_%s'%(boundary,spacing,rose)
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

    rotorDiameter = np.ones(nTurbines)*130.
    turbineZ = np.ones(nTurbines)*110.

    minSpacing = 2.0

    side_length = (np.sqrt(nTurbines)-1.)*rotorDiameter[0]*spacing
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

    nBoundaries = np.shape(boundaryVertices)[0]

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

    prob = Problem()
    root = prob.root = Group()
    root.add('turbineX', IndepVarComp('turbineX', np.zeros(nTurbines)), promotes=['*'])
    root.add('turbineY', IndepVarComp('turbineY', np.zeros(nTurbines)), promotes=['*'])

    root.add('AEP_obj', AEP_calc(nTurbines=nTurbines, nDirections=nDirections*nSpeeds), promotes=['*'])
    root.add('position_constraints', constraint_calc(nTurbines=nTurbines, nBoundaries=nBoundaries),promotes=['*'])

    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    # prob.driver.opt_settings['Major iterations limit'] = 1000
    prob.driver.opt_settings['Major optimality tolerance'] = 1.E-4
    # prob.driver.opt_settings['Major feasibility tolerance'] = 1.E-4
    prob.driver.opt_settings['Verify level'] = 3
    prob.driver.opt_settings['Scale option'] = 2
    prob.driver.opt_settings['Scale tolerance'] = .95

    prob.driver.opt_settings['Summary file'] = '%s/%s_summary.out'%(folder,run)
    prob.driver.opt_settings['Print file'] = '%s/%s_print.out'%(folder,run)

    prob.driver.add_objective('negAEP')

    prob.driver.add_desvar('turbineX')
    prob.driver.add_desvar('turbineY')

    prob.driver.add_constraint('spacing_constraint',lower=0.)
    prob.driver.add_constraint('boundary_constraint',lower=0.)

    prob.setup(check=True)

    prob['turbineZ'] = turbineZ
    prob['rotorDiameter'] = rotorDiameter

    prob['windDirections'] = windDirections
    prob['windSpeeds'] = windSpeeds
    prob['windFrequencies'] = windFrequencies

    prob['boundaryVertices'] = boundaryVertices
    prob['boundaryNormals'] = boundaryNormals


    num = 500
    for i in range(num):
        print 'iteration: ', i
        x,y = random_start(nTurbines,rotorDiameter[0],xmin,xmax,ymin,ymax)
        nCalls_obj = 0
        nCalls_con = 0

        prob['turbineX'] = x
        prob['turbineY'] = y

        prob.run()

        separation = np.min(prob['spacing_constraint'])
        boundary = np.min(prob['boundary_constraint'])
        AEP = -prob['negAEP']

        print 'AEP opt: ', AEP
        print 'function calls obj: ', nCalls_obj
        print 'function calls con: ', nCalls_con
        print 'boundary constraint: ', boundary
        print 'separation constraint: ', separation

        # if separation > -1.E-4 and boundary > -1.E-4:
        #     print 'feasible AEP opt: ', AEP
        #
        #     file = open('%s/AEP.txt'%folder, 'a')
        #     file.write('%s'%(AEP) + '\n')
        #     file.close()
        #
        #     file = open('%s/funcCalls.txt'%folder, 'a')
        #     file.write('%s'%nCalls_obj + '\n')
        #     file.close()
