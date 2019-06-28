from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
import os
from aep_calc import *
from var_reduction_exact import *
import param_fortran
import constraints
import sys
sys.path.insert(0, '/home/flowlab/PJ/reduction/')
from windRoses import *
from position_constraints import calculate_boundary
sys.dont_write_bytecode = True


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


class get_turbine_locs(Component):
    """get turbine locations from grid variables"""

    def __init__(self, nTurbines, nRows, nOuter, nGrid, nBoundaries):

        super(get_turbine_locs, self).__init__()

        self.nTurbines = nTurbines
        self.nOuter = nOuter
        self.nGrid = nGrid
        self.nRows = nRows
        self.nBoundaries = nBoundaries

        # Explicitly size input arrays
        self.add_param('dx', val=0.)
        self.add_param('dy', val=0.)
        self.add_param('shear', val=0.)
        self.add_param('rotate', val=0.)
        self.add_param('start', val=0.)
        self.add_param('y0', val=0.)

        self.add_param('turbs_per_row', val=np.zeros(nRows,dtype=int))
        self.add_param('x_start', val=np.zeros(nRows))

        self.add_param('bx', val=np.zeros(nBoundaries))
        self.add_param('by', val=np.zeros(nBoundaries))

        self.add_output('turbineX', val=np.zeros(nTurbines), pass_by_object=True)
        self.add_output('turbineY', val=np.zeros(nTurbines), pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        dx = params['dx']
        dy = params['dy']
        shear = params['shear']
        rotate = params['rotate']
        start = params['start']
        y0 = params['y0']

        turbs_per_row = params['turbs_per_row']
        x_start = params['x_start']

        bx = params['bx']
        by = params['by']

        dxd = 1.
        dyd = 0.
        sheard = 0.
        rotated = 0.
        startd = 0.
        turbineX,dx_ddx,turbineY,dy_ddx = param_fortran.turbinelocations_dv(self.nTurbines,self.nOuter,self.nGrid,dx,dxd,
                                    dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                                    start,startd,bx,by)
        dxd = 0.
        dyd = 1.
        sheard = 0.
        rotated = 0.
        startd = 0.
        turbineX,dx_ddy,turbineY,dy_ddy = param_fortran.turbinelocations_dv(self.nTurbines,self.nOuter,self.nGrid,dx,dxd,
                                    dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                                    start,startd,bx,by)

        dxd = 0.
        dyd = 0.
        sheard = 1.
        rotated = 0.
        startd = 0.
        turbineX,dx_dshear,turbineY,dy_dshear = param_fortran.turbinelocations_dv(self.nTurbines,self.nOuter,self.nGrid,dx,dxd,
                                    dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                                    start,startd,bx,by)

        dxd = 0.
        dyd = 0.
        sheard = 0.
        rotated = 1.
        startd = 0.
        turbineX,dx_drotate,turbineY,dy_drotate = param_fortran.turbinelocations_dv(self.nTurbines,self.nOuter,self.nGrid,dx,dxd,
                                    dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                                    start,startd,bx,by)

        dxd = 0.
        dyd = 0.
        sheard = 0.
        rotated = 0.
        startd = 1.
        turbineX,dx_dstart,turbineY,dy_dstart = param_fortran.turbinelocations_dv(self.nTurbines,self.nOuter,self.nGrid,dx,dxd,
                                    dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0,
                                    start,startd,bx,by)

        self.dx_ddx = dx_ddx
        self.dy_ddx = dy_ddx
        self.dx_ddy = dx_ddy
        self.dy_ddy = dy_ddy
        self.dx_dshear = dx_dshear
        self.dy_dshear = dy_dshear
        self.dx_drotate = dx_drotate
        self.dy_drotate = dy_drotate
        self.dx_dstart = dx_dstart
        self.dy_dstart = dy_dstart

        unknowns['turbineX'] = turbineX
        unknowns['turbineY'] = turbineY

    def linearize(self, params, unknowns, resids):

        # initialize Jacobian dict
        J = {}

        # populate Jacobian dict
        J[('turbineX', 'dx')] = self.dx_ddx.T
        J[('turbineY', 'dx')] = self.dy_ddx.T

        J[('turbineX', 'dy')] = self.dx_ddy.T
        J[('turbineY', 'dy')] = self.dy_ddy.T

        J[('turbineX', 'shear')] = self.dx_dshear.T
        J[('turbineY', 'shear')] = self.dy_dshear.T

        J[('turbineX', 'rotate')] = self.dx_drotate.T
        J[('turbineY', 'rotate')] = self.dy_drotate.T

        J[('turbineX', 'start')] = self.dx_dstart.T
        J[('turbineY', 'start')] = self.dy_dstart.T

        return J


if __name__ == "__main__":
    global nCalls_obj
    global nCalls_con

    run = 13
    nTurbines = 100
    frac = 0.45
    rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    spacing = 6.
    boundary = 'amalia'
    # boundary = 'circle'
    # boundary = 'square'

    folder = 'opt-results/param/%s_%s_%s'%(boundary,spacing,rose)
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
        nBounds = 100
        # circle_radius = 5280.
        xBounds = np.zeros(nBounds)
        yBounds = np.zeros(nBounds)
        theta = np.linspace(0.,2.*np.pi-2.*np.pi/float(nBounds),nBounds)
        for i in range(nBounds):
            xBounds[i] = circle_radius*np.cos(theta[i])
            yBounds[i] = circle_radius*np.sin(theta[i])

        locations = np.zeros((nBounds,2))
        locations[:, 0] = xBounds
        locations[:, 1] = yBounds
        boundaryVertices, boundaryNormals = calculate_boundary(locations)
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

    bx = boundaryVertices[:,0]
    by = boundaryVertices[:,1]
    bx = bx-(min(bx)+max(bx))/2.
    by = by-(min(by)+max(by))/2.
    bx = np.append(bx,bx[0])
    by = np.append(by,by[0])
    nBounds_con = np.shape(boundaryVertices)[0]
    nBounds_outer = len(bx)

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


    """determine how many turbines will be in the grid here"""
    lenBound = np.zeros(nBounds_outer-1)
    for i in range(nBounds_outer-1):
        lenBound[i] = np.sqrt((bx[i+1]-bx[i])**2+(by[i+1]-by[i])**2)


    circumference = sum(lenBound)
    nOuter = int(frac*nTurbines)
    spacing_outer = circumference/(nOuter*rotorDiameter[0])
    if float(spacing_outer) < 2.*np.sqrt(2.):
        nOuter = int(circumference/(2.*np.sqrt(2.)*rotorDiameter[0]))
    nGrid = int(nTurbines-nOuter)

    print 'nTurbines: ', nTurbines
    print 'nOuter: ', nOuter
    print 'nGrid: ', nGrid

    num = 100
    for i in range(num):
        rotate = float(np.random.rand(1))*360.
        print 'rotate: ', rotate

        dx_start,dy_start,shear_start,rotate_start,turbs_per_row,x_start,y0 = startGrid(nGrid,
                            boundaryVertices*0.85,boundaryNormals,rotate=rotate,my=4.0)

        nRows = len(turbs_per_row)
        print turbs_per_row
        print x_start
        prob = Problem()
        root = prob.root = Group()
        root.add('dx', IndepVarComp('dx', 0.), promotes=['*'])
        root.add('dy', IndepVarComp('dy', 0.), promotes=['*'])
        root.add('shear', IndepVarComp('shear', 0.), promotes=['*'])
        root.add('rotate', IndepVarComp('rotate', 0.), promotes=['*'])
        root.add('start', IndepVarComp('start', 0.), promotes=['*'])

        root.add('get_turbine_locs', get_turbine_locs(nTurbines=nTurbines,nRows=nRows, nOuter=nOuter,
                                        nGrid=nGrid, nBoundaries=nBounds_outer),promotes=['*'])
        root.add('AEP_obj', AEP_calc(nTurbines=nTurbines, nDirections=nDirections*nSpeeds), promotes=['*'])
        root.add('position_constraints', constraint_calc(nTurbines=nTurbines, nBoundaries=nBounds_con),promotes=['*'])

        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SNOPT'
        # prob.driver.opt_settings['Major iterations limit'] = 1000
        prob.driver.opt_settings['Major optimality tolerance'] = 1.E-4
        # prob.driver.opt_settings['Major feasibility tolerance'] = 1.E-4
        prob.driver.opt_settings['Verify level'] = 0
        prob.driver.opt_settings['Scale option'] = 2
        prob.driver.opt_settings['Scale tolerance'] = .95

        prob.driver.opt_settings['Summary file'] = '%s/%s_summary.out'%(folder,run)
        prob.driver.opt_settings['Print file'] = '%s/%s_print.out'%(folder,run)

        prob.driver.add_objective('negAEP')

        prob.driver.add_desvar('dx')
        prob.driver.add_desvar('dy')
        prob.driver.add_desvar('shear')
        prob.driver.add_desvar('rotate')
        prob.driver.add_desvar('start')

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

        prob['turbs_per_row'] = turbs_per_row
        prob['x_start'] = x_start
        prob['y0'] = y0

        prob['bx'] = bx
        prob['by'] = by

        print 'iteration: ', i
        nCalls_obj = 0
        nCalls_con = 0

        prob['dx'] = dx_start + float(np.random.rand(1))*dx_start/5.-dx_start/10.
        prob['dy'] = dy_start + float(np.random.rand(1))*dy_start/5.-dy_start/10.
        prob['shear'] = shear_start + float(np.random.rand(1))*shear_start/5.-shear_start/10.
        prob['rotate'] = rotate_start
        prob['start'] = float(np.random.rand(1))*rotorDiameter[0]*3.

        prob.run()

        separation = np.min(prob['spacing_constraint'])
        boundary = np.min(prob['boundary_constraint'])
        AEP = -prob['negAEP']

        print 'AEP opt: ', AEP
        print 'function calls obj: ', nCalls_obj
        print 'function calls con: ', nCalls_con
        print 'boundary constraint: ', boundary
        print 'separation constraint: ', separation

        # plt.figure(1)
        # plt.clf()
        # for i in range(nTurbines):
        #     circ = plt.Circle((prob['turbineX'][i],prob['turbineY'][i]),prob['rotorDiameter'][i]/2.)
        #     plt.gca().add_patch(circ)
        # bx = prob['boundaryVertices'][:,0]
        # by = prob['boundaryVertices'][:,1]
        # bx = np.append(bx,bx[0])
        # by = np.append(by,by[0])
        # plt.plot(bx,by,'--k')
        # plt.axis('equal')
        # plt.pause(0.001)

        if separation > -1.E-4 and boundary > -1.E-4:
            print 'feasible AEP opt: ', AEP

            file = open('%s/AEP.txt'%folder, 'a')
            file.write('%s'%(AEP) + '\n')
            file.close()

            file = open('%s/funcCalls.txt'%folder, 'a')
            file.write('%s'%nCalls_obj + '\n')
            file.close()
