from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
import os
from aep_calc import *
from var_reduction_exact import *
import grid_param_fortran
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
        # turbineX = turbineX-(max(turbineX)+min(turbineX))/2.
        turbineY = params['turbineY']
        # turbienY = turbineY-(max(turbineY)+min(turbineY))/2.
        rotorDiameter = params['rotorDiameter']

        boundaryVertices = params['boundaryVertices']
        # bx = params['boundaryVertices'][:,0]
        # by = params['boundaryVertices'][:,1]
        # bx = bx + (max(turbineX)+min(turbineX))/2.
        # by = by + (max(turbineY)+min(turbineY))/2.
        # boundaryVertices[:,0] = bx[:]
        # boundaryVertices[:,1] = by[:]
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


class form_grid(Component):
    """get turbine locations from grid variables"""

    def __init__(self, nTurbines, nRows):

        super(form_grid, self).__init__()

        self.nTurbines = nTurbines

        # Explicitly size input arrays
        self.add_param('dx', val=0.)
        self.add_param('dy', val=0.)
        self.add_param('shear', val=0.)
        self.add_param('rotate', val=0.)
        self.add_param('y0', val=0.)

        self.add_param('turbs_per_row', val=np.zeros(nRows,dtype=int))
        self.add_param('x_start', val=np.zeros(nRows))

        self.add_output('turbineX', val=np.zeros(nTurbines), pass_by_object=True)
        self.add_output('turbineY', val=np.zeros(nTurbines), pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        dx = params['dx']
        dy = params['dy']
        shear = params['shear']
        rotate = params['rotate']
        y0 = params['y0']
        turbs_per_row = params['turbs_per_row']
        x_start = params['x_start']

        dxd = 1.
        dyd = 0.
        sheard = 0.
        rotated = 0.
        turbineX,dx_ddx,turbineY,dy_ddx = grid_param_fortran.makegrid_fortran_dv(nTurbines,
                        dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

        dxd = 0.
        dyd = 1.
        sheard = 0.
        rotated = 0.
        turbineX,dx_ddy,turbineY,dy_ddy = grid_param_fortran.makegrid_fortran_dv(nTurbines,
                        dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

        dxd = 0.
        dyd = 0.
        sheard = 1.
        rotated = 0.
        turbineX,dx_dshear,turbineY,dy_dshear = grid_param_fortran.makegrid_fortran_dv(nTurbines,
                        dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

        dxd = 0.
        dyd = 0.
        sheard = 0.
        rotated = 1.
        turbineX,dx_drotate,turbineY,dy_drotate = grid_param_fortran.makegrid_fortran_dv(nTurbines,
                        dx,dxd,dy,dyd,shear,sheard,rotate,rotated,turbs_per_row,x_start,y0)

        self.dx_ddx = dx_ddx
        self.dy_ddx = dy_ddx
        self.dx_ddy = dx_ddy
        self.dy_ddy = dy_ddy
        self.dx_dshear = dx_dshear
        self.dy_dshear = dy_dshear
        self.dx_drotate = dx_drotate
        self.dy_drotate = dy_drotate

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


        return J


if __name__ == "__main__":
    global nCalls_obj
    global nCalls_con

    run = 4
    nTurbines = 100
    rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    spacing = 8.
    boundary = 'amalia'
    # boundary = 'circle'
    # boundary = 'square'

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
        locations = np.loadtxt('/Users/ningrsrch/Dropbox/Projects/reduction/layout_amalia.txt')
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

    num = 1
    for i in range(num):
        rotate = 0.
        dx_start,dy_start,shear_start,rotate_start,turbs_per_row,x_start,y0 = startGrid(nTurbines,
                            boundaryVertices,boundaryNormals,rotate=rotate,my=1.0)

        nRows = len(turbs_per_row)
        prob = Problem()
        root = prob.root = Group()
        root.add('dx', IndepVarComp('dx', 0.), promotes=['*'])
        root.add('dy', IndepVarComp('dy', 0.), promotes=['*'])
        root.add('shear', IndepVarComp('shear', 0.), promotes=['*'])
        root.add('rotate', IndepVarComp('rotate', 0.), promotes=['*'])

        root.add('makeGrid', form_grid(nTurbines=nTurbines,nRows=nRows),promotes=['*'])
        root.add('AEP_obj', AEP_calc(nTurbines=nTurbines, nDirections=nDirections*nSpeeds), promotes=['*'])
        root.add('position_constraints', constraint_calc(nTurbines=nTurbines, nBoundaries=nBoundaries),promotes=['*'])

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

        print 'iteration: ', i
        nCalls_obj = 0
        nCalls_con = 0

        prob['dx'] = dx_start
        prob['dy'] = dy_start
        prob['shear'] = shear_start
        prob['rotate'] = rotate_start

        prob.run()

        separation = np.min(prob['spacing_constraint'])
        boundary = np.min(prob['boundary_constraint'])
        AEP = -prob['negAEP']

        print 'AEP opt: ', AEP
        print 'function calls obj: ', nCalls_obj
        print 'function calls con: ', nCalls_con
        print 'boundary constraint: ', boundary
        print 'separation constraint: ', separation

        print 'turbineX: ', repr(prob['turbineX'])
        print 'turbineY: ', repr(prob['turbineY'])

        plt.figure(1)
        plt.clf()
        for i in range(nTurbines):
            circ = plt.Circle((prob['turbineX'][i],prob['turbineY'][i]),prob['rotorDiameter'][i]/2.)
            plt.gca().add_patch(circ)
        bx = prob['boundaryVertices'][:,0]
        by = prob['boundaryVertices'][:,1]
        bx = np.append(bx,bx[0])
        by = np.append(by,by[0])
        plt.plot(bx,by,'--k')
        plt.axis('equal')
        plt.show()
