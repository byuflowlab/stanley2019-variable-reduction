from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
import os
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients')
from aep_calc import *
import param_fortran
from var_reduction_exact import *
import grid_param_fortran
import constraints
import sys
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
    nTurbines = 100
    rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    spacing = 4.
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

    turbineX = np.array([-1055.45573535,  -582.10981314,  -108.76389094,   364.58203126,
           -1356.51783135,  -883.17190915,  -409.82598695,    63.51993525,
             536.86585745,  1010.21177966, -1657.57992736, -1184.23400516,
            -710.88808296,  -237.54216076,   235.80376144,   709.14968365,
            1182.49560585,  1655.84152805, -1958.64202337, -1485.29610117,
           -1011.95017897,  -538.60425677,   -65.25833456,   408.08758764,
             881.43350984,  1354.77943204,  1828.12535424, -1786.35819718,
           -1313.01227498,  -839.66635278,  -366.32043057,   107.02549163,
             580.37141383,  1053.71733603,  1527.06325823,  2000.40918044,
           -2087.42029319, -1614.07437099, -1140.72844879,  -667.38252658,
            -194.03660438,   279.30931782,   752.65524002,  1226.00116223,
            1699.34708443,  2172.69300663, -1915.136467  , -1441.79054479,
            -968.44462259,  -495.09870039,   -21.75277819,   451.59314401,
             924.93906622,  1398.28498842,  1871.63091062, -2216.19856301,
           -1742.8526408 , -1269.5067186 ,  -796.1607964 ,  -322.8148742 ,
             150.531048  ,   623.87697021,  1097.22289241,  1570.56881461,
            2043.91473681, -2043.91473681, -1570.56881461, -1097.22289241,
            -623.87697021,  -150.531048  ,   322.8148742 ,   796.1607964 ,
            1269.5067186 ,  1742.8526408 ,  2216.19856301, -1871.63091062,
           -1398.28498842,  -924.93906622,  -451.59314401,    21.75277819,
             495.09870039,   968.44462259,  1441.79054479,  1915.136467  ,
           -2172.69300663, -1699.34708443, -1226.00116223,  -752.65524002,
            -279.30931782,   194.03660438,   667.38252658,  1140.72844879,
           -2000.40918044, -1527.06325823, -1053.71733603,  -580.37141383,
            -107.02549163,   366.32043057,  -881.43350984,  -408.08758764])
    turbineY = np.array([-2.84007553e+03, -2.84007553e+03, -2.84007553e+03, -2.84007553e+03,
           -2.36672961e+03, -2.36672961e+03, -2.36672961e+03, -2.36672961e+03,
           -2.36672961e+03, -2.36672961e+03, -1.89338369e+03, -1.89338369e+03,
           -1.89338369e+03, -1.89338369e+03, -1.89338369e+03, -1.89338369e+03,
           -1.89338369e+03, -1.89338369e+03, -1.42003777e+03, -1.42003777e+03,
           -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -1.42003777e+03,
           -1.42003777e+03, -1.42003777e+03, -1.42003777e+03, -9.46691844e+02,
           -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
           -9.46691844e+02, -9.46691844e+02, -9.46691844e+02, -9.46691844e+02,
           -4.73345922e+02, -4.73345922e+02, -4.73345922e+02, -4.73345922e+02,
           -4.73345922e+02, -4.73345922e+02, -4.73345922e+02, -4.73345922e+02,
           -4.73345922e+02, -4.73345922e+02, -2.27373675e-13, -2.27373675e-13,
           -2.27373675e-13, -2.27373675e-13, -2.27373675e-13, -2.27373675e-13,
           -2.27373675e-13, -2.27373675e-13, -2.27373675e-13,  4.73345922e+02,
            4.73345922e+02,  4.73345922e+02,  4.73345922e+02,  4.73345922e+02,
            4.73345922e+02,  4.73345922e+02,  4.73345922e+02,  4.73345922e+02,
            4.73345922e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
            9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  9.46691844e+02,
            9.46691844e+02,  9.46691844e+02,  9.46691844e+02,  1.42003777e+03,
            1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
            1.42003777e+03,  1.42003777e+03,  1.42003777e+03,  1.42003777e+03,
            1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
            1.89338369e+03,  1.89338369e+03,  1.89338369e+03,  1.89338369e+03,
            2.36672961e+03,  2.36672961e+03,  2.36672961e+03,  2.36672961e+03,
            2.36672961e+03,  2.36672961e+03,  2.84007553e+03,  2.84007553e+03])

    prob = Problem()
    root = prob.root = Group()
    root.add('turbineX', IndepVarComp('turbineX', turbineX), promotes=['*'])
    root.add('turbineY', IndepVarComp('turbineY', turbineY), promotes=['*'])

    root.add('AEP_obj', AEP_calc(nTurbines=nTurbines, nDirections=nDirections*nSpeeds), promotes=['*'])

    prob.setup(check=True)

    prob['turbineZ'] = turbineZ
    prob['rotorDiameter'] = rotorDiameter

    prob['windDirections'] = windDirections
    prob['windSpeeds'] = windSpeeds
    prob['windFrequencies'] = windFrequencies

    num = 200
    d_x = np.linspace(min(boundaryVertices[:,0]),max(boundaryVertices[:,0]),num)
    d_y = np.linspace(min(boundaryVertices[:,1]),max(boundaryVertices[:,1]),num)

    print min(d_x)
    print max(d_x)
    print min(d_y)
    print max(d_y)
    # plt.plot(turbineX,turbineY,'o')
    # plt.plot(turbineX[50],turbineY[50],'or')
    # plt.show()
    AEP_y = np.zeros(num)
    nCalls_obj = 0
    nCalls_con = 0
    # for i in range(num):
    #         print i
    #         # x = np.zeros_like(turbineX)
    #         # x[:] = turbineX[:]
    #         # x[50] = d_x[i]
    #         y = np.zeros_like(turbineY)
    #         y[:] = turbineY[:]
    #         y[50] = d_y[i]
    #         prob['turbineY'] = y
    #         prob.run()
    #         AEP_y[i] = -prob['negAEP']
    # print 'AEP_y = np.', repr(AEP_y)

    AEP_x = np.array([835046.93570021, 835034.12179033, 835023.95147924, 834923.7585482 ,
       834755.94087527, 834654.77626259, 834687.20235653, 834835.64790684,
       834804.54073169, 834452.3463334 , 834093.59146961, 833828.56796821,
       833641.62598579, 833424.69601867, 832660.1758774 , 830928.65462942,
       828641.28938967, 826869.30316173, 826603.90269111, 827189.2541967 ,
       829031.26968982, 831141.49560383, 832475.72874029, 832890.03080746,
       833077.09848509, 833278.66557531, 833561.40347443, 833843.60039671,
       834064.36049294, 833803.84550739, 833487.44153365, 833211.85894658,
       833013.53043015, 832871.84194938, 832422.75133916, 830928.6210908 ,
       828752.91275437, 826934.0206137 , 826543.3326616 , 826872.33439781,
       828578.49926344, 830625.88201122, 832214.61129315, 832692.36705441,
       832881.70095764, 833065.57845407, 833319.33585491, 833579.90462644,
       833838.82470884, 833655.33359377, 833344.75641553, 833081.71423875,
       832896.78252953, 832795.00312138, 832451.34732224, 831229.45695749,
       829057.57405628, 827198.57909964, 826550.77149455, 826713.59245326,
       828210.78063942, 830214.91363155, 831989.67439108, 832600.29974309,
       832806.53388118, 832962.28773359, 833203.76391472, 833471.9833459 ,
       833764.28933333, 833667.06217356, 833372.84185624, 833106.28713377,
       832919.13172736, 832829.41732475, 832595.87851978, 831574.38546083,
       829490.82561818, 827497.61631804, 826621.27787608, 826699.89576916,
       827932.7356751 , 829910.17658842, 831820.94978816, 832622.68575315,
       832854.08524242, 833042.46005544, 833290.4433309 , 833587.63426004,
       833920.33031561, 833965.63338405, 833760.87555386, 833578.91939745,
       833515.07181733, 833624.95393639, 833848.74877774, 833726.04556499,
       833599.23186133, 833530.71635281, 833479.72130796, 833470.33783972,
       833592.79235726, 833673.99097768, 833736.96917887, 833754.74313798,
       833636.12347897, 833578.34793662, 833627.87162924, 833759.20555828,
       833974.19042915, 834014.35229297, 833691.7882362 , 833360.75511346,
       833089.84132711, 832936.03962289, 832783.69814852, 832081.9388962 ,
       830312.06994137, 828163.06287438, 826664.67311831, 826543.99885727,
       827230.61601532, 829098.68382814, 831151.27562229, 832349.84349043,
       832685.06908256, 832855.01502157, 833019.66467852, 833270.07760876,
       833553.90559105, 833752.15837741, 833497.9444379 , 833225.05564036,
       832997.75406761, 832879.16344111, 832795.08435621, 832245.38431569,
       830728.85012704, 828589.63091338, 826831.15123498, 826588.96011149,
       827043.51484769, 828769.46804626, 830819.84707772, 832290.03869696,
       832725.67866405, 832890.03069576, 833038.90009387, 833292.2216207 ,
       833576.3811969 , 833855.04881654, 833645.01368667, 833374.18595583,
       833135.97977338, 832993.58710521, 832907.79514484, 832490.753733  ,
       831236.16306858, 829054.70084695, 827253.25696941, 826706.35468525,
       826944.64222238, 828615.24957266, 830667.41880647, 832360.38255661,
       832951.49914136, 833133.23598871, 833255.58180236, 833482.69265145,
       833754.93723221, 834066.21185692, 833920.2842682 , 833659.04517589,
       833408.2304877 , 833243.31159936, 833151.27475924, 832865.23822847,
       831817.91632187, 829684.01855569, 827763.91803856, 827000.00241024,
       827157.91476053, 828585.00788877, 830685.37371872, 832618.37782054,
       833432.33147734, 833658.37525113, 833844.91177293, 834151.02648685,
       834520.9702012 , 834955.64151342, 835017.56591524, 834853.26451657,
       834721.73102784, 834707.35167854, 834845.13113734, 835078.83527136,
       834986.17412898, 834916.36381922, 834969.90271865, 835154.79564433])


    AEP_y = np.array([834445.40103313, 834241.19010449, 834026.47904582, 833759.93225976,
       833302.1915122 , 832871.31168452, 832587.02098724, 832500.56094421,
       832604.06068281, 832867.50304557, 833100.3932885 , 833259.25460521,
       833389.88915029, 833344.75165318, 833308.94242001, 833229.05599223,
       833063.10587822, 832903.43672873, 832673.17877426, 832530.79622176,
       832454.51078494, 832401.98437793, 832665.26883361, 832896.41352652,
       833049.2334414 , 833173.82416769, 833312.44687925, 833517.62167408,
       833605.2474247 , 833657.08077988, 833747.98626004, 833763.24986235,
       833828.17404195, 833913.76670456, 833910.42056542, 833851.5353783 ,
       833776.53827919, 833703.22433647, 833691.39041207, 833626.35617352,
       833527.18443696, 833481.78786382, 833356.35418791, 833203.19858851,
       833123.76891819, 832950.1342973 , 832731.30940418, 832449.32796225,
       832086.76762299, 831594.69186117, 830855.35440782, 829546.6280075 ,
       828671.75049599, 828585.90606981, 829499.28836127, 830507.41002165,
       831331.07624175, 831936.36457003, 832238.57620345, 832488.61499076,
       832691.14277429, 832809.93823704, 832827.60429188, 832777.46457828,
       832706.96741574, 832627.6011794 , 832629.55359582, 832713.57620871,
       832884.53245528, 833124.80131512, 833242.4054585 , 833284.50996165,
       833323.71917281, 833360.18276835, 833395.61184368, 833442.5800949 ,
       833418.71289052, 833421.99724239, 833487.76805565, 833565.98064666,
       833642.11457885, 833701.24997366, 833719.57142964, 833616.35535071,
       833481.5930475 , 833431.86919523, 833408.78759096, 833451.75729365,
       833589.99875247, 833691.65183542, 833713.61404177, 833764.61788891,
       833807.10772297, 833853.55860453, 833892.34336669, 833864.87972931,
       833814.26270869, 833680.55493317, 833563.514587  , 833470.75673457,
       833472.28056942, 833566.54225991, 833663.65620541, 833741.76628097,
       833754.8031397 , 833779.68554646, 833803.46239768, 833868.01675233,
       833923.99217372, 833939.44542339, 833948.53638958, 833841.19595277,
       833670.03087679, 833589.52141152, 833502.57599138, 833510.38090058,
       833636.30694803, 833678.75581581, 833605.99422273, 833553.17778916,
       833426.71292276, 833388.78433407, 833372.01541982, 833407.79320952,
       833438.81027885, 833414.02730893, 833386.1981717 , 833399.75770176,
       833400.41635335, 833341.33663138, 833132.48640146, 832918.17618125,
       832738.04084451, 832632.43266444, 832583.96831796, 832657.5448602 ,
       832710.55510498, 832738.40624782, 832719.80581889, 832625.62259631,
       832540.18262334, 832293.11153944, 831990.79982107, 831407.86389176,
       830480.49079945, 829339.95955732, 828529.34130521, 828667.57388457,
       829490.20657267, 830756.29779947, 831519.76227565, 832067.92225055,
       832495.46831317, 832835.51726703, 833082.67875037, 833256.7535687 ,
       833422.67689653, 833633.47270958, 833819.91602567, 833949.97091202,
       834010.52453313, 834002.51399057, 833979.1007088 , 834032.50848009,
       834011.08377693, 834010.25293321, 833972.75617685, 833921.39193577,
       833886.33834206, 833838.26459848, 833775.08183231, 833805.01694919,
       833786.87094535, 833622.63044496, 833477.07213787, 833341.37473006,
       833186.74751983, 832807.86356164, 832644.56910379, 832656.13609923,
       832835.0784661 , 833034.18398107, 833374.85718819, 833671.39947581,
       833889.6980552 , 834109.03708445, 834350.18254657, 834587.61755088,
       834707.12810903, 834782.88203969, 834821.51382438, 834891.41434737,
       835023.58907909, 835138.0722293 , 835205.99217194, 835280.88575041,
       835358.34578171, 835432.48139145, 835484.7544569 , 835584.98976928])
    fig = plt.figure(figsize=[6.5,2.5])
    ax1 = plt.subplot2grid((1, 10), (0, 3), colspan=2)
    ax2 = plt.subplot2grid((1, 10), (0, 5), colspan=2)
    # ax3 = plt.subplot2grid((1, 5), (0, 2))
    # ax4 = plt.subplot2grid((1, 5), (0, 3))
    # ax5 = plt.subplot2grid((1, 5), (0, 4))
    #
    #
    #
    #
    ax1.plot(d_x,AEP_x,color='C3',alpha=0.5)
    # ax1.set_xticks([2.*rotorDiameter[0],3.*rotorDiameter[0],4.*rotorDiameter[0]])
    # ax1.set_xticklabels(['2','3','4'])
    # ax1.set_xlabel(r'$dx$ '+'(rotor\ndiameters)')

    ax2.plot(d_y,AEP_y,color='C3',alpha=0.5)
    # ax2.set_xticks([7.*rotorDiameter[0],8.*rotorDiameter[0],9.*rotorDiameter[0]])
    # ax2.set_xticklabels(['7','8','9'])
    # ax2.set_xlabel(r'$dy$ '+'(rotor\ndiameters)')
    #
    # ax1.set_xlim(2.*rotorDiameter[0],4.*rotorDiameter[0])
    # ax2.set_xlim(7.*rotorDiameter[0],9.*rotorDiameter[0])

    ax1.set_yticks([825000.,829000.,833000.,837000.])
    ax2.set_yticks([825000.,829000.,833000.,837000.])

    #
    # ax1.set_yticklabels(['825','850','875','900'])
    # ax2.set_yticklabels(['','','',''])

    ax1.set_ylim(825000.,837000.)
    ax2.set_ylim(825000.,837000.)

    ax1.set_ylabel('AEP (GWh)')

    plt.subplots_adjust(top = 0.93, bottom = 0.3, right = 0.96, left = 0.1,
            hspace = 0., wspace = 0.3)

    plt.show()
