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
    # global nCalls_obj
    # global nCalls_con
    #
    #
    #
    #
    # run = 4
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


    frac = 0.45
    circumference = sum(lenBound)
    nOuter = int(frac*nTurbines)
    spacing_outer = circumference/(nOuter*rotorDiameter[0])
    print spacing_outer
    if float(spacing_outer) < 2.*np.sqrt(2.):
        nOuter = int(circumference/(2.*np.sqrt(2.)*rotorDiameter[0]))
    nGrid = int(nTurbines-nOuter)

    print 'nTurbines: ', nTurbines
    print 'nOuter: ', nOuter
    print 'nGrid: ', nGrid


    rotate = 0.
    # dx_start,dy_start,shear_start,rotate_start,turbs_per_row,x_start,y0 = startGrid(nGrid,
    #                     boundaryVertices*0.85,boundaryNormals,rotate=rotate,my=4.0)
    # print dx_start/rotorDiameter[0]
    # print dy_start/rotorDiameter[0]
    # nRows = len(turbs_per_row)
    # prob = Problem()
    # root = prob.root = Group()
    # root.add('dx', IndepVarComp('dx', 0.), promotes=['*'])
    # root.add('dy', IndepVarComp('dy', 0.), promotes=['*'])
    # root.add('shear', IndepVarComp('shear', 0.), promotes=['*'])
    # root.add('rotate', IndepVarComp('rotate', 0.), promotes=['*'])
    # root.add('start', IndepVarComp('start', 0.), promotes=['*'])
    #
    # root.add('get_turbine_locs', get_turbine_locs(nTurbines=nTurbines,nRows=nRows, nOuter=nOuter,
    #                                 nGrid=nGrid, nBoundaries=nBounds_outer),promotes=['*'])
    # root.add('AEP_obj', AEP_calc(nTurbines=nTurbines, nDirections=nDirections*nSpeeds), promotes=['*'])
    #
    # prob.setup(check=True)
    #
    # prob['turbineZ'] = turbineZ
    # prob['rotorDiameter'] = rotorDiameter
    #
    # prob['windDirections'] = windDirections
    # prob['windSpeeds'] = windSpeeds
    # prob['windFrequencies'] = windFrequencies
    #
    # prob['turbs_per_row'] = turbs_per_row
    # prob['x_start'] = x_start
    # prob['y0'] = y0
    #
    # prob['bx'] = bx
    # prob['by'] = by
    #
    # print 'iteration: ', i
    # nCalls_obj = 0
    # nCalls_con = 0

    num = 200
    d_dx = np.linspace(2.*rotorDiameter[0],4.*rotorDiameter[0],num)
    d_dy = np.linspace(7.*rotorDiameter[0],9.*rotorDiameter[0],num)
    d_shear = np.linspace(-2.*rotorDiameter[0],2.*rotorDiameter[0],num)
    d_rot = np.linspace(-30.,30.,num)
    d_start = np.linspace(-1.5*rotorDiameter[0],1.5*rotorDiameter[0],num)
    # AEP_dy = np.zeros(num)
    # prob['dx'] = dx_start
    # prob['dy'] = dy_start
    # prob['shear'] = shear_start
    # prob['rotate'] = rotate_start
    # prob['start'] = 0.
    # for i in range(num):
    #         print i
    #         prob['start'] = d_dy[i]
    #         prob.run()
    #         AEP_dy[i] = -prob['negAEP']
    # print 'AEP_dy = np.', repr(AEP_dy)


    AEP_dx = np.array([872013.77177027, 872097.9712696 , 872135.10907682, 872158.17025491,
       872161.09712706, 872191.50498914, 872239.58959301, 872239.55357416,
       872200.32060321, 872116.94474279, 872058.18386638, 872023.61911056,
       872020.717229  , 872035.19221311, 872063.15048729, 872110.74891972,
       872177.80170844, 872217.30869576, 872202.95446425, 872134.55304576,
       871980.21524026, 871804.1488036 , 871596.53472124, 871270.38039426,
       870930.34925711, 870606.8915653 , 870288.71393352, 869918.9514929 ,
       869564.08323191, 869263.65594723, 869026.68482608, 868715.15470098,
       868289.77553671, 867823.24707929, 867322.99540916, 866667.48455893,
       865958.30251514, 865113.99910316, 864055.56848943, 862789.70478458,
       861493.52459549, 860237.92178724, 859009.72273284, 857740.90991   ,
       856724.69038897, 855909.95801353, 855454.18204325, 855009.37724546,
       854609.30401324, 854258.44858512, 853966.59424969, 853838.01542408,
       853665.16394116, 853638.28355374, 853391.29809797, 853214.83797938,
       853245.89378424, 853355.16520284, 853573.10067412, 853844.03141942,
       854411.44595287, 855113.44766342, 855684.43538105, 856325.13180838,
       856966.15251958, 857653.29624854, 858153.70449236, 858470.36415653,
       858664.5120081 , 858772.94421802, 858817.03960112, 858792.90148445,
       858812.83112288, 858793.06013815, 858702.86999464, 858586.40305233,
       858424.98333315, 858095.06933839, 857493.14284635, 856545.15823422,
       855713.91436821, 854991.58625937, 854250.94963085, 853472.84797652,
       852716.1983296 , 852065.53312541, 851513.51806606, 851238.73485983,
       850861.27739076, 850500.37481926, 850275.82033805, 850335.85225092,
       850541.44812089, 850926.98608588, 851333.22800558, 852053.43703657,
       852898.60824725, 853776.29441296, 854723.70065669, 855532.25706068,
       856306.64983104, 856975.67930712, 857683.86857301, 858376.02136006,
       859019.12286362, 859793.07947072, 860638.31035682, 861491.42587122,
       862459.54286441, 863386.99363226, 864159.3446376 , 864789.06791403,
       865227.21777805, 865640.02591345, 865963.35210569, 866183.4326627 ,
       866356.05274413, 866479.30703793, 866546.73295503, 866453.8769284 ,
       866455.45833264, 866523.69657767, 866582.27537972, 866660.3410654 ,
       866586.80201974, 866357.93269475, 866158.52773287, 866024.96638451,
       865820.27058867, 865631.32839789, 865387.68069323, 865108.83297107,
       864708.61996016, 864224.71516602, 863676.93776236, 862922.16354616,
       862211.91503775, 861567.30524082, 860877.0859316 , 860162.53325693,
       859545.26802875, 859084.20200575, 858769.96938508, 858473.38290835,
       858413.84893795, 858311.00265432, 858260.05110844, 858331.43844983,
       858436.9284864 , 858777.54790799, 859226.29656268, 859968.25291111,
       860892.80690549, 861958.57349209, 863226.39394221, 864484.75132249,
       865705.46805264, 866911.27069226, 868136.55154708, 869247.20512146,
       870090.92915435, 870904.85227379, 871599.87577231, 872242.55963822,
       872901.87811498, 873565.45032735, 874235.36406264, 874957.73230676,
       875655.04055249, 876293.34402788, 876797.43229336, 877264.63051449,
       877706.4733415 , 878135.99465876, 878578.68082505, 879047.82688077,
       879486.90217716, 879801.76022379, 880067.52228739, 880337.28287423,
       880619.1513073 , 880896.27012738, 881223.44097007, 881629.49365559,
       882038.34778366, 882391.16127155, 882692.02077896, 883037.44474067,
       883406.4483937 , 883751.21713674, 884082.6392115 , 884394.34009763,
       884657.50667987, 884935.9562132 , 885181.21004139, 885315.11242042,
       885359.78906504, 885362.06552253, 885371.08844523, 885391.32409746])

    # AEP_dy = np.array([782291.26031571, 782677.76010315, 783067.69241748, 783460.62566254,
       # 783856.10565685, 784253.24742731, 784651.45886541, 785050.38581619,
       # 785449.76149182, 785849.74210737, 786250.7285979 , 786653.62225489,
       # 787059.54522465, 787468.73338878, 787880.77569057, 788295.2441446 ,
       # 788711.30524281, 789125.24128657, 789534.94748627, 789938.58708433,
       # 790335.18362595, 790723.73880823, 791105.49110558, 791482.23397073,
       # 791854.16609888, 792222.22700674, 792587.79177765, 792951.6145401 ,
       # 793313.60088884, 793674.04393393, 794033.30637793, 794392.57090163,
       # 794753.21513745, 795115.30261246, 795479.12582671, 795844.69346962,
       # 796211.80965876, 796580.17300738, 796949.81746031, 797320.54431564,
       # 797692.38547125, 798065.35588584, 798440.44476151, 798817.74511487,
       # 799197.02585869, 799578.15609227, 799961.16492854, 800346.05605732,
       # 800732.74589287, 801120.89350471, 801509.71507983, 801898.9158562 ,
       # 802288.40729543, 802678.06001277, 803067.48924749, 803456.62718712,
       # 803844.9302408 , 804231.98507012, 804617.69170009, 805001.74069934,
       # 805384.11079001, 805764.58524024, 806143.02137845, 806519.55429524,
       # 806893.98708699, 807265.87754675, 807634.60564864, 808000.22183315,
       # 808362.9039888 , 808722.68131491, 809079.75529361, 809434.76322498,
       # 809788.37568243, 810141.0372933 , 810493.26905094, 810845.1590372 ,
       # 811196.76672762, 811548.18143827, 811899.46037767, 812250.80905695,
       # 812602.74969717, 812955.35883627, 813308.33133948, 813661.66692693,
       # 814015.36930754, 814369.14515533, 814722.77404087, 815076.1980614 ,
       # 815429.42153682, 815782.337893  , 816134.91160866, 816487.28938771,
       # 816839.19085244, 817190.55214277, 817541.86508596, 817893.07165111,
       # 818243.8581286 , 818593.76934028, 818942.49058615, 819289.9208725 ,
       # 819635.54645501, 819977.4272332 , 820314.11770876, 820644.40170781,
       # 820967.41597925, 821282.14753492, 821588.28236314, 821884.89776528,
       # 822171.59701172, 822448.86482814, 822717.294923  , 822978.30499903,
       # 823233.30314728, 823483.34103076, 823728.93799381, 823970.58654488,
       # 824209.5966598 , 824446.80067661, 824682.31442951, 824916.45871363,
       # 825149.30378022, 825380.88224905, 825611.27475175, 825840.60867834,
       # 826069.23606352, 826297.43983427, 826525.59806045, 826754.02934073,
       # 826982.99339966, 827212.6461523 , 827443.14582298, 827674.29990259,
       # 827906.0769874 , 828138.50917554, 828371.80174269, 828606.11729343,
       # 828841.56426013, 829078.20271969, 829316.11704242, 829555.26904414,
       # 829795.64718843, 830037.19305578, 830279.91637583, 830523.84442981,
       # 830769.05115692, 831015.59313245, 831263.46011389, 831512.62840262,
       # 831763.06963945, 832014.7070178 , 832267.35919619, 832521.01348795,
       # 832775.71462201, 833031.47622714, 833288.33273015, 833546.27061816,
       # 833805.26884619, 834065.49130621, 834326.94459296, 834589.24780501,
       # 834851.82084461, 835114.42762975, 835376.71540954, 835638.36058077,
       # 835899.14360875, 836158.76560948, 836416.85459295, 836673.1673366 ,
       # 836927.83869102, 837180.85086011, 837431.45767983, 837679.47835217,
       # 837924.8411329 , 838167.55913556, 838407.98985439, 838646.01892621,
       # 838881.35821294, 839114.15794179, 839344.61292929, 839572.88830821,
       # 839798.70544127, 840021.88143068, 840242.25692774, 840459.20595999,
       # 840671.88992336, 840879.53489988, 841082.03295648, 841279.20945304,
       # 841471.89789272, 841660.26549142, 841844.30500064, 842023.98911594,
       # 842199.39575347, 842370.57098596, 842537.37874875, 842699.81888753,
       # 842857.79391917, 843011.18036403, 843160.12529796, 843304.85538386])

    AEP_dy = np.array([871881.90756159, 871888.65664618, 871895.08096209, 871901.2258427 ,
       871907.17889694, 871912.99901273, 871918.63167876, 871924.03533885,
       871928.74346515, 871932.55709069, 871935.12437476, 871936.15079652,
       871935.76799494, 871933.42736795, 871929.23521088, 871923.26032672,
       871915.04657958, 871905.82065216, 871894.85600743, 871880.13365683,
       871861.41691267, 871837.26555336, 871809.62661466, 871781.21074181,
       871753.55572253, 871727.49525653, 871703.11357414, 871680.47582658,
       871659.84333372, 871641.15322326, 871624.61908816, 871610.59322295,
       871599.37143653, 871591.14953022, 871585.95259016, 871583.64180735,
       871582.97450679, 871583.39128135, 871584.42241831, 871585.26528471,
       871585.77648013, 871586.36176451, 871587.8242126 , 871591.44246714,
       871597.07115249, 871602.08968633, 871606.26652647, 871609.44569818,
       871611.64277182, 871612.86255181, 871613.32418121, 871614.20574071,
       871616.25038801, 871619.65215962, 871624.70510849, 871630.75451153,
       871636.79911682, 871642.26610839, 871647.4203765 , 871652.21173404,
       871655.25426086, 871658.08726207, 871663.23788242, 871670.52454322,
       871680.06535064, 871692.05975169, 871706.39156644, 871723.00358967,
       871741.72486868, 871762.28834409, 871784.55634431, 871808.65829473,
       871834.58310048, 871862.35533485, 871891.56749061, 871922.28512329,
       871954.53423878, 871988.92499498, 872025.32505977, 872063.38671968,
       872102.43842461, 872142.22428643, 872182.42923941, 872222.67222745,
       872263.10430828, 872304.06381945, 872345.20351489, 872386.19692583,
       872427.25699483, 872468.09611069, 872508.21910233, 872547.49350785,
       872585.48977875, 872621.84837135, 872655.78873572, 872686.89382852,
       872715.37362284, 872741.55991935, 872766.74197191, 872791.29617208,
       872815.1111962 , 872838.65387512, 872862.32389652, 872886.23735175,
       872908.72027165, 872928.86525858, 872947.15206869, 872961.03278531,
       872973.57121833, 872984.76831791, 872994.03670988, 873002.25017301,
       873007.99742385, 873010.93835644, 873010.84327121, 873007.88551665,
       873003.90399997, 873000.17259556, 872996.66809831, 872993.71990768,
       872991.49445612, 872989.98107762, 872988.77532309, 872988.01881009,
       872987.88909408, 872988.1414376 , 872988.96713025, 872990.34444295,
       872992.18313873, 872993.29851063, 872991.27777561, 872989.54064074,
       872987.66122482, 872985.90668206, 872983.41281855, 872978.08310908,
       872971.34047413, 872963.64944792, 872955.19936788, 872946.60276263,
       872937.83306464, 872928.37275333, 872917.51829876, 872905.49025199,
       872892.47219194, 872878.85400854, 872864.36792096, 872849.1187518 ,
       872833.33699719, 872816.8901988 , 872799.78718384, 872782.23267743,
       872764.51631681, 872746.43694548, 872728.05605713, 872710.04206939,
       872692.51062011, 872674.56084613, 872655.88653422, 872636.74271644,
       872617.18336388, 872597.57134874, 872578.39547872, 872559.87635143,
       872542.13338166, 872524.96982202, 872506.98578833, 872487.98411613,
       872468.32395465, 872447.24153627, 872423.80257368, 872399.22814299,
       872373.9022862 , 872348.75911289, 872324.24484767, 872300.82769333,
       872278.94781581, 872258.67845921, 872239.80556321, 872221.53031617,
       872202.92793129, 872184.02767875, 872164.84236701, 872145.80137351,
       872127.53570401, 872110.37405369, 872095.71979345, 872087.72067235,
       872081.56111764, 872076.92366316, 872073.43569108, 872070.93787066,
       872068.67347283, 872065.95772373, 872062.76603866, 872059.91693478,
       872058.15567524, 872057.32913353, 872057.50056811, 872057.29813509])

    AEP_shear = np.array([854888.94973399, 854778.89480857, 854581.15380882, 854316.09062374,
       854036.06479524, 853725.28606694, 853349.94268443, 852865.57284355,
       852395.19328695, 851912.86904061, 851399.69957052, 850898.20688247,
       850420.86195154, 849911.6833759 , 849450.88253376, 849070.57053049,
       848674.51558293, 848283.07998212, 848012.69212723, 847711.48794904,
       847430.21960869, 847082.97299325, 847001.31049596, 847012.38283823,
       847090.10528233, 847244.85068406, 847431.14340073, 847769.20083145,
       848140.68782109, 848438.76205815, 848851.83451354, 849319.14893161,
       849747.59999164, 850162.57272217, 850511.86642602, 850905.76931147,
       851349.18270065, 851788.7141425 , 852265.17321272, 852803.40069222,
       853301.26574391, 853778.22217083, 854229.45147832, 854590.92919966,
       854860.10204862, 855096.58700341, 855268.37512209, 855398.48594627,
       855509.31872894, 855655.83424658, 855825.83649541, 855995.52156137,
       856119.44244557, 856133.21296284, 856097.32583638, 856018.83803755,
       855912.19407639, 855747.24942137, 855488.69886874, 855122.51561453,
       854677.62275674, 854223.56892389, 853758.66987226, 853290.30715028,
       852852.33747313, 852399.13745495, 851961.85634839, 851626.98745425,
       851339.15072598, 851075.58484889, 850945.32342413, 850815.62133733,
       850709.40972502, 850521.73083834, 850589.01076564, 850766.90230088,
       850973.19502593, 851197.44655924, 851467.07515269, 851759.47792288,
       852041.05641809, 852332.74580355, 852690.93070476, 853125.61986729,
       853535.55077112, 853966.51262893, 854390.63455768, 854796.78832539,
       855240.86228047, 855728.00988693, 856233.23156379, 856785.70713623,
       857251.85440289, 857679.1932566 , 858049.2267975 , 858377.82829033,
       858672.02487213, 858937.99109291, 859201.52721204, 859452.11322564,
       859691.17307056, 859937.17663341, 860181.58165972, 860398.77000756,
       860578.45725586, 860664.88647445, 860723.54436814, 860782.84216661,
       860850.55770521, 860927.3502428 , 860987.89424304, 860997.32616908,
       861005.00111052, 861049.91716897, 861122.59793165, 861208.97454119,
       861319.25709503, 861437.48625765, 861569.19234955, 861696.77104065,
       861831.18589704, 862006.01812402, 862201.47735478, 862413.21354389,
       862669.30382332, 862962.36161818, 863292.54071448, 863654.94688815,
       864016.93433504, 864343.25238196, 864658.71748375, 864968.66957804,
       865248.48746063, 865426.50741985, 865587.31551692, 865763.46922106,
       865935.22891644, 866121.72653513, 866327.33414787, 866520.40076342,
       866723.40123003, 866951.31478237, 867198.38310995, 867462.02698339,
       867749.88980014, 868056.63216043, 868331.94013494, 868578.57732908,
       868794.35984688, 869009.63854116, 869235.58601963, 869472.2652324 ,
       869715.56597816, 869961.26951568, 870202.01147335, 870423.34246415,
       870606.44554165, 870759.54528541, 870886.92819751, 870993.27397851,
       871064.95822186, 871092.94375347, 871112.53050792, 871076.9646414 ,
       870993.70417194, 870889.44142771, 870770.36604475, 870651.63002103,
       870546.67261658, 870461.82339418, 870399.63280967, 870361.44484766,
       870347.99187512, 870356.74265115, 870387.20553745, 870433.54689114,
       870489.96462197, 870552.69580333, 870620.27844539, 870688.7411848 ,
       870721.83932268, 870730.88317817, 870734.24745501, 870693.68042676,
       870577.35502438, 870422.86165691, 870268.34309896, 870120.36645501,
       869955.98320694, 869790.40225792, 869635.1332222 , 869499.14863443,
       869384.34036914, 869286.07554406, 869202.94365287, 869129.22424506,
       869069.06439745, 869010.93951501, 868954.9160098 , 868911.32687447])

    AEP_rot = np.array([873768.65782628, 873022.63712664, 871831.95135273, 870431.17843137,
       868988.01638073, 867660.53617707, 866551.0242035 , 865613.10618605,
       865049.66561994, 864645.04391848, 864616.47398717, 865158.01152127,
       866025.25434822, 866919.65691068, 867921.73243242, 868686.80651486,
       869277.49338624, 869774.20984535, 870315.75414203, 870724.66321325,
       871247.05418813, 871790.06642786, 872356.69748383, 872846.16155251,
       873387.87828216, 873940.32094928, 874514.48898196, 874778.24834634,
       874812.55374419, 874771.82779494, 874718.99124647, 874779.97624783,
       874925.1467554 , 875071.0167031 , 875043.90748252, 875149.07680445,
       875362.85280517, 875614.64468255, 876040.86040101, 876286.5488708 ,
       876479.15938349, 876691.83610037, 876876.33364587, 877080.82535391,
       877429.40203639, 877519.98854363, 877557.9413307 , 877577.70930255,
       877537.08728778, 877455.60425047, 877101.82128216, 876612.10702866,
       876207.17897849, 875844.96841554, 875470.4925852 , 875157.04961903,
       874818.21740905, 874524.09794043, 874308.92323699, 874160.62953296,
       873948.99614824, 873514.10406754, 873036.89232558, 872464.37934259,
       871889.58369916, 871491.47261885, 871191.8270639 , 871055.09501075,
       871022.27545525, 871035.52951877, 871123.19414014, 871349.88832856,
       871656.4558848 , 871912.30963228, 872168.04757455, 872451.06830018,
       872721.40181023, 872874.56488708, 873096.24354549, 873360.44080252,
       873709.57580347, 874124.61200049, 874555.50279157, 874963.23567615,
       875303.81482575, 875397.62175743, 875448.81750571, 875366.13974072,
       875165.60341718, 874975.92450562, 874849.34686773, 874691.80254275,
       874456.18629397, 874094.6972375 , 873719.14733913, 873360.42677619,
       873035.94333472, 872732.06877393, 872483.65770601, 872221.23248142,
       871989.24048577, 871786.31394185, 871654.62852308, 871517.33106448,
       871369.84626601, 871231.2092207 , 871067.85103037, 870771.47460618,
       870275.86460185, 869736.44652198, 869196.36995513, 868612.19552362,
       868100.57470104, 867743.51487925, 867496.90304623, 867336.47204352,
       867229.57928444, 867121.85483704, 866950.90124086, 866778.34609668,
       866650.26234568, 866544.22720024, 866521.39505478, 866565.52371536,
       866666.54128241, 866876.60092826, 867255.63432016, 867771.73037102,
       868337.05099458, 868883.23996964, 869471.15386083, 870035.90673243,
       870466.15102137, 870812.02117804, 871125.82800088, 871369.18430337,
       871608.76267218, 871819.50594716, 871951.56216041, 872066.32871301,
       872060.72757633, 871964.52647606, 871889.46410108, 871884.33042602,
       871859.65381422, 871823.87985677, 871703.34513565, 871627.66803207,
       871559.85485894, 871473.29867431, 871347.27377061, 871211.47804765,
       871140.51538845, 871055.22589511, 871018.23772914, 871011.5277911 ,
       871004.74795592, 871089.30089694, 871196.35322109, 871217.99551696,
       871007.07396673, 870730.66331429, 870419.41267022, 870073.2599742 ,
       869652.4074736 , 869218.89966305, 868801.57459156, 868383.42463264,
       868060.38475334, 867872.0604492 , 867733.97077057, 867490.35611035,
       867306.37475263, 867253.91175496, 867222.33007501, 867249.90854458,
       867352.0870523 , 867539.22832809, 867718.91716121, 867867.47150854,
       867942.06023116, 868066.04799961, 868254.39538814, 868394.77821193,
       868509.77244328, 868628.68863598, 868754.29487907, 868978.96027858,
       869187.41385861, 869245.4284808 , 869338.38567111, 869455.09613936,
       869539.77383458, 869509.5612271 , 869322.29139943, 868999.1124424 ,
       868571.67265439, 868092.48949756, 867523.28556836, 866962.18950727])

    AEP_start = np.array([871583.12099784, 871584.51390924, 871585.591814  , 871586.42938119,
       871589.67829701, 871597.57559309, 871604.67511439, 871609.68292216,
       871612.46836832, 871613.36968914, 871615.24669755, 871620.03269488,
       871628.18082231, 871637.30520772, 871645.34121297, 871652.61579629,
       871656.6413847 , 871663.80672652, 871675.86939361, 871693.24933402,
       871715.92439195, 871743.50450564, 871775.19539008, 871810.90787041,
       871850.74734325, 871894.25308196, 871941.10875626, 871992.11003773,
       872047.61216649, 872105.97646929, 872165.90894887, 872226.27710223,
       872287.20517448, 872348.87313434, 872410.41591524, 872471.72709591,
       872531.48328703, 872588.82072434, 872642.2015623 , 872689.55227871,
       872731.00821459, 872768.97410733, 872805.40429618, 872840.76471645,
       872876.42104258, 872910.61741774, 872939.86643207, 872962.19380162,
       872980.41717132, 872994.81270029, 873005.9504577 , 873011.0524993 ,
       873009.43585994, 873003.55798756, 872998.08153797, 872993.48647112,
       872990.55085037, 872988.67633512, 872987.89904795, 872988.18462865,
       872989.72308304, 872992.36978722, 872992.01799857, 872989.37576414,
       872986.60803639, 872982.98895552, 872974.26991684, 872962.91731298,
       872950.12564684, 872937.0331046 , 872922.14579745, 872904.35364294,
       872884.53048048, 872863.02780271, 872839.88174561, 872815.38450243,
       872789.45550882, 872762.91270359, 872735.60990426, 872708.45424908,
       872682.04222278, 872654.18806619, 872625.26771115, 872595.83166046,
       872567.3801928 , 872540.58794016, 872514.54649755, 872486.25350116,
       872456.17258291, 872421.63481321, 872384.36625495, 872346.53235117,
       872310.2633055 , 872277.06551979, 872247.43239285, 872219.8774201 ,
       872191.82144217, 872163.11164403, 872134.91627843, 872108.90297053,
       872090.75474142, 872081.08546923, 872074.74101912, 872070.74140682,
       872067.16725003, 872062.48263419, 872058.76026564, 872057.3058295 ,
       872057.49844256, 872055.24954372, 872046.07516391, 872029.1785628 ,
       872009.63013349, 871991.01885381, 871974.64156046, 871959.55966252,
       871946.17125102, 871934.57379764, 871925.96629582, 871920.58193017,
       871915.16135731, 871909.33927598, 871902.11893558, 871894.06517258,
       871884.05505942, 871871.61165697, 871858.20983671, 871846.94754579,
       871837.5518232 , 871828.14889772, 871819.10671078, 871810.51701821,
       871804.03194265, 871802.7530591 , 871806.59840714, 871813.17209055,
       871821.34924369, 871831.41237858, 871844.26943505, 871858.0982434 ,
       871870.8846817 , 871881.93240425, 871890.53857482, 871899.1508912 ,
       871907.32875113, 871915.03265197, 871922.49977152, 871930.743358  ,
       871940.85892139, 871952.31023952, 871964.00756965, 871973.95009774,
       871979.88753934, 871982.92875477, 871983.86782949, 871980.15816905,
       871971.53070852, 871962.00825723, 871953.48665805, 871944.49859603,
       871934.68690523, 871923.77466274, 871912.87416758, 871901.30106795,
       871891.11621426, 871882.98030141, 871877.03847482, 871869.84697382,
       871866.4604707 , 871866.57652809, 871869.68049314, 871875.84906186,
       871885.12580482, 871894.93391369, 871904.08560383, 871912.8630731 ,
       871921.22100584, 871928.64281312, 871933.93942632, 871936.14529754,
       871934.92492541, 871929.35288091, 871919.562262  , 871906.0511027 ,
       871888.33572931, 871861.90509443, 871824.3603199 , 871781.87880947,
       871740.92202417, 871703.66976682, 871670.39587972, 871641.56983998,
       871617.60130722, 871599.60159022, 871588.29507166, 871583.67182581,
       871583.0681342 , 871584.39781063, 871585.538006  , 871586.34525381])

    fig = plt.figure(figsize=[6.5,2.5])
    ax1 = plt.subplot2grid((1, 10), (0, 0),colspan=2)
    ax2 = plt.subplot2grid((1, 10), (0, 2),colspan=2)
    ax3 = plt.subplot2grid((1, 10), (0, 4),colspan=2)
    ax4 = plt.subplot2grid((1, 10), (0, 6),colspan=2)
    ax5 = plt.subplot2grid((1, 10), (0, 8),colspan=2)




    ax1.plot(d_dx,AEP_dx,color='C0',alpha=0.5)
    ax1.set_xticks([2.*rotorDiameter[0],3.*rotorDiameter[0],4.*rotorDiameter[0]])
    ax1.set_xticklabels(['2','3','4'])
    ax1.set_xlabel(r'$dx$ '+'(rotor\ndiameters)')

    ax2.plot(d_dy,AEP_dy,color='C0',alpha=0.5)
    ax2.set_xticks([7.*rotorDiameter[0],8.*rotorDiameter[0],9.*rotorDiameter[0]])
    ax2.set_xticklabels(['7','8','9'])
    ax2.set_xlabel(r'$dy$ '+'(rotor\ndiameters)')

    ax3.plot(d_shear,AEP_shear,color='C0',alpha=0.5)
    ax3.set_xticks([-2.*rotorDiameter[0],0.,2.*rotorDiameter[0]])
    ax3.set_xticklabels(['-2 D','0','2 D'])
    ax3.set_xlabel(r'$b$ '+'(rotor\ndiameters)')

    ax4.plot(d_rot,AEP_rot,color='C0',alpha=0.5)
    ax4.set_xticks([-30.,0.,30.])
    ax4.set_xticklabels(['-30','0','30'])
    ax4.set_xlabel(r'$\theta$ '+'(degrees)')

    ax5.plot(d_start,AEP_start,color='C0',alpha=0.5)
    ax5.set_xticks([-1.5*rotorDiameter[0],0.,1.5*rotorDiameter[0]])
    ax5.set_xticklabels(['-1.5','0','1.5'])
    ax5.set_xlabel(r'$s$ '+'(rotor\ndiameters)')

    ax1.set_xlim(2.*rotorDiameter[0],4.*rotorDiameter[0])
    ax2.set_xlim(7.*rotorDiameter[0],9.*rotorDiameter[0])
    ax3.set_xlim(-2.*rotorDiameter[0],2.*rotorDiameter[0])
    ax4.set_xlim(-30.,30.)
    ax5.set_xlim(-1.5*rotorDiameter[0],1.5*rotorDiameter[0])

    ax1.set_yticks([825000.,850000.,875000.,900000.])
    ax2.set_yticks([825000.,850000.,875000.,900000.])
    ax3.set_yticks([825000.,850000.,875000.,900000.])
    ax4.set_yticks([825000.,850000.,875000.,900000.])
    ax5.set_yticks([825000.,850000.,875000.,900000.])

    ax1.set_yticklabels(['825','850','875','900'])
    ax2.set_yticklabels(['','','',''])
    ax3.set_yticklabels(['','','',''])
    ax4.set_yticklabels(['','','',''])
    ax5.set_yticklabels(['','','',''])

    ax1.set_ylim(825000.,900000.)
    ax2.set_ylim(825000.,900000.)
    ax3.set_ylim(825000.,900000.)
    ax4.set_ylim(825000.,900000.)
    ax5.set_ylim(825000.,900000.)

    ax1.set_ylabel('AEP (GWh)')


    plt.subplots_adjust(top = 0.93, bottom = 0.3, right = 0.96, left = 0.1,
            hspace = 0., wspace = 0.3)

    plt.show()
