from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
import os
from optComponents import *
import full_aep
import sys
from windRoses import *
from position_constraints import calculate_boundary
sys.dont_write_bytecode = True


def sunflower_points(n, alpha=1.0):
    # this function generates n points within a circle in a sunflower seed pattern
    # the code is based on the example found at
    # https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

    def radius(k, n, b):
        if (k + 1) > n - b:
            r = 1. # put on the boundary
        else:
            r = np.sqrt((k + 1.) - 1. / 2.) / np.sqrt(n - (b + 1.) / 2.)  # apply squareroot

        return r

    x = np.zeros(n)
    y = np.zeros(n)

    b = np.round(alpha * np.sqrt(n)) # number of boundary points

    phi = (np.sqrt(5.) + 1.) / 2.  # golden ratio

    for k in np.arange(0, n):

        r = radius(k, n, b)

        theta = 2. * np.pi * (k+1) / phi**2

        x[k] = r * np.cos(theta)
        y[k] = r * np.sin(theta)

    return x, y


class AEP_obj_ncalls(Component):
    """ Calculates aep with a gaussian wake model """

    def __init__(self, nTurbines, nDirections, nRotorPoints, nCtPoints):

        # print 'entering windframe __init__ - analytic'

        super(AEP_obj_ncalls, self).__init__()

        self.nTurbines = nTurbines

        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines))
        self.add_param('turbineY', val=np.zeros(nTurbines))
        self.add_param('turbineZ', val=np.ones(nTurbines)*110.)
        self.add_param('rotorDiameter', val=np.ones(nTurbines)*130.)
        self.add_param('Ct', val=np.ones(nTurbines)*8./9.)
        self.add_param('yaw', val=np.zeros(nTurbines))

        self.add_param('windDirections', val=np.zeros(nDirections))
        self.add_param('windSpeeds', val=np.zeros(nDirections))
        self.add_param('windFrequencies', val=np.zeros(nDirections))

        self.add_param('RotorPointsY', val=np.zeros(nRotorPoints))
        self.add_param('RotorPointsZ', val=np.zeros(nRotorPoints))

        self.add_param('ct_curve_wind_speed', val=np.zeros(nCtPoints))
        self.add_param('ct_curve_ct', val=np.zeros(nCtPoints))

        self.add_param('shear_exp', val=0.1)
        self.add_param('rated_ws', val=9.8)
        self.add_param('rated_power', val=3.6)
        self.add_param('cut_in_speed', val=4.0)
        self.add_param('cut_out_speed', val=25.0)
        self.add_param('zref', val=50.0)
        self.add_param('z0', val=0.0)
        self.add_param('ky', val=0.022)
        self.add_param('kz', val=0.022)
        self.add_param('alpha', val=2.32)
        self.add_param('beta', val=0.154)
        self.add_param('TI', val=0.11)
        self.add_param('relaxationFactor', val=1.0)
        self.add_param('sm_smoothing', val=700.)
        self.add_param('generator_efficiency', val=0.936)

        self.add_param('wake_combination_method', val=1)
        self.add_param('ti_calculation_method', val=2)
        self.add_param('wake_model_version', val=2016)
        self.add_param('interp_type', val=1)
        self.add_param('calc_k_star', val=True)
        self.add_param('print_ti', val=False)
        self.add_param('use_ct_curve', val=True)

        self.add_output('negAEP', val=0.0, pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        nTurbines = self.nTurbines

        global nCalls_obj
        nCalls_obj += 1

        turbineX = params['turbineX']
        turbineY = params['turbineY']
        turbineZ = params['turbineZ']
        rotorDiameter = params['rotorDiameter']
        Ct = params['Ct']
        yaw = params['yaw']

        windDirections = params['windDirections']
        windSpeeds = params['windSpeeds']
        windFrequencies = params['windFrequencies']

        RotorPointsY = params['RotorPointsY']
        RotorPointsZ = params['RotorPointsZ']

        ct_curve_wind_speed = params['ct_curve_wind_speed']
        ct_curve_ct = params['ct_curve_ct']

        shear_exp = params['shear_exp']
        rated_ws = params['rated_ws']
        rated_power = params['rated_power']
        cut_in_speed = params['cut_in_speed']
        cut_out_speed = params['cut_out_speed']
        zref = params['zref']
        z0 = params['z0']
        ky = params['ky']
        kz = params['kz']
        alpha = params['alpha']
        beta = params['beta']
        TI = params['TI']
        relaxationFactor = params['relaxationFactor']
        sm_smoothing = params['sm_smoothing']
        generator_efficiency = params['generator_efficiency']

        wake_combination_method = params['wake_combination_method']
        ti_calculation_method = params['ti_calculation_method']
        wake_model_version = params['wake_model_version']
        interp_type = params['interp_type']
        calc_k_star = params['calc_k_star']
        print_ti = params['print_ti']
        use_ct_curve = params['use_ct_curve']


        turbineYd = np.zeros((nTurbines,nTurbines))
        turbineXd = np.eye(nTurbines)
        _,daep_dx = full_aep.calcaep_dv(turbineX,turbineXd,turbineY,turbineYd,turbineZ,rotorDiameter,
                    Ct,yaw,windDirections,
                    windSpeeds,windFrequencies,shear_exp,rated_ws,rated_power,
                    cut_in_speed,cut_out_speed,zref,z0,ky,kz,alpha,beta,TI,
                    relaxationFactor,RotorPointsY, RotorPointsZ,ct_curve_wind_speed,ct_curve_ct,
                    sm_smoothing,wake_combination_method,ti_calculation_method,wake_model_version,
                    interp_type,calc_k_star,print_ti,use_ct_curve)

        turbineXd = np.zeros((nTurbines,nTurbines))
        turbineYd = np.eye(nTurbines)
        aep,daep_dy = full_aep.calcaep_dv(turbineX,turbineXd,turbineY,turbineYd,turbineZ,rotorDiameter,Ct,yaw,windDirections,
                    windSpeeds,windFrequencies,shear_exp,rated_ws,rated_power,
                    cut_in_speed,cut_out_speed,zref,z0,ky,kz,alpha,beta,TI,
                    relaxationFactor,RotorPointsY, RotorPointsZ,ct_curve_wind_speed,ct_curve_ct,
                    sm_smoothing,wake_combination_method,ti_calculation_method,wake_model_version,
                    interp_type,calc_k_star,print_ti,use_ct_curve)

        self.daep_dx = daep_dx*generator_efficiency
        self.daep_dy = daep_dy*generator_efficiency
        unknowns['negAEP'] = -aep*generator_efficiency

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


if __name__ == "__main__":
    global nCalls_obj
    global nCalls_con

    # run = 1
    nTurbines = 100
    # rose = 'northIslandRose'
    # rose = 'ukiahRose'
    # rose = 'victorvilleRose'
    # spacing = 4.
    # boundary = 'amalia'
    # boundary = 'circle'
    # boundary = 'square'
    spacing = float(sys.argv[1])
    boundary = '%s'%sys.argv[2]
    rose = '%s'%sys.argv[3]
    run = int(sys.argv[4])

    print 'spacing: ', spacing
    print 'boundary: ', boundary
    print 'rose: ', rose
    print 'run: ', run

    folder = 'opt-results-revision/direct/%s_%s_%s'%(boundary,spacing,rose)
    if not os.path.exists(folder):
        os.makedirs(folder)

    # nDirections = 30
    # nSpeeds = 8
    nDirections = 24
    nSpeeds = 5
    if rose == 'northIslandRose':
        windDirections, windFrequencies, windSpeeds = northIslandRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = northIslandRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle

        windDirections_eval, windFrequencies_eval, windSpeeds_eval = northIslandRose(360,nSpeeds=50)
        windDirections_eval -= wind_angle
    elif rose == 'ukiahRose':
        windDirections, windFrequencies, windSpeeds = ukiahRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = ukiahRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle

        windDirections_eval, windFrequencies_eval, windSpeeds_eval = ukiahRose(360,nSpeeds=50)
        windDirections_eval -= wind_angle
    elif rose == 'victorvilleRose':
        windDirections, windFrequencies, windSpeeds = victorvilleRose(nDirections)
        wind_angle = windDirections[np.argmax(windFrequencies*windSpeeds**3)]
        windDirections, windFrequencies, windSpeeds = victorvilleRose(nDirections,nSpeeds=nSpeeds)
        windDirections -= wind_angle

        windDirections_eval, windFrequencies_eval, windSpeeds_eval = victorvilleRose(360,nSpeeds=50)
        windDirections_eval -= wind_angle

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
        locations = np.loadtxt('/fslhome/pjstanle/compute/reduction/layout_amalia.txt')
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

    ct_speeds = np.array([ 0.000001, 0.1       ,  0.60816327,  1.11632653,  1.6244898 ,  2.13265306,
           2.64081633,  3.14897959,  3.65714286,  4.16530612,  4.67346939,
           5.18163265,  5.68979592,  6.19795918,  6.70612245,  7.21428571,
           7.72244898,  8.23061224,  8.73877551,  9.24693878,  9.75510204,
          10.26326531, 10.77142857, 11.27959184, 11.7877551 , 12.29591837,
          12.80408163, 13.3122449 , 13.82040816, 14.32857143, 14.83673469,
          15.34489796, 15.85306122, 16.36122449, 16.86938776, 17.37755102,
          17.88571429, 18.39387755, 18.90204082, 19.41020408, 19.91836735,
          20.42653061, 20.93469388, 21.44285714, 21.95102041, 22.45918367,
          22.96734694, 23.4755102 , 23.98367347, 24.49183673, 25.        ])
    ct_ct = np.array([0.74988552, 0.74988552, 0.74988552, 0.74988552, 0.74988552, 0.74988552,
          0.74988552, 0.74945275, 0.74736838, 0.74578062, 0.74452166,
          0.7432327 , 0.74240891, 0.74171844, 0.74113119, 0.74062551,
          0.7401854 , 0.7397988 , 0.73945643, 0.73915104, 0.71535516,
          0.50902345, 0.42264255, 0.36002829, 0.31616439, 0.27728908,
          0.2449473 , 0.2179915 , 0.19464155, 0.17388996, 0.15676952,
          0.14116089, 0.12769325, 0.11564223, 0.104593  , 0.09546578,
          0.08765315, 0.08043937, 0.07409357, 0.06822311, 0.06322334,
          0.05887784, 0.05481244, 0.05114998, 0.0474271 , 0.04415572,
          0.04104199, 0.0383636 , 0.03582949, 0.03401271, 0.03235028])

    nCtPoints = len(ct_ct)
    RotorPointsYopt = np.array([0, 0, -0.69, 0.69])
    RotorPointsZopt = np.array([0.69, -0.69, 0, 0])
    nRotorPoints_opt = 4

    nRotorPoints_eval = 100
    RotorPointsYeval, RotorPointsZeval = sunflower_points(nRotorPoints_eval)

    nCalls_obj = 0

    prob = Problem()
    root = prob.root = Group()
    root.add('turbineX', IndepVarComp('turbineX', np.zeros(nTurbines)), promotes=['*'])
    root.add('turbineY', IndepVarComp('turbineY', np.zeros(nTurbines)), promotes=['*'])

    root.add('AEP_obj', AEP_obj_ncalls(nTurbines=nTurbines, nDirections=nDirections*nSpeeds,
            nRotorPoints=nRotorPoints_opt, nCtPoints=nCtPoints), promotes=['*'])
    root.add('position_constraints', constraint_calc(nTurbines=nTurbines, nBoundaries=nBoundaries),promotes=['*'])

    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    # prob.driver.opt_settings['Major iterations limit'] = 1000
    prob.driver.opt_settings['Major optimality tolerance'] = 1.E-4
    # prob.driver.opt_settings['Major feasibility tolerance'] = 1.E-4
    prob.driver.opt_settings['Verify level'] = 3
    prob.driver.opt_settings['Scale option'] = 2
    prob.driver.opt_settings['Scale tolerance'] = .95

    prob.driver.opt_settings['Verify level'] = 0

    prob.driver.opt_settings['Summary file'] = '%s/%s_summary.out'%(folder,run)
    prob.driver.opt_settings['Print file'] = '%s/%s_print.out'%(folder,run)

    prob.driver.add_objective('negAEP',scaler=1E-4)

    prob.driver.add_desvar('turbineX')
    prob.driver.add_desvar('turbineY')

    prob.driver.add_constraint('spacing_constraint',lower=0.)
    prob.driver.add_constraint('boundary_constraint',lower=0.)

    prob.setup(check=True)


    prob['turbineZ'] = turbineZ
    prob['rotorDiameter'] = rotorDiameter

    prob['boundaryVertices'] = boundaryVertices
    prob['boundaryNormals'] = boundaryNormals

    prob['RotorPointsY'] = RotorPointsYopt
    prob['RotorPointsZ'] = RotorPointsZopt

    prob['ct_curve_wind_speed'] = ct_speeds
    prob['ct_curve_ct'] = ct_ct

    prob['windDirections'] = windDirections
    prob['windSpeeds'] = windSpeeds
    prob['windFrequencies'] = windFrequencies


    num = 100
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
        AEPopt = -prob['negAEP']

        xopt = prob['turbineX']
        yopt = prob['turbineY']

        """evaluate"""

        Ct = np.ones(nTurbines)*8./9.
        yaw = np.zeros(nTurbines)
        shear_exp = 0.1
        rated_ws = 9.8
        rated_power = 3.6
        cut_in_speed = 4.0
        cut_out_speed = 25.
        zref = 50.
        z0 = 0.
        ky = 0.022
        kz = 0.022
        alpha = 2.32
        beta = 0.154
        TI = 0.11
        relaxationFactor = 1.0
        sm_smoothing = 700.
        generator_efficiency = 0.936
        wake_combination_method = 1
        ti_calculation_method = 2
        wake_model_version = 2016
        interp_type = 1
        calc_k_star = True
        print_ti = False
        use_ct_curve = True

        if separation > -1.E-2 and boundary > -1.E-2:

            AEPtrue = full_aep.calcaep(xopt,yopt,turbineZ,rotorDiameter,Ct,yaw,windDirections_eval,
                        windSpeeds_eval,windFrequencies_eval,shear_exp,rated_ws,rated_power,
                        cut_in_speed,cut_out_speed,zref,z0,ky,kz,alpha,beta,TI,
                        relaxationFactor,RotorPointsYeval, RotorPointsZeval,ct_speeds,ct_ct,
                        sm_smoothing,wake_combination_method,ti_calculation_method,wake_model_version,
                        interp_type,calc_k_star,print_ti,use_ct_curve)

            AEPtrue = AEPtrue*0.936

            file = open('%s/AEPopt.txt'%folder, 'a')
            file.write('%s'%(AEPopt) + '\n')
            file.close()

            file = open('%s/AEPtrue.txt'%folder, 'a')
            file.write('%s'%(AEPtrue) + '\n')
            file.close()

            file = open('%s/funcCalls.txt'%folder, 'a')
            file.write('%s'%nCalls_obj + '\n')
            file.close()
