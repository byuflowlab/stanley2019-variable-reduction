from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# from position_constraints import *
from windRoses import *
# from grid_param import *
# from aep_calc import *
import os
from openmdao.api import pyOptSparseDriver, Problem, ScipyOptimizer, Group, IndepVarComp
from plantenergy.OptimizationGroups import *
from plantenergy.GeneralWindFarmComponents import calculate_boundary, calculate_distance
from plantenergy.gauss import gauss_wrapper, add_gauss_params_IndepVarComps
from plantenergy.floris import floris_wrapper, add_floris_params_IndepVarComps
from plantenergy import config
import sys
sys.dont_write_bytecode = True

import cPickle as pickle

def random_start(N,D,xmin,xmax,ymin,ymax):
    x = np.zeros(N)
    y = np.zeros(N)
    i = 0
    while i < N:
        print i
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


if __name__ == "__main__":
    test = False

    run = 5

    nTurbs = 100
    rose = 'victorvilleRose'
    spacing = 4.
    boundary = 'amalia'

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

    if test == True:
        nTurbs = 50
        windDirections = np.array([0.])
        windFrequencies = np.array([1.])
        windSpeeds = np.array([8.5])
        spacing = 4.
        boundary = 'circle'
        nDirections = 1
        nSpeeds = 1

    hubHeight = np.ones(nTurbs)*110.
    rotorDiameter = np.ones(nTurbs)*130.

    shearExp = 0.15
    minSpacing = 2.0

    side_length = (np.sqrt(nTurbs)-1.)*rotorDiameter[0]*spacing
    a = side_length**2
    circle_radius = np.sqrt(a/np.pi)

    if boundary == 'circle':
        nVertices = 1
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
        nVertices = len(boundaryVertices)
    elif boundary == 'amalia':
        locations = np.loadtxt('layout_amalia.txt')
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
        nVertices = len(boundaryVertices)


    axialInduction = np.zeros(nTurbs)
    generatorEfficiency = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)

    for turbI in range(0, nTurbs):
        axialInduction[turbI] = 1.0/3.0
        generatorEfficiency[turbI] = 1.0
        Ct[turbI] = 4.0 * axialInduction[turbI] * (1.0 - axialInduction[turbI])
        Cp[turbI] = 4.0 * 1.0 / 3.0 * np.power((1 - 1.0 / 3.0), 2)

    # Define flow properties
    air_density = 1.18    # kg/m^3

    wake_model_version = 2016
    wake_combination_method = 1
    ti_calculation_method = 5 #0 with WEC (original paper uses 4)
    sm_smoothing = 700.
    calc_k_star = True #false for 0 ti, true for 4 and 5
    sort = True
    z_ref = hubHeight[0]
    z_0 = 0.
    TI = 0.1
    k_calc = 0.3837 * TI + 0.003678
    shear_exp = 0.15

    config.obj_func_calls_array[:] = 0.
    config.sens_func_calls_array[:] = 0.

    """OpenMDAO"""
    # prob = Problem()
    # root = prob.root = Group()

    # set up problem
    wake_model_options = {'nSamples': 0,
                         'nRotorPoints': 1,
                         'use_ct_curve': False,
                         'interp_type': 1,
                         'differentiable': True}

    prob = Problem(root=OptAEP(nTurbines=nTurbs, nDirections=nDirections*nSpeeds, wake_model=gauss_wrapper, minSpacing=minSpacing,
                                         wake_model_options=wake_model_options, nVertices=nVertices, datasize=0, use_rotor_components=False,
                                         params_IdepVar_func=add_gauss_params_IndepVarComps, differentiable=True,
                                         params_IndepVar_args={'nRotorPoints': 1},rec_func_calls=True))

    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    # prob.driver.opt_settings['Major iterations limit'] = 1000
    prob.driver.opt_settings['Major optimality tolerance'] = 1.E-4
    prob.driver.opt_settings['Major feasibility tolerance'] = 1.E-4
    prob.driver.opt_settings['Verify level'] = -1
    # prob.driver.opt_settings['Scale option'] = 2
    # prob.driver.opt_settings['Scale tolerance'] = .95

    prob.driver.opt_settings['Summary file'] = 'exact_summary__scaled_%s_%s_%s_%s.out'%(boundary,spacing,rose,run)
    prob.driver.opt_settings['Print file'] = 'exact_print_scaled_%s_%s_%s_%s.out'%(boundary,spacing,rose,run)

    prob.driver.add_objective('obj', scaler=1E-3)

    if boundary == 'circle':
        prob.driver.add_desvar('turbineX', scaler=1E1,lower=np.ones(nTurbs)*(-circle_radius),upper=np.ones(nTurbs)*(circle_radius))
        prob.driver.add_desvar('turbineY', scaler=1E1,lower=np.ones(nTurbs)*(-circle_radius),upper=np.ones(nTurbs)*(circle_radius))
    else:
        prob.driver.add_desvar('turbineX', scaler=1E1,lower=np.ones(nTurbs)*np.min(xBounds),upper=np.ones(nTurbs)*np.max(xBounds))
        prob.driver.add_desvar('turbineY', scaler=1E1,lower=np.ones(nTurbs)*np.min(yBounds),upper=np.ones(nTurbs)*np.max(yBounds))

    prob.driver.add_constraint('sc', lower=np.zeros(int(((nTurbs-1.)*nTurbs/2.))),scaler=1E-2,
                            active_tol=(2.*rotorDiameter[0])**2)
    # prob.driver.add_constraint('boundaryDistances', lower=np.zeros((nTurbs,nVertices)),scaler=1E2,
    #                         active_tol=2.*rotorDiameter)
    prob.driver.add_constraint('boundaryDistances', lower=np.zeros(nTurbs*nVertices),scaler=1E-2,
                            active_tol=2.*rotorDiameter[0])


    prob.root.ln_solver.options['mode'] = 'rev'
    prob.root.ln_solver.options['single_voi_relevance_reduction'] = True

    print 'setup'
    prob.setup(check=False)
    print 'setup finished'

    # voi reduction, active set, first test with verify 2 or 3
    # assign values to constant inputs (not design variables)
    prob['hubHeight'] = hubHeight
    for direction_id in range(nDirections*nSpeeds):
        prob['yaw%i' % direction_id] = np.zeros(nTurbs)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generatorEfficiency'] = generatorEfficiency
    prob['windSpeeds'] = windSpeeds
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies
    prob['air_density'] = air_density
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp

    prob['use_power_curve_definition'] = True
    prob['rated_wind_speed'] = np.ones(nTurbs)*10.
    prob['rated_power'] = np.ones(nTurbs)*3350.
    prob['cut_in_speed'] = np.ones(nTurbs)*3.0
    prob['cut_out_speed'] = np.ones(nTurbs)*25.

    prob['model_params:wake_combination_method'] = wake_combination_method
    prob['model_params:ti_calculation_method'] = ti_calculation_method
    prob['model_params:wake_model_version'] = wake_model_version
    prob['model_params:wec_factor'] = 1.0
    prob['model_params:calc_k_star'] = calc_k_star
    prob['model_params:sort'] = sort
    prob['model_params:z_ref'] = z_ref
    prob['model_params:z_0'] = z_0
    prob['model_params:ky'] = k_calc
    prob['model_params:kz'] = k_calc
    prob['model_params:print_ti'] = False
    prob['model_params:shear_exp'] = shear_exp
    prob['model_params:I'] = TI
    prob['model_params:sm_smoothing'] = sm_smoothing

    if boundary == 'circle':
        prob['boundary_center'] = np.array([0.,0.])
        prob['boundary_radius'] = circle_radius
    else:
        prob['boundaryVertices'] = boundaryVertices
        prob['boundaryNormals'] = boundaryNormals


    if test == True:
        num = 1
    else:
        num = 500

    for i in range(num):

        config.obj_func_calls_array[:] = 0.
        config.sens_func_calls_array[:] = 0.
        print 'iteration: ', i
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
        # plt.figure(1)
        # plt.clf()
        # for i in range(nTurbs):
        #     circ = plt.Circle((x[i],y[i]), rotorDiameter[i]/2.,facecolor="blue",edgecolor="blue",alpha=0.2)
        #     plt.gca().add_patch(circ)
        # circ = plt.Circle((0.,0.), circle_radius,facecolor="None",edgecolor="black",alpha=0.8)
        # plt.gca().add_patch(circ)
        # plt.axis('equal')
        # plt.show()

        prob['turbineX'] = x
        prob['turbineY'] = y

        """Optimization"""
        print 'run'
        prob.run()
        print 'run finished'

        opt_AEP = prob['AEP']
        nCalls = config.obj_func_calls_array[0]
        separation = np.min(prob['sc'])
        boundary = np.min(prob['boundaryDistances'])

        if test == True:
            print 'separation: ', separation
            print 'wtSeparationSquared: ', prob['wtSeparationSquared']
            print 'boundary: ', boundary
            print 'AEP: ', opt_AEP
            print 'nCalls: ', nCalls
            print 'turbineX: ', prob['turbineX']
            print 'turbineY: ', prob['turbineY']
            print 'AEP: ', prob['AEP']
            print 'function calls: ', config.obj_func_calls_array[0] + config.sens_func_calls_array[0]

            plt.figure(1)
            plt.clf()
            for i in range(nTurbs):
                circ = plt.Circle((prob['turbineX'][i],prob['turbineY'][i]), rotorDiameter[i]/2.,facecolor="blue",edgecolor="blue",alpha=0.2)
                plt.gca().add_patch(circ)
            circ = plt.Circle((0.,0.), circle_radius,facecolor="None",edgecolor="black",alpha=0.8)
            plt.gca().add_patch(circ)
            plt.axis('equal')
            plt.show()

        else:

            print 'AEP opt: ', opt_AEP
            print 'boundary constraint: ', boundary
            print 'separation constraint: ', separation

            # plt.figure(1)
            # plt.clf()
            # for i in range(nTurbs):
            #     circ = plt.Circle((prob['turbineX'][i],prob['turbineY'][i]), rotorDiameter[i]/2.,facecolor="blue",edgecolor="blue",alpha=0.2)
            #     plt.gca().add_patch(circ)
            # circ = plt.Circle((0.,0.), circle_radius,facecolor="None",edgecolor="black",alpha=0.8)
            # plt.gca().add_patch(circ)
            # plt.axis('equal')
            # plt.show()

            if separation > -1.E-4 and boundary > -1.E-4:
                print 'AEP opt: ', opt_AEP

                file = open('%s/AEP.txt'%folder, 'a')
                file.write('%s'%(opt_AEP) + '\n')
                file.close()

                file = open('%s/funcCalls.txt'%folder, 'a')
                file.write('%s'%nCalls + '\n')
                file.close()
