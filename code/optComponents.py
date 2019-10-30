from openmdao.api import Group, Component
import numpy as np
import full_aep
import grid_param_fortran
import param_fortran
import constraints
import sys
sys.path.insert(0, '/home/flowlab/PJ/older_projects/reduction/')
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


class AEP_obj(Component):
    """ Calculates aep with a gaussian wake model """

    def __init__(self, nTurbines, nDirections, nRotorPoints, nCtPoints):

        # print 'entering windframe __init__ - analytic'

        super(AEP_obj, self).__init__()

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


class constraint_calc(Component):
    """ Calculates spacing and boundary constraints"""

    def __init__(self, nTurbines, nBoundaries):

        super(constraint_calc, self).__init__()

        self.nTurbines = nTurbines
        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines))
        self.add_param('turbineY', val=np.zeros(nTurbines))
        self.add_param('rotorDiameter', val=np.ones(nTurbines)*130.)

        self.add_param('boundaryVertices', val=np.zeros((nBoundaries,2)))
        self.add_param('boundaryNormals', val=np.zeros((nBoundaries,2)))

        # self.add_param('nCalls_con', val=0)

        self.add_output('spacing_constraint', val=np.zeros((nTurbines-1)*nTurbines/2), pass_by_object=True)
        # self.add_output('spacing_constraint', val=np.zeros(nTurbines), pass_by_object=True)
        self.add_output('boundary_constraint', val=np.zeros(nTurbines), pass_by_object=True)


    def solve_nonlinear(self, params, unknowns, resids):

        # params['nCalls_con'] += 1

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

        bounds = np.zeros(self.nTurbines)
        bounds_index = np.zeros(self.nTurbines, dtype=int)
        for i in range(self.nTurbines):
            bounds[i] = np.min(bd[i])
            bounds_index[i] = np.argmin(bd[i])

        # seps = np.zeros(self.nTurbines)
        # seps_index = np.zeros(self.nTurbines)
        # for i in range(self.nTurbines):
        #     seps[i] = np.min(ss[i])
        #     seps_index[i] = np.argmin(ss[i])

        self.bounds_index = bounds_index
        # self.seps_index = seps_index
        self.ss_dx = ss_dx
        self.ss_dy = ss_dy
        self.bd_dx = bd_dx
        self.bd_dy = bd_dy

        unknowns['spacing_constraint'] = ss-(2.*rotorDiameter[0])**2
        # unknowns['spacing_constraint'] = seps-(2.*rotorDiameter[0])**2
        unknowns['boundary_constraint'] = bounds

    def linearize(self, params, unknowns, resids):

        # initialize Jacobian dict
        J = {}

        # populate Jacobian dict

        # ds_dx = np.zeros((self.nTurbines,self.nTurbines))
        # ds_dy = np.zeros((self.nTurbines,self.nTurbines))
        # for i in range(self.nTurbines):
        #     ds_dx[i][i] = self.ss_dx[i][i][self.seps_index[i]]
        #     ds_dy[i][i] = self.ss_dy[i][i][self.seps_index[i]]
        # J[('spacing_constraint', 'turbineX')] = ds_dx
        # J[('spacing_constraint', 'turbineY')] = ds_dy
        J[('spacing_constraint', 'turbineX')] = self.ss_dx.T
        J[('spacing_constraint', 'turbineY')] = self.ss_dy.T

        db_dx = np.zeros((self.nTurbines,self.nTurbines))
        db_dy = np.zeros((self.nTurbines,self.nTurbines))
        for i in range(self.nTurbines):
            db_dx[i][i] = self.bd_dx[i][i][self.bounds_index[i]]
            db_dy[i][i] = self.bd_dy[i][i][self.bounds_index[i]]
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

        nTurbines = self.nTurbines
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
