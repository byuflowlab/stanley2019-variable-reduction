from pyoptsparse import Optimization, SNOPT, pyOpt_solution, NSGA2
from openmdao.api import Group, Component, Problem, IndepVarComp, pyOptSparseDriver
import numpy as np
import scipy as sp
import os
import sys
sys.path.insert(0, '/Users/ningrsrch/Dropbox/Projects/reduction/exact-gradients')
from aep_calc import *
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

    AEP_dx = np.array([758571.3925761 , 759372.10573822, 760088.16732161, 760730.99624786,
           761317.54518145, 761843.8579298 , 762309.47377813, 762718.66548757,
           763075.87814329, 763382.25940078, 763640.7468331 , 763856.88353666,
           764035.47678197, 764181.22906932, 764297.6851575 , 764389.38101883,
           764461.2356392 , 764518.18219062, 764563.70094388, 764601.10281197,
           764633.3674478 , 764663.47638735, 764695.50840206, 764732.92793811,
           764777.65100966, 764831.39157733, 764896.00788439, 764974.24032095,
           765068.97210194, 765181.23234166, 765308.40604053, 765452.9967108 ,
           765616.26605424, 765795.30333841, 765984.99253261, 766190.77099916,
           766418.23710763, 766669.14222595, 766944.49539416, 767246.52991211,
           767577.35920254, 767939.02711944, 768333.84996461, 768765.04234321,
           769237.4168748 , 769753.23635545, 770314.01038069, 770920.12765146,
           771571.56184228, 772267.87862689, 773008.41585255, 773792.66937365,
           774620.04324629, 775489.95247538, 776401.75238255, 777354.45768818,
           778347.15820033, 779378.85918476, 780448.36192623, 781554.59685428,
           782696.41430879, 783872.50090445, 785081.67222754, 786322.71322285,
           787594.35636051, 788895.33466086, 790224.34324109, 791580.02749336,
           792961.0750319 , 794366.10402616, 795793.57207333, 797241.86890709,
           798709.3121687 , 800193.80398444, 801681.03496405, 803159.08884433,
           804622.93859255, 806061.47592962, 807444.0103751 , 808763.72104362,
           810016.88321802, 811203.56298454, 812346.89635152, 813443.19483319,
           814482.6770299 , 815448.74698822, 816343.54775614, 817173.05558846,
           817944.69819308, 818673.06399443, 819370.47021726, 820047.87047238,
           820710.97661648, 821361.7746969 , 821998.89238913, 822619.69000483,
           823222.49248808, 823805.77532558, 824366.71078854, 824902.87146902,
           825409.51081268, 825884.81478832, 826328.33582082, 826722.92777   ,
           827058.91010457, 827336.69662241, 827558.7893403 , 827730.38432472,
           827871.49766246, 827987.69553828, 828079.92338031, 828149.5151447 ,
           828199.19793051, 828235.23806995, 828259.30430408, 828271.99406416,
           828274.08613788, 828266.90022418, 828251.4224739 , 828228.30839484,
           828198.03476405, 828161.42065027, 828119.39796599, 828072.90287891,
           828022.48271932, 827969.0533218 , 827913.67147856, 827857.080951  ,
           827800.2229523 , 827743.79052862, 827688.93006523, 827637.1317384 ,
           827589.66589025, 827547.77634116, 827512.75494682, 827485.73345763,
           827467.86115884, 827460.26935666, 827464.21865005, 827480.7589441 ,
           827510.86140965, 827555.38865605, 827615.19198115, 827691.38173094,
           827784.99405357, 827896.63534958, 828026.83051136, 828176.03705997,
           828344.67782593, 828532.96282955, 828741.1399415 , 828969.507329  ,
           829218.29646312, 829487.67835633, 829777.76845474, 830088.62662112,
           830420.25707558, 830772.60820639, 831145.57235849, 831538.98776296,
           831952.638475  , 832386.25659889, 832839.52209701, 833312.06424987,
           833803.46355383, 834313.25391748, 834840.92513487, 835385.92560853,
           835947.66344592, 836525.49275198, 837118.72748698, 837726.66207661,
           838348.56288639, 838983.66830919, 839631.14959241, 840290.16389512,
           840959.9537792 , 841639.77797826, 842328.89328863, 843026.57919814,
           843732.13517521, 844444.8496525 , 845163.98208441, 845888.80014092,
           846618.59900835, 847352.69452271, 848090.42812002, 848831.1766485 ,
           849574.32470528, 850319.27465124, 851065.45705934, 851812.21581137,
           852554.39928813, 853288.89177138, 854014.94693326, 854731.86239126,
           855438.97960279, 856132.7296434 , 856801.34336775, 857426.29177933])
    AEP_dy = np.array([736020.81802209, 737046.49890055, 738086.73287083, 739140.08369086,
           740205.22450226, 741281.75010871, 742370.09628149, 743468.97421599,
           744578.01678061, 745697.75678834, 746828.34497151, 747970.09814933,
           749122.5711138 , 750285.21599537, 751459.81849691, 752647.79248249,
           753849.53318381, 755065.50582137, 756295.3616589 , 757538.84127884,
           758796.54195638, 760068.4147818 , 761355.03932556, 762657.10309369,
           763974.5619423 , 765306.35299026, 766651.63908462, 768009.69096114,
           769380.10799101, 770762.93910077, 772157.53885125, 773563.26391135,
           774977.65954018, 776393.24105019, 777807.46602013, 779217.55381378,
           780620.9888032 , 782010.5206962 , 783382.87608976, 784734.79389936,
           786064.08077907, 787367.07737475, 788645.78004949, 789905.78097738,
           791158.28979886, 792403.15727387, 793638.2806534 , 794860.01772233,
           796067.33581719, 797261.28817294, 798443.78048233, 799614.41446759,
           800769.86241415, 801909.60468707, 803034.33584273, 804145.40488384,
           805243.81401916, 806329.53508896, 807401.70885862, 808460.23758928,
           809505.01814816, 810535.92750415, 811553.79887588, 812560.75688061,
           813556.87742022, 814542.01829087, 815516.85153124, 816482.21795321,
           817438.79572647, 818387.28487242, 819327.83354487, 820258.38322252,
           821175.50583383, 822078.79381148, 822967.17957331, 823831.91838903,
           824665.99644879, 825471.99372086, 826248.15604501, 826996.92095305,
           827725.47831479, 828436.67587151, 829129.57245922, 829801.86259268,
           830453.65584904, 831085.48813892, 831697.78929187, 832291.92393442,
           832868.92068766, 833431.21966931, 833980.58112166, 834518.08408569,
           835044.29195742, 835559.69805445, 836060.59945886, 836541.42775257,
           837001.67927283, 837440.94992058, 837857.44891613, 838244.5392461 ,
           838598.05296147, 838918.91476905, 839205.92824192, 839444.53247451,
           839631.08970574, 839763.99064254, 839847.85335715, 839880.52228107,
           839868.30876014, 839826.44836434, 839758.05141546, 839660.21019358,
           839536.40840915, 839390.59590821, 839226.84558772, 839049.54464844,
           838861.69127884, 838665.52181191, 838465.445822  , 838264.81991419,
           838064.510622  , 837865.77507992, 837669.54405313, 837475.9486646 ,
           837284.71895719, 837095.63517498, 836909.48837944, 836726.21397415,
           836545.82805883, 836368.11658232, 836192.87534276, 836019.85056742,
           835849.14287372, 835681.16744974, 835516.31155625, 835354.71809314,
           835196.46177514, 835041.94202476, 834891.66648263, 834746.05656812,
           834605.74125274, 834471.21981831, 834343.08099596, 834221.74621091,
           834107.41743429, 834000.26337344, 833900.44602628, 833808.08675846,
           833723.26785839, 833646.10772363, 833576.72327521, 833515.2303247 ,
           833461.73397367, 833416.37414791, 833379.33671631, 833350.76760131,
           833330.78177715, 833319.49088558, 833316.99238238, 833323.3848873 ,
           833338.75835078, 833363.17520626, 833396.68140806, 833439.31272499,
           833491.09495108, 833552.04316412, 833622.16326545, 833701.45292103,
           833789.8961977 , 833887.4646539 , 833994.12065545, 834110.33743263,
           834236.54380952, 834372.77349119, 834519.03297197, 834674.15937295,
           834837.2068258 , 835008.01114686, 835186.34366346, 835370.87764263,
           835560.23257712, 835752.95029975, 835947.36632917, 836142.3237844 ,
           836337.5209641 , 836534.68186357, 836734.29087386, 836936.65199159,
           837143.070662  , 837354.98660062, 837573.54462443, 837800.02881967,
           838034.92038703, 838278.16213811, 838529.49580292, 838788.81734955,
           839056.04541737, 839331.09199858, 839613.86263551, 839904.25964808])

    AEP_shear = np.array([866614.20791648, 867083.31508584, 867607.69679304, 868174.5506042 ,
           868601.40038776, 868773.82551025, 868756.76134137, 868727.47085432,
           868717.75331526, 868799.30319957, 868989.36476399, 869261.96533619,
           869562.3894521 , 869859.71339524, 870052.31316873, 870086.29322821,
           869955.49152229, 869737.3174354 , 869484.33457198, 869244.15713215,
           869001.35575567, 868751.57519814, 868511.21821923, 868282.47302541,
           868059.43200178, 867789.01691814, 867347.57330904, 866408.08985545,
           864950.01798671, 863317.55039412, 861631.29974855, 859961.17577074,
           858334.18936351, 856767.12936921, 855275.54700604, 853872.90161019,
           852569.97332648, 851375.17813766, 850293.33416896, 849317.41240315,
           848421.35470454, 847603.88139134, 846883.03053211, 846284.02963094,
           845814.556511  , 845479.30926029, 845289.00632737, 845252.59062152,
           845363.82200523, 845587.65938075, 845920.56242634, 846317.4832543 ,
           846656.72148124, 846903.85600497, 847194.98724136, 847631.48093562,
           848227.61083406, 848976.13104177, 849865.59650795, 850883.81380879,
           852014.63928463, 853236.09018763, 854519.11029063, 855838.0787497 ,
           857187.53395122, 858564.21445779, 859962.86859457, 861375.70208355,
           862793.87989206, 864200.1680851 , 865529.7474374 , 866679.7033053 ,
           867348.08121422, 867494.20366501, 867364.76381785, 867106.25967859,
           866761.76602596, 866178.42457937, 864874.14101795, 862750.42856122,
           860361.4681185 , 857897.44264418, 855445.70367739, 853023.66901962,
           850641.04859266, 848299.44593573, 846013.35647149, 843806.53766511,
           841705.58859056, 839721.78945297, 837854.47957625, 836114.37795172,
           834529.4575602 , 833121.63706559, 831903.05032179, 830880.8924476 ,
           830058.12035995, 829432.93594366, 829004.13829758, 828774.47391829,
           828745.37853624, 828914.56421701, 829274.17439011, 829809.1310179 ,
           830522.93565998, 831426.51584773, 832517.9297852 , 833790.42652673,
           835236.80657923, 836848.25236749, 838614.17023679, 840522.22508061,
           842556.81865843, 844694.93628412, 846903.47181903, 849151.24072424,
           851410.29896538, 853672.7354191 , 855930.18636588, 858170.45652632,
           860316.76632907, 862209.69741001, 863388.45147206, 863859.46353613,
           864010.04694957, 864032.35381979, 863935.75414871, 863557.8198174 ,
           862476.9856212 , 860672.03863586, 858611.79608039, 856447.21673484,
           854258.41483717, 852064.1569366 , 849874.15378958, 847699.10159217,
           845540.21994986, 843400.42079447, 841315.1528563 , 839327.73059401,
           837477.91000247, 835789.5137938 , 834268.47211587, 832908.90075439,
           831704.93922925, 830631.70589502, 829645.90758107, 828661.31255484,
           827639.30254848, 826654.56913448, 825794.76262363, 825088.16975071,
           824540.2293204 , 824169.64935152, 823974.78686248, 823934.72379926,
           824048.31453   , 824333.99897519, 824797.6255898 , 825438.57170204,
           826253.87240821, 827238.63340921, 828387.0265534 , 829691.68151625,
           831143.26846214, 832730.09241152, 834438.12466096, 836251.29475062,
           838152.00329615, 840121.82333294, 842129.43767492, 844087.72200397,
           845868.13256698, 847113.01726278, 847775.27794732, 848172.03989302,
           848480.84121342, 848803.66286792, 849146.23494219, 849511.44139605,
           849902.68300399, 850325.61919522, 850783.84285338, 851235.98349272,
           851626.7925979 , 851818.66766951, 851728.90078917, 851453.79275458,
           851182.63660523, 850984.36652604, 850936.92856183, 851070.69020007,
           851360.88201646, 851729.89790467, 852109.38525478, 852331.39265405,
           852289.25805578, 852052.42968445, 851762.40756602, 851464.4437076 ])
    AEP_rot = np.array([838515.74534606, 840540.47973276, 842665.39775764, 844881.58248964,
           847169.77225837, 849490.39842834, 851731.53530875, 853933.82367575,
           856145.36089982, 858366.7610357 , 860363.54816002, 861802.54772357,
           862841.81223273, 863734.26731151, 864573.00923015, 865292.53551165,
           865860.27234124, 866340.89871928, 866750.48896302, 867098.8312258 ,
           867429.94922581, 867785.81341382, 868195.27497049, 868668.48389518,
           869142.07618188, 869422.37120325, 869168.80503178, 868357.70133152,
           867307.0176534 , 866157.22701589, 864954.26818819, 863690.5614184 ,
           862231.00653949, 860580.85649476, 858888.91262567, 857200.82152759,
           855504.77719346, 853794.62030658, 852036.81118587, 850239.84158866,
           848418.89648681, 846415.62807951, 844201.90052325, 842022.99715444,
           840017.21173995, 838247.5842605 , 836741.46897901, 835510.38164783,
           834561.47881698, 833897.15962989, 833515.64635747, 833417.35889849,
           833573.00744677, 833930.8948609 , 834445.31773577, 835126.11886526,
           835980.18736796, 836985.77564353, 838083.4734594 , 839299.13823637,
           840635.51825252, 842076.21755323, 843407.71734792, 844389.60863703,
           845183.28787712, 845984.95939623, 846877.69589363, 847829.71777412,
           848824.57690897, 849875.36867476, 850958.40758453, 852074.90549183,
           853263.23750307, 854558.29726489, 855981.64566338, 857534.21530306,
           859151.17316378, 860677.24330454, 861854.28844937, 862677.2432988 ,
           863355.21376761, 863965.67671855, 864491.9959609 , 864796.37809657,
           864532.89792238, 863772.74885875, 862829.36342961, 861810.16653338,
           860704.53485582, 859482.18505201, 858051.67446146, 856397.24802014,
           854504.85029536, 852082.88520026, 849120.25066455, 846013.0730422 ,
           842943.85810395, 840008.65066007, 837242.81515848, 834665.72112031,
           832294.73499783, 830144.73496466, 828229.32061359, 826565.54364729,
           825138.54720254, 823908.60135775, 822850.03205043, 822000.38986602,
           821400.28911455, 821056.79040034, 820959.10163121, 821129.64427017,
           821553.13945478, 822186.69875858, 822850.07053322, 823413.29706024,
           824035.55079465, 824859.39638467, 825935.47443579, 827242.72822309,
           828771.0203439 , 830517.44556574, 832453.93236591, 834576.98124301,
           836909.5919445 , 839468.61866111, 842259.12407854, 845267.38936353,
           848435.05780439, 851657.24974102, 854757.88176985, 857695.80360148,
           860547.59223136, 863346.98358641, 866036.59039002, 868332.08047109,
           869702.87676486, 870321.34798081, 870645.51671431, 870858.06101198,
           870968.84814116, 870936.73948829, 870670.2839728 , 870123.15326621,
           869219.25033325, 867631.89731735, 865382.30193105, 862868.12690106,
           860238.94240722, 857586.55151572, 854940.88102373, 852323.19100396,
           849755.78292411, 847261.33747121, 844863.50019132, 842587.81829321,
           840442.21961503, 838421.33659865, 836532.55172129, 834807.47415653,
           833277.61661319, 831962.27759206, 830877.86654401, 830061.34300752,
           829498.05672195, 829081.58168272, 828647.82950267, 828213.54989055,
           827964.16479669, 828023.08350089, 828408.03245143, 829114.39345928,
           830147.09574645, 831512.52219164, 833197.15340071, 835198.66666996,
           837527.63496443, 840185.39588221, 843161.49773622, 846428.27865123,
           849928.80165496, 853580.16334811, 857198.1448885 , 860645.72495751,
           863914.4999426 , 867077.66745917, 870144.79760283, 872800.60627493,
           874578.08348891, 875706.16740794, 876575.2991063 , 877347.94832104,
           878017.83195756, 878513.05872416, 878707.84451544, 878561.54617675,
           878028.71963642, 876849.645879  , 875073.38132547, 873010.00430608])


    run = 4
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


    rotate = 0.
    # dx_start,dy_start,shear_start,rotate_start,turbs_per_row,x_start,y0 = startGrid(nTurbines,
    #             boundaryVertices,boundaryNormals,rotate=rotate,my=1.0)
    #
    # nRows = len(turbs_per_row)
    # prob = Problem()
    # root = prob.root = Group()
    # root.add('dx', IndepVarComp('dx', 0.), promotes=['*'])
    # root.add('dy', IndepVarComp('dy', 0.), promotes=['*'])
    # root.add('shear', IndepVarComp('shear', 0.), promotes=['*'])
    # root.add('rotate', IndepVarComp('rotate', 0.), promotes=['*'])
    #
    # root.add('makeGrid', form_grid(nTurbines=nTurbines,nRows=nRows),promotes=['*'])
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

    nCalls_obj = 0
    nCalls_con = 0

    fig = plt.figure(figsize=[6.5,3.])
    ax1 = plt.subplot2grid((1, 10), (0, 1),colspan=2)
    ax2 = plt.subplot2grid((1, 10), (0, 3),colspan=2)
    ax3 = plt.subplot2grid((1, 10), (0, 5),colspan=2)
    ax4 = plt.subplot2grid((1, 10), (0, 7),colspan=2)

    num = 200
    d_dx = np.linspace(2.*rotorDiameter[0],4.*rotorDiameter[0],num)
    # AEP_dx = np.zeros(num)
    # for i in range(num):
    #         print i
    #         prob['dx'] = d[i]
    #         prob['dy'] = dy_start
    #         prob['shear'] = shear_start
    #         prob['rotate'] = rotate_start
    #
    #         prob.run()
    #         AEP_dx[i] = -prob['negAEP']


    d_dy = np.linspace(2.*rotorDiameter[0],4.*rotorDiameter[0],num)
    # AEP_dy = np.zeros(num)
    # for i in range(num):
    #         print i
    #         prob['dx'] = dx_start
    #         prob['dy'] = d[i]
    #         prob['shear'] = shear_start
    #         prob['rotate'] = rotate_start
    #
    #         prob.run()
    #         AEP_dy[i] = -prob['negAEP']


    d_shear = np.linspace(-2.*rotorDiameter[0],2.*rotorDiameter[0],num)
    # AEP_shear = np.zeros(num)
    # for i in range(num):
    #         print i
    #         prob['dx'] = dx_start
    #         prob['dy'] = dy_start
    #         prob['shear'] = d[i]
    #         prob['rotate'] = rotate_start
    #
    #         prob.run()
    #         AEP_shear[i] = -prob['negAEP']


    d_rot = np.linspace(-30.,30.,num)
    # AEP_rot = np.zeros(num)
    # for i in range(num):
    #         print i
    #         prob['dx'] = dx_start
    #         prob['dy'] = dy_start
    #         prob['shear'] = shear_start
    #         prob['rotate'] = d[i]
    #
    #         prob.run()
    #         AEP_rot[i] = -prob['negAEP']



    ax1.plot(d_dx,AEP_dx,color='C2',alpha=0.5)
    ax1.set_xticks([2.*rotorDiameter[0],3.*rotorDiameter[0],4.*rotorDiameter[0]])
    ax1.set_xticklabels(['2','3','4'])
    ax1.set_xlabel(r'$dx$ '+'(rotor\ndiameters)')

    ax2.plot(d_dy,AEP_dy,color='C2',alpha=0.5)
    ax2.set_xticks([2.*rotorDiameter[0],3.*rotorDiameter[0],4.*rotorDiameter[0]])
    ax2.set_xticklabels(['2','3','4'])
    ax2.set_xlabel(r'$dy$ '+'(rotor\ndiameters)')

    ax3.plot(d_shear,AEP_shear,color='C2',alpha=0.5)
    ax3.set_xticks([-2.*rotorDiameter[0],0.,2.*rotorDiameter[0]])
    ax3.set_xticklabels(['-2','0','2'])
    ax3.set_xlabel(r'$b$ '+'(rotor\ndiameters)')

    ax4.plot(d_rot,AEP_rot,color='C2',alpha=0.5)
    ax4.set_xticks([-30.,0.,30.])
    ax4.set_xticklabels(['-30','0','30'])
    ax4.set_xlabel(r'$\theta$ (degrees)')

    ax1.set_xlim(2.*rotorDiameter[0],4.*rotorDiameter[0])
    ax2.set_xlim(2.*rotorDiameter[0],4.*rotorDiameter[0])
    ax3.set_xlim(-2.*rotorDiameter[0],2.*rotorDiameter[0])
    ax4.set_xlim(-30.,30.)

    ax1.set_yticks([725000.,775000.,825000.,875000.])
    ax2.set_yticks([725000.,775000.,825000.,875000.])
    ax3.set_yticks([725000.,775000.,825000.,875000.])
    ax4.set_yticks([725000.,775000.,825000.,875000.])

    ax1.set_yticklabels(['725','775','825','875'])
    ax2.set_yticklabels(['','','',''])
    ax3.set_yticklabels(['','','',''])
    ax4.set_yticklabels(['','','',''])

    ax1.set_ylim(725000.,875000.)
    ax2.set_ylim(725000.,875000.)
    ax3.set_ylim(725000.,875000.)
    ax4.set_ylim(725000.,875000.)

    ax1.set_ylabel('AEP (GWh)')



    # ax2.set_yticklabels([''])
    # ax3.set_yticklabels([''])
    # ax4.set_yticklabels([''])
    plt.subplots_adjust(top = 0.93, bottom = 0.3, right = 0.96, left = 0.1,
            hspace = 0., wspace = 0.3)

    plt.show()
