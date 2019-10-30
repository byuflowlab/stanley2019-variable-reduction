from ccblade import *
import numpy as np
import scipy.interpolate as sp


r = np.array([0.000,0.010,0.020,0.030,0.040,0.050,0.060,0.070,0.080,0.100,0.120,0.140,0.160,0.180,0.200,0.220,
              0.240,0.260,0.280,0.300,0.320,0.340,0.360,0.380,0.400,0.420,0.440,0.460,0.480,0.500,0.520,0.540,
              0.560,0.580,0.600,0.620,0.640,0.660,0.680,0.700,0.720,0.740,0.760,0.780,0.800,0.820,0.840,0.860,
              0.880,0.900,0.920,0.940,0.960,0.980,1.000])*62.4+2.6

chord = np.array([2600.0,2600.0,2600.0,2620.6,2680.1,2776.7,2910.6,3032.4,3155.4,3400.0,3618.0,3803.9,3958.7,4083.2,4178.5,4245.5,
                  4285.1,4298.4,4286.7,4252.0,4196.1,4121.1,4028.7,3921.1,3800.0,3667.7,3527.8,3384.0,3240.1,3100.0,2966.8,2841.2,
                  2723.1,2612.7,2509.8,2414.5,2326.9,2246.8,2174.3,2109.5,2052.3,2002.8,1960.9,1926.6,1900.0,1879.2,1854.9,1816.1,
                  1751.5,1650.0,1500.5,1291.9,1013.1,652.8 ,200.0 ])/1000.

theta = np.array([20.00,19.87,19.73,19.56,19.39,19.20,18.99,18.76,18.51,17.94,17.07,15.81,14.37,12.95,11.75,10.73,
                  9.71 ,8.71 ,7.74 ,6.81 ,5.93 ,5.12 ,4.39 ,3.75 ,3.21 ,2.76 ,2.38 ,2.05 ,1.77 ,1.53 ,1.32 ,1.13 ,
                  0.95 ,0.78 ,0.60 ,0.42 ,0.27 ,0.13 ,0.01 ,-0.11,-0.22,-0.34,-0.47,-0.60,-0.75,-0.91,-1.07,-1.24,
                  -1.45,-1.70,-2.05,-2.54,-3.14,-3.84 ,-4.62 ])


# r = np.array([0.0, 0.02, 0.1371, 0.2118, 0.3188, 0.5236, 0.6781, 0.8967, 1.0])
#
# chord = np.zeros_like(r)
# theta = np.zeros_like(r)
#
# for i in range(len(r)):
#     index = 0
#     on = True
#     while on:
#         if r[i] >= r1[index] and r[i] <= r1[index+1]:
#             on = False
#         else:
#             index += 1
#     chord[i] = (r[i]-r1[index])/(r1[index+1]-r1[index])*(chord1[index+1]-chord1[index])+chord1[index]
#     theta[i] = (r[i]-r1[index])/(r1[index+1]-r1[index])*(theta1[index+1]-theta1[index])+theta1[index]

# r = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
#               28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
#               56.1667, 58.9000, 61.6333])
# chord = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
#                   3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
# theta = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
#                   6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])

import os
afinit = CCAirfoil.initFromAerodynFile  # just for shorthand
basepath = '335MW_AFFiles' + os.path.sep

# load all airfoils
airfoil_types = [0]*7
airfoil_types[0] = afinit(basepath + 'Cylinder1.dat')
airfoil_types[1] = afinit(basepath + 'FX77-W-500.dat')
airfoil_types[2] = afinit(basepath + 'FX77-W-400.dat')
airfoil_types[3] = afinit(basepath + 'DU00-W2-350.dat')
airfoil_types[4] = afinit(basepath + 'DU97-W-300.dat')
airfoil_types[5] = afinit(basepath + 'DU91-W2-250.dat')
airfoil_types[6] = afinit(basepath + 'DU08-W-210.dat')

# basepath = '5MW_AFFiles' + os.path.sep

# airfoil_types = [0]*8
# airfoil_types[0] = afinit(basepath + 'Cylinder1.dat')
# airfoil_types[1] = afinit(basepath + 'Cylinder2.dat')
# airfoil_types[2] = afinit(basepath + 'DU40_A17.dat')
# airfoil_types[3] = afinit(basepath + 'DU35_A17.dat')
# airfoil_types[4] = afinit(basepath + 'DU30_A17.dat')
# airfoil_types[5] = afinit(basepath + 'DU25_A17.dat')
# airfoil_types[6] = afinit(basepath + 'DU21_A17.dat')
# airfoil_types[7] = afinit(basepath + 'NACA64_A17.dat')


# r = np.array([0.0, 0.02, 0.1371, 0.2118, 0.3188, 0.5236, 0.6781, 0.8967, 1.0])
# af_idx = [0, 0, 1, 2, 3, 4, 5, 6, 6]
af_idx = [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4,
            4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 6, 6, 6, 6, 6, 6]

# af_idx = [0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6]

# af_idx = [0, 0, 0, 2, 3, 4, 5, 6, 6]

af = [0]*len(r)
for i in range(len(r)):
    af[i] = airfoil_types[af_idx[i]]

tilt = 5.0
precone = 2.5
hubHt = 110.
Rhub = 2.6
Rtip = 65.

rotor = CCBlade(r, chord, theta, af, Rhub, Rtip, precone, tilt, hubHt)

omega_rpm = np.array([0.,3.9114391143911433,11.75,11.75])
omega_ws = np.array([0.,3.060344827586208,9.331896551724139, 25.1,])
f_omega = sp.interp1d(omega_ws,omega_rpm)



pitch_angle = np.array([1.09,1.09,1.09,1.09,0.8645195667798404,
            3.102803043445931,4.465909935810053,5.732360157624726,6.706981585665069,
            7.486740687501545,8.461981709583881,9.339636670053778,10.510049815360972,
            11.583186696076726,12.949081761629778,14.705630870653547,16.95042008476046,
            18.707898584847207,20.17261890009665,21.83313093261295,23.493642965129244,
            25.447532776524817,26.424322783712107])
pitch_ws = np.array([0.,3.1528662420382156,7.67515923566879,  9.26751592356688,  9.777070063694268,
            9.968152866242038, 10.222929936305732,10.605095541401274,10.987261146496815,
            11.305732484076435,11.815286624203821,12.261146496815286,12.898089171974522,
            13.535031847133759,14.363057324840764,15.509554140127388,17.038216560509554,
            18.37579617834395, 19.522292993630575,20.923566878980893,22.32484076433121,
            24.044585987261147,25.])
f_pitch = sp.interp1d(pitch_ws,pitch_angle)



Uinf = np.linspace(0.1,25.,50)

import matplotlib.pyplot as plt
plt.plot(Uinf,f_omega(Uinf),'o')
plt.show()

print 'Uinf: ', repr(Uinf)
# print 'f_omega(Uinf): ', f_omega(Uinf)
# print 'f_pitch(Uinf): ', f_pitch(Uinf)

P, T, Q, M, CP, CT, CQ, CM = rotor.evaluate(Uinf, f_omega(Uinf), f_pitch(Uinf), coefficients=True)
# print 'P: ', P
# print 'T: ', T
# print 'Q: ', Q
# print 'M: ', M
# print 'CP: ', CP
# print 'CQ: ', CQ

print 'Uinf: ', repr(Uinf)
print 'CT: ', repr(CT)
# print 'CM: ', CM

import matplotlib.pyplot as plt
plt.plot(Uinf,CT,'o')
plt.show()
