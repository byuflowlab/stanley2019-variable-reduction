import numpy as np
from akima import Akima, akima_interp
from math import radians
import sys
import scipy.special as spec
sys.dont_write_bytecode = True


def novatoRose(nDirections,nSpeeds=1):
    data = np.array([[0.839,     0.710,     0.562,     0.591,     0.041,     0.003],
                     [0.882,     0.656,     0.426,     0.426,     0.027,     0.002],
                     [0.738,     0.427,     0.246,     0.243,     0.027,     0.010],
                     [0.563,     0.279,     0.130,     0.125,     0.020,     0.010],
                     [0.491,     0.180,     0.074,     0.063,     0.008,     0.005],
                     [0.445,     0.154,     0.065,     0.041,     0.001,     0.001],
                     [0.420,     0.147,     0.041,     0.028,     0.002,     0.001],
                     [0.443,     0.148,     0.044,     0.024,     0.000,     0.000],
                     [0.448,     0.177,     0.041,     0.022,     0.001,     0.000],
                     [0.487,     0.201,     0.054,     0.017,     0.002,     0.001],
                     [0.497,     0.269,     0.092,     0.032,     0.003,     0.000],
                     [0.510,     0.344,     0.147,     0.045,     0.001,     0.000],
                     [0.676,     0.562,     0.231,     0.052,     0.000,     0.000],
                     [1.021,     0.894,     0.293,     0.064,     0.001,     0.000],
                     [1.114,     0.654,     0.173,     0.069,     0.006,     0.001],
                     [0.897,     0.483,     0.240,     0.133,     0.012,     0.002],
                     [0.699,     0.419,     0.225,     0.140,     0.015,     0.003],
                     [0.681,     0.422,     0.218,     0.192,     0.033,     0.014],
                     [0.758,     0.551,     0.262,     0.189,     0.058,     0.027],
                     [0.855,     0.728,     0.335,     0.218,     0.048,     0.040],
                     [0.938,     0.667,     0.335,     0.264,     0.044,     0.022],
                     [0.757,     0.554,     0.356,     0.367,     0.053,     0.020],
                     [0.584,     0.444,     0.409,     0.655,     0.143,     0.014],
                     [0.471,     0.429,     0.472,     1.026,     0.350,     0.023],
                     [0.406,     0.383,     0.427,     1.088,     0.297,     0.030],
                     [0.334,     0.340,     0.376,     0.878,     0.243,     0.068],
                     [0.263,     0.264,     0.315,     0.517,     0.172,     0.064],
                     [0.215,     0.230,     0.260,     0.515,     0.229,     0.124],
                     [0.216,     0.207,     0.285,     0.593,     0.271,     0.124],
                     [0.221,     0.196,     0.274,     0.525,     0.162,     0.035],
                     [0.261,     0.214,     0.271,     0.437,     0.116,     0.018],
                     [0.331,     0.276,     0.321,     0.542,     0.116,     0.020],
                     [0.503,     0.359,     0.419,     0.870,     0.203,     0.016],
                     [0.763,     0.542,     0.478,     0.875,     0.122,     0.010],
                     [1.097,     0.923,     0.742,     1.059,     0.105,     0.005],
                     [1.132,     1.016,     0.841,     1.061,     0.083,     0.003]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)

    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)
    return dirs, frequencies, speeds


def ukiahRose(nDirections,nSpeeds=1):
    data = np.array([[0.479,     0.639,     0.516,     0.471,     0.025,     0.002],
                     [0.374,     0.524,     0.426,     0.339,     0.015,     0.001],
                     [0.337,     0.489,     0.356,     0.201,     0.009,     0.000],
                     [0.243,     0.353,     0.184,     0.068,     0.002,     0.002],
                     [0.156,     0.166,     0.065,     0.024,     0.000,     0.000],
                     [0.154,     0.115,     0.034,     0.006,     0.000,     0.000],
                     [0.180,     0.129,     0.025,     0.005,     0.000,     0.000],
                     [0.166,     0.093,     0.021,     0.004,     0.000,     0.000],
                     [0.188,     0.095,     0.018,     0.005,     0.000,     0.000],
                     [0.219,     0.149,     0.034,     0.015,     0.002,     0.000],
                     [0.197,     0.142,     0.053,     0.035,     0.003,     0.002],
                     [0.255,     0.195,     0.115,     0.081,     0.004,     0.001],
                     [0.369,     0.457,     0.340,     0.292,     0.022,     0.003],
                     [0.374,     0.443,     0.408,     0.408,     0.033,     0.005],
                     [0.441,     0.584,     0.537,     0.580,     0.037,     0.005],
                     [0.529,     0.809,     0.763,     0.857,     0.073,     0.014],
                     [0.642,     0.773,     0.641,     0.642,     0.066,     0.006],
                     [0.676,     0.788,     0.564,     0.499,     0.055,     0.008],
                     [0.481,     0.549,     0.360,     0.313,     0.035,     0.009],
                     [0.295,     0.241,     0.169,     0.136,     0.013,     0.003],
                     [0.215,     0.134,     0.093,     0.064,     0.007,     0.000],
                     [0.140,     0.091,     0.067,     0.076,     0.005,     0.001],
                     [0.103,     0.059,     0.060,     0.052,     0.002,     0.001],
                     [0.094,     0.058,     0.077,     0.071,     0.005,     0.002],
                     [0.130,     0.120,     0.122,     0.158,     0.012,     0.001],
                     [0.146,     0.110,     0.119,     0.153,     0.012,     0.002],
                     [0.252,     0.148,     0.160,     0.202,     0.022,     0.003],
                     [0.374,     0.267,     0.261,     0.408,     0.059,     0.009],
                     [0.567,     0.240,     0.245,     0.453,     0.078,     0.012],
                     [0.587,     0.288,     0.314,     0.603,     0.104,     0.015],
                     [0.560,     0.381,     0.461,     0.890,     0.155,     0.022],
                     [0.472,     0.317,     0.389,     0.719,     0.108,     0.014],
                     [0.478,     0.349,     0.390,     0.700,     0.094,     0.009],
                     [0.548,     0.537,     0.516,     0.809,     0.099,     0.010],
                     [0.531,     0.451,     0.387,     0.533,     0.055,     0.005],
                     [0.486,     0.549,     0.494,     0.553,     0.047,     0.003]])
    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)

    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)
    return dirs, frequencies, speeds


def victorvilleRose(nDirections,nSpeeds=1):
    data = np.array([[0.484,     0.289,     0.176,     0.136,     0.017,     0.015],
                     [0.551,     0.352,     0.218,     0.204,     0.043,     0.017],
                     [0.457,     0.334,     0.259,     0.301,     0.053,     0.024],
                     [0.378,     0.271,     0.237,     0.345,     0.072,     0.019],
                     [0.290,     0.196,     0.170,     0.296,     0.070,     0.017],
                     [0.235,     0.156,     0.120,     0.225,     0.055,     0.009],
                     [0.192,     0.108,     0.075,     0.125,     0.029,     0.005],
                     [0.140,     0.073,     0.049,     0.066,     0.015,     0.004],
                     [0.152,     0.069,     0.035,     0.037,     0.010,     0.001],
                     [0.140,     0.067,     0.028,     0.032,     0.004,     0.002],
                     [0.151,     0.067,     0.031,     0.027,     0.006,     0.003],
                     [0.183,     0.098,     0.050,     0.045,     0.011,     0.003],
                     [0.341,     0.220,     0.140,     0.091,     0.019,     0.004],
                     [0.517,     0.364,     0.246,     0.122,     0.024,     0.011],
                     [0.855,     0.620,     0.362,     0.218,     0.053,     0.046],
                     [1.307,     0.996,     0.588,     0.371,     0.144,     0.129],
                     [1.915,     1.322,     0.699,     0.507,     0.260,     0.195],
                     [2.154,     1.532,     0.803,     0.766,     0.488,     0.462],
                     [2.014,     1.122,     0.617,     1.074,     0.714,     0.602],
                     [1.792,     1.020,     0.614,     1.191,     0.666,     0.427],
                     [1.475,     0.838,     0.599,     1.072,     0.395,     0.139],
                     [1.373,     0.802,     0.528,     0.668,     0.137,     0.039],
                     [1.160,     0.642,     0.426,     0.418,     0.077,     0.023],
                     [1.045,     0.644,     0.405,     0.424,     0.113,     0.044],
                     [0.917,     0.598,     0.477,     0.574,     0.180,     0.089],
                     [0.775,     0.536,     0.519,     0.843,     0.311,     0.157],
                     [0.692,     0.532,     0.571,     1.277,     0.610,     0.365],
                     [0.519,     0.462,     0.543,     1.345,     0.766,     0.543],
                     [0.414,     0.344,     0.360,     0.981,     0.560,     0.363],
                     [0.344,     0.247,     0.264,     0.572,     0.293,     0.197],
                     [0.281,     0.179,     0.139,     0.227,     0.094,     0.075],
                     [0.240,     0.123,     0.073,     0.095,     0.032,     0.018],
                     [0.273,     0.110,     0.057,     0.051,     0.012,     0.006],
                     [0.283,     0.136,     0.057,     0.038,     0.007,     0.001],
                     [0.329,     0.153,     0.071,     0.048,     0.007,     0.001],
                     [0.431,     0.218,     0.108,     0.074,     0.009,     0.006]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)

    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)
    return dirs, frequencies, speeds


#newhall is one dominant wind direction (North) with the semi-dominant direction exactly opposite (South)
def newhallRose(nDirections,nSpeeds=1):
    data = np.array([[8.192,     2.952,     1.761,     1.511,     0.367,     0.108],
                     [0.629,     0.312,     0.217,     0.244,     0.059,     0.025],
                     [0.805,     0.377,     0.166,     0.092,     0.010,     0.003],
                     [0.537,     0.245,     0.149,     0.080,     0.004,     0.001],
                     [0.762,     0.231,     0.180,     0.162,     0.019,     0.005],
                     [0.621,     0.149,     0.116,     0.207,     0.049,     0.008],
                     [0.643,     0.118,     0.104,     0.252,     0.067,     0.014],
                     [0.547,     0.062,     0.041,     0.082,     0.012,     0.005],
                     [0.693,     0.067,     0.023,     0.030,     0.001,     0.000],
                     [0.539,     0.057,     0.027,     0.007,     0.000,     0.000],
                     [0.867,     0.105,     0.032,     0.015,     0.000,     0.000],
                     [0.703,     0.052,     0.019,     0.005,     0.000,     0.000],
                     [0.918,     0.055,     0.023,     0.005,     0.000,     0.000],
                     [0.705,     0.046,     0.008,     0.017,     0.000,     0.000],
                     [0.925,     0.067,     0.017,     0.003,     0.000,     0.000],
                     [0.821,     0.071,     0.018,     0.002,     0.000,     0.000],
                     [1.093,     0.099,     0.036,     0.008,     0.003,     0.001],
                     [0.824,     0.112,     0.084,     0.088,     0.013,     0.004],
                     [8.345,     2.182,     1.119,     0.746,     0.220,     0.070],
                     [1.321,     0.261,     0.116,     0.080,     0.014,     0.010],
                     [1.674,     0.550,     0.250,     0.157,     0.023,     0.008],
                     [0.643,     0.196,     0.083,     0.022,     0.000,     0.000],
                     [0.806,     0.221,     0.079,     0.009,     0.000,     0.000],
                     [0.613,     0.107,     0.028,     0.008,     0.001,     0.000],
                     [0.908,     0.211,     0.094,     0.074,     0.006,     0.002],
                     [1.098,     0.463,     0.098,     0.023,     0.000,     0.000],
                     [1.337,     0.869,     0.357,     0.096,     0.003,     0.000],
                     [1.035,     0.642,     0.242,     0.024,     0.001,     0.000],
                     [1.578,     1.097,     0.560,     0.068,     0.000,     0.000],
                     [1.229,     0.594,     0.230,     0.019,     0.001,     0.000],
                     [1.047,     0.345,     0.150,     0.044,     0.013,     0.001],
                     [0.830,     0.206,     0.117,     0.057,     0.007,     0.002],
                     [1.252,     0.474,     0.210,     0.153,     0.029,     0.005],
                     [0.724,     0.375,     0.164,     0.112,     0.023,     0.010],
                     [1.247,     0.567,     0.351,     0.363,     0.087,     0.024],
                     [0.723,     0.326,     0.231,     0.315,     0.122,     0.041]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)

    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)
    return dirs, frequencies, speeds


#beale is one dominant wind directions with the semi-dominant direction exactly opposite
#a little less extreme than newhall though
def bealeRose(nDirections,nSpeeds=1):
    data = np.array([[0.856,     0.347,     0.149,     0.099,     0.027,     0.009],
                     [0.611,     0.205,     0.078,     0.063,     0.016,     0.007],
                     [0.790,     0.182,     0.074,     0.060,     0.020,     0.007],
                     [0.487,     0.127,     0.046,     0.040,     0.009,     0.006],
                     [0.492,     0.116,     0.039,     0.028,     0.007,     0.002],
                     [0.557,     0.125,     0.035,     0.028,     0.003,     0.002],
                     [0.598,     0.134,     0.032,     0.020,     0.003,     0.000],
                     [0.657,     0.167,     0.045,     0.015,     0.001,     0.000],
                     [0.795,     0.261,     0.081,     0.033,     0.001,     0.000],
                     [0.929,     0.437,     0.177,     0.064,     0.003,     0.002],
                     [0.882,     0.571,     0.291,     0.107,     0.007,     0.003],
                     [0.886,     0.519,     0.324,     0.218,     0.030,     0.008],
                     [0.843,     0.634,     0.487,     0.431,     0.094,     0.031],
                     [0.937,     0.801,     0.635,     0.738,     0.270,     0.110],
                     [1.227,     1.002,     0.888,     1.435,     0.514,     0.286],
                     [1.044,     0.928,     0.906,     1.893,     0.742,     0.482],
                     [1.189,     1.017,     1.067,     1.992,     0.767,     0.603],
                     [1.483,     1.131,     1.015,     1.617,     0.513,     0.306],
                     [1.223,     0.952,     0.853,     1.169,     0.269,     0.113],
                     [0.989,     0.763,     0.709,     0.825,     0.130,     0.042],
                     [1.023,     0.698,     0.555,     0.497,     0.058,     0.022],
                     [0.736,     0.480,     0.299,     0.248,     0.032,     0.017],
                     [0.601,     0.331,     0.181,     0.108,     0.018,     0.014],
                     [0.654,     0.275,     0.118,     0.060,     0.016,     0.006],
                     [0.467,     0.193,     0.073,     0.036,     0.009,     0.005],
                     [0.454,     0.169,     0.055,     0.025,     0.008,     0.002],
                     [0.556,     0.190,     0.051,     0.021,     0.005,     0.002],
                     [0.526,     0.224,     0.069,     0.018,     0.003,     0.001],
                     [0.508,     0.219,     0.079,     0.032,     0.005,     0.001],
                     [0.714,     0.305,     0.132,     0.077,     0.016,     0.005],
                     [0.785,     0.398,     0.214,     0.184,     0.048,     0.020],
                     [0.780,     0.482,     0.303,     0.316,     0.112,     0.055],
                     [1.107,     0.696,     0.480,     0.586,     0.255,     0.139],
                     [1.032,     0.692,     0.512,     0.674,     0.289,     0.132],
                     [1.093,     0.686,     0.503,     0.567,     0.175,     0.066],
                     [1.153,     0.591,     0.310,     0.266,     0.066,     0.024]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)

    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)
    return dirs, frequencies, speeds


def northIslandRose(nDirections,nSpeeds=1):
    data = np.array([[1.194,     0.560,     0.216,     0.090,     0.009,     0.003],
                     [0.727,     0.303,     0.090,     0.032,     0.003,     0.001],
                     [0.554,     0.200,     0.046,     0.022,     0.004,     0.003],
                     [0.460,     0.129,     0.035,     0.016,     0.004,     0.002],
                     [0.243,     0.069,     0.022,     0.014,     0.005,     0.003],
                     [0.223,     0.055,     0.014,     0.013,     0.006,     0.003],
                     [0.207,     0.043,     0.012,     0.011,     0.003,     0.002],
                     [0.271,     0.064,     0.020,     0.015,     0.004,     0.001],
                     [0.386,     0.118,     0.026,     0.014,     0.001,     0.001],
                     [0.690,     0.232,     0.038,     0.018,     0.001,     0.002],
                     [0.708,     0.223,     0.046,     0.020,     0.005,     0.001],
                     [0.633,     0.209,     0.064,     0.033,     0.003,     0.001],
                     [0.678,     0.228,     0.109,     0.062,     0.011,     0.004],
                     [0.496,     0.208,     0.125,     0.083,     0.020,     0.004],
                     [0.524,     0.258,     0.143,     0.110,     0.024,     0.010],
                     [0.611,     0.314,     0.157,     0.134,     0.041,     0.023],
                     [0.676,     0.397,     0.228,     0.217,     0.090,     0.072],
                     [0.736,     0.508,     0.320,     0.323,     0.110,     0.089],
                     [1.115,     0.730,     0.467,     0.435,     0.101,     0.049],
                     [0.891,     0.709,     0.500,     0.401,     0.050,     0.025],
                     [0.930,     0.829,     0.707,     0.448,     0.043,     0.020],
                     [1.023,     0.880,     0.947,     0.716,     0.043,     0.014],
                     [0.675,     0.468,     0.395,     0.335,     0.030,     0.011],
                     [0.646,     0.431,     0.369,     0.278,     0.036,     0.012],
                     [0.661,     0.485,     0.476,     0.352,     0.029,     0.008],
                     [0.685,     0.574,     0.552,     0.477,     0.035,     0.011],
                     [0.723,     0.693,     0.763,     0.681,     0.054,     0.019],
                     [1.024,     0.994,     1.082,     1.099,     0.117,     0.044],
                     [0.866,     0.903,     1.089,     1.363,     0.172,     0.069],
                     [0.848,     0.887,     1.091,     1.720,     0.227,     0.082],
                     [1.206,     1.140,     1.219,     2.052,     0.287,     0.070],
                     [1.181,     1.204,     1.209,     1.986,     0.285,     0.039],
                     [1.227,     1.195,     1.159,     1.743,     0.318,     0.028],
                     [1.243,     1.007,     0.807,     1.062,     0.164,     0.017],
                     [1.239,     0.906,     0.567,     0.491,     0.058,     0.008],
                     [1.099,     0.678,     0.330,     0.193,     0.019,     0.004]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    wd = np.linspace(0.,360.-360/float(m), m)
    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)

    return dirs, frequencies, speeds


def alturasRose(nDirections,nSpeeds=1):
    data = np.array([[0.447,     0.281,     0.233,     0.292,     0.054,     0.015],
                     [0.462,     0.338,     0.293,     0.447,     0.097,     0.017],
                     [0.556,     0.388,     0.264,     0.295,     0.054,     0.011],
                     [0.609,     0.416,     0.207,     0.167,     0.029,     0.008],
                     [0.616,     0.404,     0.159,     0.114,     0.010,     0.003],
                     [0.665,     0.415,     0.165,     0.113,     0.009,     0.002],
                     [0.631,     0.359,     0.141,     0.093,     0.010,     0.001],
                     [0.600,     0.319,     0.116,     0.060,     0.007,     0.001],
                     [0.584,     0.246,     0.067,     0.028,     0.006,     0.000],
                     [0.547,     0.187,     0.052,     0.031,     0.002,     0.001],
                     [0.535,     0.163,     0.044,     0.016,     0.003,     0.000],
                     [0.509,     0.171,     0.048,     0.013,     0.001,     0.001],
                     [0.559,     0.205,     0.058,     0.014,     0.001,     0.000],
                     [0.559,     0.230,     0.067,     0.027,     0.001,     0.000],
                     [0.572,     0.287,     0.113,     0.067,     0.009,     0.003],
                     [0.566,     0.307,     0.186,     0.192,     0.038,     0.008],
                     [0.568,     0.362,     0.312,     0.510,     0.117,     0.030],
                     [0.546,     0.426,     0.442,     0.922,     0.345,     0.155],
                     [0.569,     0.449,     0.438,     0.937,     0.352,     0.221],
                     [0.528,     0.367,     0.309,     0.547,     0.205,     0.147],
                     [0.496,     0.311,     0.219,     0.283,     0.112,     0.096],
                     [0.470,     0.306,     0.189,     0.293,     0.118,     0.079],
                     [0.443,     0.295,     0.214,     0.347,     0.138,     0.073],
                     [0.494,     0.367,     0.278,     0.472,     0.205,     0.088],
                     [0.559,     0.489,     0.389,     0.675,     0.292,     0.123],
                     [0.660,     0.653,     0.662,     1.076,     0.438,     0.173],
                     [0.669,     0.819,     0.813,     1.501,     0.625,     0.186],
                     [0.645,     0.736,     0.752,     1.318,     0.530,     0.129],
                     [0.598,     0.550,     0.485,     0.827,     0.295,     0.071],
                     [0.518,     0.426,     0.384,     0.564,     0.136,     0.033],
                     [0.452,     0.345,     0.346,     0.551,     0.144,     0.039],
                     [0.426,     0.289,     0.303,     0.567,     0.196,     0.080],
                     [0.433,     0.290,     0.291,     0.541,     0.175,     0.075],
                     [0.432,     0.252,     0.275,     0.396,     0.102,     0.034],
                     [0.406,     0.248,     0.211,     0.246,     0.040,     0.009],
                     [0.391,     0.237,     0.161,     0.203,     0.029,     0.007]])

    m,n = np.shape(data)
    wf = np.zeros(m)
    for i in range(m):
        wf[i] = sum(data[i])/100.

    speeds = np.array([2.+(4.9-2.)/2.,5.+(6.9-5.)/2.,7.+(9.9-7.)/2.,10.+(14.9-10.)/2.,15.+(19.9-15.)/2.,22.])
    ws = np.zeros(m)
    for i in range(m):
        for j in range(n):
            ws[i] += data[i][j]/100.*speeds[j]
        ws[i] = ws[i]/wf[i]

    import matplotlib.pyplot as plt

    # wf = np.array([78.,99.,93.,87.,78.,82.,73.,66.,52.,49.,44.,44.,49.,53.,62.,78.,
    #         116.,173.,181.,128.,91.,90.,92.,40.,153.,220.,280.,247.,170.,124.,112.,112.,
    #         109.,90.,70.,61.])
    # wf = wf/sum(wf)
    # ws = np.ones(nDirections)*8.

    wd = np.linspace(0.,360.-360/float(m), m)
    dirs, frequencies, speeds = process_rose_data(ws,wd,wf,nDirections,nSpeeds)

    return dirs, frequencies, speeds


def amaliaRose(nDirections,nSpeeds=1):
    windSpeeds = np.array([6.53163342, 6.11908394, 6.13415514, 6.0614625,  6.21344602,
                     5.87000793, 5.62161519, 5.96779107, 6.33589422, 6.4668016,
                     7.9854581,  7.6894432,  7.5089221,  7.48638098, 7.65764618,
                     6.82414044, 6.36728201, 5.95982999, 6.05942132, 6.1176321,
                     5.50987893, 4.18461796, 4.82863115, 0.,         0.,         0.,
                     5.94115843, 5.94914252, 5.59386528, 6.42332524, 7.67904937,
                     7.89618066, 8.84560463, 8.51601497, 8.40826823, 7.89479475,
                     7.86194762, 7.9242645,  8.56269962, 8.94563889, 9.82636368,
                     10.11153102, 9.71402212, 9.95233636,  10.35446959, 9.67156182,
                     9.62462527, 8.83545158, 8.18011771, 7.9372492,  7.68726143,
                     7.88134508, 7.31394723, 7.01839896, 6.82858346, 7.06213432,
                     7.01949894, 7.00575122, 7.78735165, 7.52836352, 7.21392201,
                     7.4356621,  7.54099962, 7.61335262, 7.90293531, 7.16021596,
                     7.19617087, 7.5593657,  7.03278586, 6.76105501, 6.48004694,
                     6.94716392])

    windFrequencies = np.array([1.17812570e-02, 1.09958570e-02, 9.60626600e-03, 1.21236860e-02,
                                1.04722450e-02, 1.00695140e-02, 9.68687400e-03, 1.00090550e-02,
                                1.03715390e-02, 1.12172280e-02, 1.52249700e-02, 1.56279300e-02,
                                1.57488780e-02, 1.70577560e-02, 1.93535770e-02, 1.41980570e-02,
                                1.20632100e-02, 1.20229000e-02, 1.32111160e-02, 1.74605400e-02,
                                1.72994400e-02, 1.43993790e-02, 7.87436000e-03, 0.00000000e+00,
                                2.01390000e-05, 0.00000000e+00, 3.42360000e-04, 3.56458900e-03,
                                7.18957000e-03, 8.80068000e-03, 1.13583200e-02, 1.41576700e-02,
                                1.66951900e-02, 1.63125500e-02, 1.31709000e-02, 1.09153300e-02,
                                9.48553000e-03, 1.01097900e-02, 1.18819700e-02, 1.26069900e-02,
                                1.58895900e-02, 1.77021600e-02, 2.04208100e-02, 2.27972500e-02,
                                2.95438600e-02, 3.02891700e-02, 2.69861000e-02, 2.21527500e-02,
                                2.12465500e-02, 1.82861400e-02, 1.66147400e-02, 1.90111800e-02,
                                1.90514500e-02, 1.63932050e-02, 1.76215200e-02, 1.65341460e-02,
                                1.44597600e-02, 1.40370300e-02, 1.65745000e-02, 1.56278200e-02,
                                1.53459200e-02, 1.75210100e-02, 1.59702700e-02, 1.51041500e-02,
                                1.45201100e-02, 1.34527800e-02, 1.47819600e-02, 1.33923300e-02,
                                1.10562900e-02, 1.04521380e-02, 1.16201970e-02, 1.10562700e-02])

    windDirections = np.linspace(0.,360.-360./float(len(windSpeeds)), len(windSpeeds))
    dirs, frequencies, speeds = process_rose_data(windSpeeds,windDirections,windFrequencies,nDirections,nSpeeds)

    return dirs, frequencies, speeds


def setup_weibull(windDirections,windFrequencies,windSpeeds,nSpeeds):

    def Weibull(x,L):
        k = 2.0
        if L < 0.0001:
            L = 0.0001
        return (k/L)*(x/L)**(k-1)*np.exp(-(x/L)**k)

    if nSpeeds == 1:
        return windDirections, windFrequencies, windSpeeds

    nDirections = len(windDirections)
    dirs = np.zeros(nDirections*nSpeeds)
    freqs = np.zeros(nDirections*nSpeeds)
    speeds = np.zeros(nDirections*nSpeeds)
    print np.shape(speeds)
    #direction loops
    for i in range(nDirections):
        for j in range(nSpeeds):
            dirs[i*nSpeeds+j] = windDirections[i]

    #speed and frequency loops
    for i in range(nDirections):
        avg_speed = windSpeeds[i]
        speed_dist = np.linspace((25.)/(2.*float(nSpeeds))+0.001,25.-(25.)/(2.*float(nSpeeds)),nSpeeds)
        # print 'speed_dist: ', speed_dist
        dspeed = speed_dist[1]-speed_dist[0]
        num_int = 1000
        for j in range(nSpeeds):
            speed_int = np.linspace(speed_dist[j]-dspeed/2.,speed_dist[j]+dspeed/2.,num_int)

            k = 2.0
            scale = avg_speed/(spec.gamma(1.0+1./k))

            freq_int = Weibull(speed_int,scale)
            speed_freq = np.trapz(freq_int,speed_int)
            speeds[i*nSpeeds+j] = speed_dist[j]
            freqs[i*nSpeeds+j] = speed_freq*windFrequencies[i]

    return dirs, freqs, speeds


def process_rose_data(windSpeeds,windDirections,windFrequencies,nDirections,nSpeeds):
    """take arbitrary wind rose data, and return data equally spaced for a desired
    number of directions"""
    spline_freq = Akima(windDirections, windFrequencies)
    spline_speed = Akima(windDirections, windSpeeds)
    num = nDirections
    dirs = np.linspace(0.,360.-360./float(num), num)
    ddir = dirs[1]-dirs[0]

    frequencies = np.zeros(num)
    speeds = np.zeros(num)

    num_int = 100

    dir_int1 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int1 = np.zeros(num_int/2)
    speed_freq_int1 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int1[j],_,_,_ = spline_freq.interp(dir_int1[j])
        ws,_,_,_ = spline_speed.interp(dir_int1[j])
        speed_freq_int1[j] = freq_int1[j]*ws

    dir_int2 = np.linspace(dirs[0],dirs[0]+ddir/2.,num_int/2)
    freq_int2 = np.zeros(num_int/2)
    speed_freq_int2 = np.zeros(num_int/2)
    for j in range(num_int/2):
        freq_int2[j],_,_,_ = spline_freq.interp(dir_int2[j])
        ws,_,_,_ = spline_speed.interp(dir_int2[j])
        speed_freq_int2[j] = freq_int2[j]*ws

    frequencies[0] = np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2)
    speeds[0] = (np.trapz(speed_freq_int1,dir_int1)+np.trapz(speed_freq_int2,dir_int2))/\
        (np.trapz(freq_int1,dir_int1)+np.trapz(freq_int2,dir_int2))

    for i in range(1,num):
        dir_int = np.linspace(dirs[i]-ddir/2.,dirs[i]+ddir/2.,num_int)
        freq_int = np.zeros(num_int)
        speed_freq_int = np.zeros(num_int)
        for j in range(num_int):
            freq_int[j],_,_,_ = spline_freq.interp(dir_int[j])
            ws,_,_,_ = spline_speed.interp(dir_int[j])
            speed_freq_int[j] = freq_int[j]*ws
        frequencies[i] = np.trapz(freq_int,dir_int)
        speeds[i] = np.trapz(speed_freq_int,dir_int)/np.trapz(freq_int,dir_int)

    frequencies = frequencies/sum(frequencies)
    for i in range(len(frequencies)):
        if speeds[i] < 0.:
            speeds[i] = 0.
        if frequencies[i] < 0.:
            frequencies[i] = 0.

    dirs, frequencies, speeds = setup_weibull(dirs,frequencies,speeds,nSpeeds)

    return dirs, frequencies, speeds


if __name__=="__main__":
    import matplotlib.pyplot as plt
    nDirections = 30

    wd, wf, ws = ukiahRose(nDirections)
    print 'directions: ', wd
    print 'frequencies: ', wf
    print sum(wf)
    print 'speeds: ', ws

    wd, wf, ws = ukiahRose(nDirections,nSpeeds=10)
    print 'directions: ', wd
    print 'frequencies: ', wf
    print sum(wf)
    print 'speeds: ', ws



    bottom = 0
    width = (2*np.pi) / nDirections

    wd, wf, ws = ukiahRose(nDirections)
    wd -= wd[np.argmax(wf)]
    wd += 270.
    wd +=180./float(nDirections)
    for i in range(nDirections):
        wd[i] = radians(wd[i])*-1.
    """wind rose"""
    plt.figure(1)
    max_height = max(wf)
    ax = plt.subplot(121, polar=True)
    bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    plt.title('Ukiah Rose', y=1.09,fontsize=18,family='serif')
    plt.tight_layout()
    """wind speeds"""
    max_height = max(ws)
    ax = plt.subplot(122, polar=True)
    bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    plt.title('Ukiah Speeds', y=1.09,fontsize=18,family='serif')
    plt.tight_layout()
    # plt.savefig('ukiahRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = victorvilleRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(2)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Victorville Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Victorville Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('victorvilleRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = newhallRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(3)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Newhall Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Newhall Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('newhallRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = bealeRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(4)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Beale Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Beale Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('bealeRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = northIslandRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(5)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('North Island Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('North Island Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('northIslandRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = alturasRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(6)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Alturas Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Alturas Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('alturasRose_r.pdf',transparent=True)
    #
    #
    # wd, wf, ws = amaliaRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(7)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Amalia Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Amalia Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('amaliaRose_r.pdf',transparent=True)
    #
    # wd, wf, ws = novatoRose(nDirections)
    # wd -= wd[np.argmax(wf)]
    # wd += 270.
    # wd +=180./float(nDirections)
    # for i in range(nDirections):
    #     wd[i] = radians(wd[i])*-1.
    # """wind rose"""
    # plt.figure(8)
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf, width=width, bottom=bottom, color='blue',alpha = 0.25)
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Novato Rose', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, alpha = 0.25, color='red')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=18,family='serif')
    # plt.title('Novato Speeds', y=1.09,fontsize=18,family='serif')
    # plt.tight_layout()
    # plt.savefig('novatoRose_r.pdf',transparent=True)
    #
    #
    #
    #
    # plt.show()
