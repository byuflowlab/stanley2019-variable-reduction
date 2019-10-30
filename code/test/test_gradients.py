import numpy as np
import sys
sys.dont_write_bytecode = True
sys.path.insert(1, '/fslhome/pjstanle/compute/reduction')
from aep_calc import *
import time

nTurbines = 25
turbineX = np.random.rand(nTurbines)*1200.
turbineY = np.random.rand(nTurbines)*1200.
# turbineX = np.array([0.,0.,0.,0.,500.,500.,500.,500.,1000.,1000.,1000.,1000.,1500.,1500.,1500.,1500.])*2./3.
# turbineY = np.array([0.,500.,1000.,1500.,0.,500.,1000.,1500.,0.,500.,1000.,1500.,0.,500.,1000.,1500.])*2./3.
turbineZ = np.random.rand(nTurbines)*100.
# turbineZ = np.ones(nTurbines)*100.
rotorDiameter = np.random.rand(nTurbines)*100.
# rotorDiameter = np.ones(nTurbines)*100.
# windDirections = np.array([10.])
# windSpeeds = np.array([10.])
# windFrequencies = np.array([1.])
windDirections = np.random.rand(16)*360.
windSpeeds = np.random.rand(16)*25.
windFrequencies = np.random.rand(16)

# N = 20
# start = time.time()
# for i in range(N):
AEP, d_dx, d_dy = fast_calc_gradients(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
            windSpeeds, windFrequencies)
# print (time.time()-start)/float(N)
# FD
fd_x = np.zeros(nTurbines)
fd_y = np.zeros(nTurbines)
d = 1E-5
for i in range(nTurbines):
    x = np.zeros(nTurbines)
    y = np.zeros(nTurbines)
    x[:] = turbineX
    y[:] = turbineY
    x[i] += d
    y[i] += d
    aep_x,_,_ = fast_calc_gradients(x, turbineY, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies)
    aep_y,_,_ = fast_calc_gradients(turbineX, y, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies)
    fd_x[i] = (aep_x-AEP)/d
    fd_y[i] = (aep_y-AEP)/d

print 'AEP: ', AEP
print 'd_dx: ', d_dx
print 'fd_x: ', fd_x
print 'd_dy: ', d_dy
print 'fd_y: ', fd_y

print d_dx-fd_x
print d_dy-fd_y
