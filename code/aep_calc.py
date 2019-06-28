import numpy as np
import gaus_aep
import sys
sys.dont_write_bytecode = True


def fast_calc_AEP(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
            windSpeeds, windFrequencies, rated_ws=10., rated_power=3.35,
            cut_in_speed=3., cut_out_speed=25., shearExp=0.15, relaxationFactor=1.0, zref=99999., z0=0.):

    if zref == 99999.:
        zref = turbineZ[0]

    AEP = gaus_aep.calcaep(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies, shearExp, relaxationFactor, rated_ws, rated_power,
                cut_in_speed, cut_out_speed, zref, z0)

    return AEP


def fast_calc_gradients(turbineX, turbineY, turbineZ, rotorDiameter, windDirections,
            windSpeeds, windFrequencies, rated_ws=10., rated_power=3.35,
            cut_in_speed=3., cut_out_speed=25., shearExp=0.15, relaxationFactor=1.0, zref=99999., z0=0.):

    if zref == 99999.:
        zref = turbineZ[0]

    nTurbs = len(turbineX)

    dx = np.eye(nTurbs)
    dy = np.zeros((nTurbs,nTurbs))
    aep, daep_dx = gaus_aep.calcaep_dv(turbineX, dx, turbineY, dy, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies, shearExp, relaxationFactor, rated_ws, rated_power,
                cut_in_speed, cut_out_speed, zref, z0)

    dx = np.zeros((nTurbs,nTurbs))
    dy = np.eye(nTurbs)
    aep, daep_dy = gaus_aep.calcaep_dv(turbineX, dx, turbineY, dy, turbineZ, rotorDiameter, windDirections,
                windSpeeds, windFrequencies, shearExp, relaxationFactor, rated_ws, rated_power,
                cut_in_speed, cut_out_speed, zref, z0)

    return aep, daep_dx, daep_dy


# def fast_calc_grid_gradients(dx, dy, shear, rotate, turbs_per_row, x_start, y0,
#             turbineZ, rotorDiameter, windDirections,
#             windSpeeds, windFrequencies, rated_ws=10., rated_power=3.35,
#             cut_in_speed=3., cut_out_speed=25., shearExp=0.15, relaxationFactor=1.0, zref=99999., z0=0.):
#
#     dxd = 1.
#     dyd = 0.
#     sheard = 0.
#     rotated = 0.
#     aep, daep_ddx = gaus_aep_grid.calcaep_grid_dv(turbineZ,rotorDiameter,windDirections,dx,dxd,dy,dyd,shear,sheard,rotate,
#         rotated,turbs_per_row,x_start,y0,windSpeeds,windFrequencies,shearExp,relaxationFactor,
#         rated_ws,rated_power,cut_in_speed,cut_out_speed,zref,z0)
#
#     dxd = 0.
#     dyd = 1.
#     sheard = 0.
#     rotated = 0.
#     aep, daep_ddy = gaus_aep_grid.calcaep_grid_dv(turbineZ,rotorDiameter,windDirections,dx,dxd,dy,dyd,shear,sheard,rotate,
#         rotated,turbs_per_row,x_start,y0,windSpeeds,windFrequencies,shearExp,relaxationFactor,
#         rated_ws,rated_power,cut_in_speed,cut_out_speed,zref,z0)
#
#     dxd = 0.
#     dyd = 0.
#     sheard = 1.
#     rotated = 0.
#     aep, daep_ddshear = gaus_aep_grid.calcaep_grid_dv(turbineZ,rotorDiameter,windDirections,dx,dxd,dy,dyd,shear,sheard,rotate,
#         rotated,turbs_per_row,x_start,y0,windSpeeds,windFrequencies,shearExp,relaxationFactor,
#         rated_ws,rated_power,cut_in_speed,cut_out_speed,zref,z0)
#
#     dxd = 0.
#     dyd = 0.
#     sheard = 0.
#     rotated = 1.
#     aep, daep_drotate = gaus_aep_grid.calcaep_grid_dv(turbineZ,rotorDiameter,windDirections,dx,dxd,dy,dyd,shear,sheard,rotate,
#         rotated,turbs_per_row,x_start,y0,windSpeeds,windFrequencies,shearExp,relaxationFactor,
#         rated_ws,rated_power,cut_in_speed,cut_out_speed,zref,z0)
#
#     return aep, daep_ddx, daep_ddy, daep_ddshear, daep_drotate
