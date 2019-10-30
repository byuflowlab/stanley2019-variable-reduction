
!global functions
subroutine WindFrame(nTurbines, wind_direction, turbineX, turbineY, turbineXw, turbineYw)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: wind_direction
    real(dp), dimension(nTurbines), intent(in) :: turbineX, turbineY

    ! out
    real(dp), dimension(nTurbines), intent(out) :: turbineXw, turbineYw

    ! local
    real(dp) :: windDirectionDeg, windDirectionRad
    real(dp), parameter :: pi = 3.141592653589793_dp, tol = 0.000001_dp

    windDirectionDeg = 270. - wind_direction
    if (windDirectionDeg < 0.) then
        windDirectionDeg = windDirectionDeg + 360.
    end if
    windDirectionRad = pi*windDirectionDeg/180.0

    turbineXw = turbineX*cos(-windDirectionRad)-turbineY*sin(-windDirectionRad)
    turbineYw = turbineX*sin(-windDirectionRad)+turbineY*cos(-windDirectionRad)

end subroutine WindFrame


subroutine Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1, y)
    !    This function produces the y and dy values for a hermite cubic spline
    !    interpolating between two end points with known slopes
    !
    !    :param x: x position of output y
    !    :param x0: x position of upwind endpoint of spline
    !    :param x1: x position of downwind endpoint of spline
    !    :param y0: y position of upwind endpoint of spline
    !    :param dy0: slope at upwind endpoint of spline
    !    :param y1: y position of downwind endpoint of spline
    !    :param dy1: slope at downwind endpoint of spline
    !
    !    :return: y: y value of spline at location x

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x, x0, x1, y0, dy0, y1, dy1

    ! out
    real(dp), intent(out) :: y !, dy_dx

    ! local
    real(dp) :: c3, c2, c1, c0

    ! initialize coefficients for parametric cubic spline
    c3 = (2.0_dp*(y1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         (2.0_dp*(y0))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + &
         (dy0)/(x0**2 - 2.0_dp*x0*x1 + x1**2) + &
         (dy1)/(x0**2 - 2.0_dp*x0*x1 + x1**2)

    c2 = (3.0_dp*(y0)*(x0 + x1))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) - &
         ((dy1)*(2.0_dp*x0 + x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - ((dy0)*(x0 + &
         2.0_dp*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (3.0_dp*(y1)*(x0 + x1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)

    c1 = ((dy0)*(x1**2 + 2.0_dp*x0*x1))/(x0**2 - 2.0_dp*x0*x1 + x1**2) + ((dy1)*(x0**2 + &
         2.0_dp*x1*x0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - (6.0_dp*x0*x1*(y0))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3) + (6.0_dp*x0*x1*(y1))/(x0**3 - &
         3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - x1**3)

    c0 = ((y0)*(- x1**3 + 3.0_dp*x0*x1**2))/(x0**3 - 3.0_dp*x0**2*x1 + 3.0_dp*x0*x1**2 - &
         x1**3) - ((y1)*(- x0**3 + 3.0_dp*x1*x0**2))/(x0**3 - 3.0_dp*x0**2*x1 + &
         3.0_dp*x0*x1**2 - x1**3) - (x0*x1**2*(dy0))/(x0**2 - 2.0_dp*x0*x1 + x1**2) - &
         (x0**2*x1*(dy1))/(x0**2 - 2.0_dp*x0*x1 + x1**2)

    ! Solve for y and dy values at the given point
    y = c3*x**3 + c2*x**2 + c1*x + c0
    !dy_dx = c3*3*x**2 + c2*2*x + c1

end subroutine Hermite_Spline



!
! !power calculations"
subroutine PowWind(nTurbines, Uref, turbineZ, shearExp, zref, z0, &
                    &turbineSpeeds)

    implicit none
    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)
    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: Uref, shearExp, zref, z0
    real(dp), dimension(nTurbines), intent(in) :: turbineZ
    ! out
    real(dp), dimension(nTurbines), intent(out) :: turbineSpeeds
    ! local
    integer :: n

    do n = 1, nTurbines
        turbineSpeeds(n)= Uref*((turbineZ(n)-z0)/(zref-z0))**shearExp
    end do

end subroutine PowWind


subroutine DirPower(nTurbines, wtVelocity, rated_ws, rated_power, cut_in_speed, cut_out_speed,&
                        &dir_power)
    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: rated_ws, rated_power, cut_in_speed, cut_out_speed
    real(dp), dimension(nTurbines), intent(in) :: wtVelocity

    ! out
    real(dp), intent(out) :: dir_power

    ! local
    real(dp), dimension(nTurbines) :: wtPower
    real(dp) :: buffer, x0, x1, y0, y1, dy0, dy1
    integer :: n

    buffer = 0.1

    do n = 1, nTurbines
        ! If we're below cut-in
        if (wtVelocity(n) < (cut_in_speed-buffer)) then
            wtPower(n) = 0.
        ! If we're at the spline of cut-in
        else if (wtVelocity(n) > (cut_in_speed-buffer) .and. (wtVelocity(n) < (cut_in_speed+buffer))) then
            x0 = cut_in_speed-buffer
            x1 = cut_in_speed+buffer
            y0 = 0.
            y1 = rated_power*((cut_in_speed+buffer)/rated_ws)**3
            dy0 = 0.
            dy1 = 3.*rated_power*(cut_in_speed+buffer)**2/(rated_ws**3)
            call Hermite_Spline(wtVelocity(n), x0, x1, y0, dy0, y1, dy1, wtPower(n))
        ! If we're between cut-in and rated
        else if ((wtVelocity(n) > (cut_in_speed+buffer)) .and. (wtVelocity(n) < (rated_ws-buffer))) then
            wtPower(n) = rated_power*(wtVelocity(n)/rated_ws)**3
        ! If we're at the spline of rated
        else if ((wtVelocity(n) > (rated_ws-buffer)) .and. (wtVelocity(n) < (rated_ws+buffer))) then
            x0 = rated_ws-buffer
            x1 = rated_ws+buffer
            y0 = rated_power*((rated_ws-buffer)/rated_ws)**3
            y1 = rated_power
            dy0 = 3.*rated_power*(rated_ws-buffer)**2/(rated_ws**3)
            dy1 = 0.
            call Hermite_Spline(wtVelocity(n), x0, x1, y0, dy0, y1, dy1, wtPower(n))
        ! If we're between rated and cut-out
        else if ((wtVelocity(n) > (rated_ws+buffer)) .and. (wtVelocity(n) < (cut_out_speed-buffer))) then
            wtPower(n) = rated_power
        ! If we're at the spline of cut-out
        else if ((wtVelocity(n) > (cut_out_speed-buffer)) .and. (wtVelocity(n) < (cut_out_speed+buffer))) then
            x0 = cut_out_speed-buffer
            x1 = cut_out_speed+buffer
            y0 = rated_power
            y1 = 0.
            dy0 = 0.
            dy1 = 0.
            call Hermite_Spline(wtVelocity(n), x0, x1, y0, dy0, y1, dy1, wtPower(n))
        ! If we're above cut-out
        else if (wtVelocity(n) > (cut_out_speed+buffer)) then
            wtPower(n) = 0.
        end if

    end do

    dir_power = sum(wtPower)

end subroutine DirPower


subroutine calcAEP(nTurbines, nDirections, nRotorPoints, nCtPoints,&
            &turbineX, turbineY, turbineZ, rotorDiameter, Ct, yawDeg, windDirections,&
            &windSpeeds, windFrequencies, shearExp, rated_ws, rated_power,&
            &cut_in_speed, cut_out_speed, zref, z0,&
            &ky, kz, alpha, beta, TI, wec_factor, RotorPointsY, RotorPointsZ,&
            &ct_curve_wind_speed, ct_curve_ct, sm_smoothing, wake_combination_method, TI_calculation_method,&
            &wake_model_version, interp_type, calc_k_star, print_ti, use_ct_curve,AEP)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines, nDirections, nRotorPoints, nCtPoints
    real(dp), dimension(nTurbines), intent(in) :: turbineX, turbineY, turbineZ, rotorDiameter, Ct, yawDeg
    real(dp), dimension(nDirections), intent(in) :: windDirections, windSpeeds, windFrequencies
    real(dp), intent(in) :: shearExp, rated_ws, rated_power, cut_in_speed, cut_out_speed, zref, z0
    real(dp), intent(in) :: ky, kz, alpha, beta, TI, wec_factor
    real(dp), dimension(nRotorPoints), intent(in) :: RotorPointsY, RotorPointsZ
    real(dp), dimension(nCtPoints), intent(in) :: ct_curve_wind_speed, ct_curve_ct
    real(dp), intent(in) :: sm_smoothing
    integer, intent(in) :: wake_combination_method, TI_calculation_method, &
                        &  wake_model_version, interp_type
    logical, intent(in) :: calc_k_star, print_ti, use_ct_curve

    ! out
    real(dp), intent(out) :: AEP

    ! local
    real(dp), dimension(nDirections) :: dir_powers
    real(dp), dimension(nTurbines) :: turbineXw, turbineYw, Vinf_floris, wtVelocity, loss
    real(dp) :: hrs_per_year, pwrDir, Vinf
    integer :: n, i
    integer, dimension(nTurbines) :: sorted_x_idx

    do n = 1, nDirections
        call WindFrame(nTurbines, windDirections(n), turbineX, turbineY, turbineXw, turbineYw)
        call PowWind(nTurbines, windSpeeds(n), turbineZ, shearExp, zref, z0, Vinf_floris)
        Vinf = Vinf_floris(1)
        call sort_turbs(nTurbines,turbineXw,sorted_x_idx)
        ! call GaussianWake(nTurbines, turbineXw, turbineYw, rotorDiameter(1), relaxationFactor, loss)
        ! wtVelocity = Vinf*(1.0_dp-loss)

        if (Vinf > cut_in_speed) then
            call porteagel_analyze(nTurbines, nRotorPoints, nCtPoints, turbineXw, &
                                       sorted_x_idx, turbineYw, turbineZ, &
                                       rotorDiameter, Ct, Vinf, &
                                       yawDeg, ky, kz, alpha, beta, TI, RotorPointsY, RotorPointsZ, &
                                       zref, z0, shearExp, wake_combination_method, &
                                       TI_calculation_method, calc_k_star, wec_factor, print_ti, &
                                       wake_model_version, interp_type, &
                                       use_ct_curve, ct_curve_wind_speed, ct_curve_ct, sm_smoothing, &
                                       wtVelocity)
        else
            do i = 1, nTurbines
                wtVelocity(i) = 0.0
            end do
        end if
        ! print *, wtVelocity
        call DirPower(nTurbines, wtVelocity, rated_ws, rated_power, cut_in_speed, cut_out_speed, pwrDir)
        dir_powers(n) = pwrDir
    end do

    hrs_per_year = 365.*24.
    AEP = hrs_per_year * (sum(windFrequencies * dir_powers))

end subroutine calcAEP


subroutine sort_turbs(nTurbines,turbineX,sorted_x_idx)

      implicit none
      ! define precision to be the standard for a double precision ! on local system
      integer, parameter :: dp = kind(0.d0)

      ! in
      integer, intent(in) :: nTurbines
      real(dp), dimension(nTurbines), intent(in) :: turbineX

      ! out
      integer, dimension(nTurbines), intent(out) :: sorted_x_idx

      ! local
      integer :: i, j, index
      real(dp) :: min
      real(dp), dimension(nTurbines) :: sorting_turbines


      do i = 1, nTurbines
          sorting_turbines(i) = turbineX(i)
      end do

      do i = 1, nTurbines
          index = 1
          min = 9999999.0
          do j = 1, nTurbines
              if (sorting_turbines(j) < min) then
                  min = sorting_turbines(j)
                  index = j
              end if
          end do
          sorted_x_idx(i) = index-1
          sorting_turbines(index) = 99999999.0
      end do

end subroutine sort_turbs



! implementation of the Bastankhah and Porte Agel (BPA) wake model for analysis
subroutine porteagel_analyze(nTurbines, nRotorPoints, nCtPoints, turbineXw, &
                             sorted_x_idx, turbineYw, turbineZ, &
                             rotorDiameter, Ct, wind_speed, &
                             yawDeg, ky, kz, alpha, beta, TI, RotorPointsY, RotorPointsZ, &
                             z_ref, z_0, shear_exp, wake_combination_method, &
                             TI_calculation_method, calc_k_star, wec_factor, print_ti, &
                             wake_model_version, interp_type, &
                             use_ct_curve, ct_curve_wind_speed, ct_curve_ct, sm_smoothing, &
                             wtVelocity)

    ! independent variables: turbineXw turbineYw turbineZ rotorDiameter Ct yawDeg

    ! dependent variables: wtVelocity


    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines, nRotorPoints, nCtPoints
    integer, intent(in) :: wake_combination_method, TI_calculation_method, &
                        &  wake_model_version, interp_type
    logical, intent(in) :: calc_k_star, print_ti, use_ct_curve
    real(dp), dimension(nTurbines), intent(in) :: turbineXw, turbineYw, turbineZ
    integer, dimension(nTurbines), intent(in) :: sorted_x_idx
    real(dp), dimension(nTurbines), intent(in) :: rotorDiameter, yawDeg
    real(dp), dimension(nTurbines) :: Ct
    real(dp), intent(in) :: ky, kz, alpha, beta, TI, wind_speed, z_ref, z_0, shear_exp, wec_factor
    real(dp), dimension(nRotorPoints), intent(in) :: RotorPointsY, RotorPointsZ
    real(dp), dimension(nCtPoints), intent(in) :: ct_curve_wind_speed, ct_curve_ct
    real(dp), intent(in) :: sm_smoothing

    ! local (General)
    real(dp), dimension(nTurbines) :: yaw, TIturbs, Ct_local, ky_local, kz_local
    real(dp) :: x0, deltax0, deltay, theta_c_0, sigmay, sigmaz, wake_offset, k_star
    real(dp) :: x, deltav, deltaz, sigmay_dp, sigmaz_dp, deltax0_dp, deficit_sum
    real(dp) :: tol, discontinuity_point, TI_area_ratio
    real(dp) :: TI_area_ratio_tmp, TI_dst_tmp, TI_ust_tmp, rpts
    real(dp) :: LocalRotorPointY, LocalRotorPointZ, point_velocity, point_z, point_velocity_with_shear
    Integer :: u, d, turb, turbI, p
    real(dp), parameter :: pi = 3.141592653589793_dp

    ! model out
    real(dp), dimension(nTurbines), intent(out) :: wtVelocity

    intrinsic sin, cos, atan, max, sqrt, log

    ! bastankhah and porte agel 2016 define yaw to be positive clockwise, this is reversed
    yaw = - yawDeg*pi/180.0_dp

    ! set tolerance for location checks
    tol = 0.1_dp

    ! initialize wind turbine velocities to 0.0
    wtVelocity = 0.0_dp

    ! initialize TI of all turbines to free-stream value
    !print *, "start TIturbs: ", TIturbs
    TIturbs = TI


    ! initialize the local wake factors
    ky_local(:) = ky
    kz_local(:) = kz
    Ct_local(:) = Ct


    !print *, 'wake model version: ', wake_model_version

    !print *, "ky_local: ", ky_local
    !print *, "kz_local: ", kz_local
    !print *, "TIturbs init: ", TIturbs

    do, d=1, nTurbines

        ! get index of downstream turbine
        turbI = sorted_x_idx(d) + 1

        do, p=1, nRotorPoints

            ! initialize the TI_area_ratio to 0.0 for each turbine
            TI_area_ratio = 0.0_dp

            ! initialize deficit summation term to zero
            deficit_sum = 0.0_dp

            ! scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            LocalRotorPointY = RotorPointsY(p)*0.5_dp*rotorDiameter(turbI)
            LocalRotorPointZ = RotorPointsZ(p)*0.5_dp*rotorDiameter(turbI)
!             !print *, "rotorDiameter after local rotor points", rotorDiameter
!             !print *, "local rotor points Y,Z: ", LocalRotorPointY, LocalRotorPointZ

            do, u=1, nTurbines ! at turbineX-locations

                ! get index of upstream turbine
                turb = sorted_x_idx(u) + 1

                ! skip this loop if turb = turbI (turbines impact on itself)
                if (turb .eq. turbI) cycle

                ! downstream distance between upstream turbine and point
                x = turbineXw(turbI) - turbineXw(turb) + LocalRotorPointY*sin(yaw(turbI))

                ! set this iterations velocity deficit to 0
                deltav = 0.0_dp

                ! check turbine relative locations
                if (x > (0.0_dp + tol)) then

                    !print *, "rotorDiameter before x0 ", rotorDiameter

                    ! determine the onset location of far wake
                    call x0_func(rotorDiameter(turb), yaw(turb), Ct_local(turb), alpha, &
                                & TIturbs(turb), beta, x0)
!                     call x0_func(rotorDiameter(turb), yaw(turb), Ct(turb), alpha, &
!                                 & TI, beta, x0)

                    ! downstream distance from far wake onset to downstream turbine
                    deltax0 = x - x0

                    ! calculate wake spreading parameter at each turbine if desired
                    if (calc_k_star .eqv. .true.) then
                        call k_star_func(TIturbs(turb), k_star)
                        ky_local(turb) = k_star
                        kz_local(turb) = k_star
                    end if

                    !print *, "ky_local ", ky_local
                    !print *, "deltax0 ", deltax0
                    !print *, "turbineZ ", turbineZ
                    !print *, "rotorDiameter after x0 ", rotorDiameter
                    !print *, "Ct ", Ct
                    !print *, "yaw ", yaw

                    ! determine the initial wake angle at the onset of far wake
                    call theta_c_0_func(yaw(turb), Ct_local(turb), theta_c_0)
                    !print *, "theta_c_0 ", theta_c_0
                    ! horizontal spread
                    call sigmay_func(ky_local(turb), deltax0, rotorDiameter(turb), yaw(turb), sigmay)
                    !print *, "sigmay ", sigmay
                    !print *, "rotorDiameter after sigmay", rotorDiameter
                    ! vertical spread
                    call sigmaz_func(kz_local(turb), deltax0, rotorDiameter(turb), sigmaz)
                    !print *, "sigmaz ", sigmaz
                    !print *, "rotorDiameter after sigmaz ", rotorDiameter
                    ! horizontal cross-wind wake displacement from hub
                    call wake_offset_func(rotorDiameter(turb), theta_c_0, x0, yaw(turb), &
                                         & ky_local(turb), kz_local(turb), Ct_local(turb), sigmay, sigmaz, wake_offset)


                    !print *, "wake_offset ", wake_offset
                    ! cross wind distance from downstream point location to wake center
                    deltay = LocalRotorPointY*cos(yaw(turbI)) + turbineYw(turbI) - (turbineYw(turb) + wake_offset)

                    ! cross wind distance from hub height to height of point of interest
                    deltaz = LocalRotorPointZ + turbineZ(turbI) - turbineZ(turb)

                    !print *, "dx, dy, dz: ", x, deltay, deltaz
                    !print *, "local y,z : ", LocalRotorPointY, LocalRotorPointZ, turb, turbI, p
                    !print *, deltaz, deltay
                    ! far wake region

                    ! find the final point where the original model is undefined
                    call discontinuity_point_func(x0, rotorDiameter(turb), ky_local(turb), &
                                                 & kz_local(turb), yaw(turb), Ct_local(turb), &
                                                 & discontinuity_point)
                    !print *, "discontinuity point is: ", discontinuity_point
                    if (x > discontinuity_point) then

                        !print *, x

                        ! velocity difference in the wake
                        call deltav_func(deltay, deltaz, Ct_local(turb), yaw(turb), &
                                        & sigmay, sigmaz, rotorDiameter(turb), &
                                        & wake_model_version, kz_local(turb), x, &
                                        & wec_factor, deltav)
                        !print *, "rotorDiameter after far deltav ", rotorDiameter
                    ! near wake region (linearized)
                    else

                        ! determine distance from discontinuity point to far wake onset
                        deltax0_dp = discontinuity_point - x0

                        ! horizontal spread at far wake onset
                        call sigmay_func(ky_local(turb), deltax0_dp, rotorDiameter(turb), yaw(turb), sigmay_dp)

                        ! vertical spread at far wake onset
                        call sigmaz_func(kz_local(turb), deltax0_dp, rotorDiameter(turb), sigmaz_dp)

                       !  print *, "inputs in parent: ", deltay, deltaz, Ct(turb), yaw(turb), sigmay_dp, sigmaz_dp, &
!                                          & rotorDiameter(turb), x, discontinuity_point, sigmay_dp, sigmaz_dp, &
!                                          & wake_model_version, kz_local, x0, &
!                                          & wec_factor

                        ! velocity deficit in the nearwake (linear model)
                        call deltav_near_wake_lin_func(deltay, deltaz, &
                                         & Ct_local(turb), yaw(turb), sigmay_dp, sigmaz_dp, &
                                         & rotorDiameter(turb), x, discontinuity_point, sigmay_dp, sigmaz_dp, &
                                         & wake_model_version, kz_local(turb), x0, &
                                         & wec_factor, deltav)

                        !print *, "rotorDiameter after deltav near ", rotorDiameter
                    end if

                    ! combine deficits according to selected method wake combination method
                    call wake_combination_func(wind_speed, wtVelocity(turb), deltav,         &
                                               wake_combination_method, deficit_sum)

                    if ((x > 0.0_dp) .and. (TI_calculation_method > 0)) then
                        !print *, "turbI, turb: ", turbI, turb
                        ! calculate TI value at each turbine
!                         print *, "turb, turbI: ", turb, turbI

                        ! save ti_area_ratio and ti_dst to new memory locations to avoid
                        ! aliasing during differentiation
                        TI_area_ratio_tmp = TI_area_ratio
                        TI_dst_tmp = TIturbs(turbI)
                        TI_ust_tmp = TIturbs(turb)

                        call added_ti_func(TI, Ct_local(turb), x, ky_local(turb), rotorDiameter(turb), &
                                           & rotorDiameter(turbI), deltay, turbineZ(turb), &
                                           & turbineZ(turbI), sm_smoothing, TI_ust_tmp, &
                                           & TI_calculation_method, TI_area_ratio_tmp, &
                                           & TI_dst_tmp, TI_area_ratio, TIturbs(turbI))

                        !print *, "rotorDiameter after TI calcs", rotorDiameter
                    end if

!                     !print *, "deficit_sum, turbI, p, turb: ", deficit_sum, turbI, p, turb

                end if

            end do

            ! print *, deficit_sum

            ! find velocity at point p due to the wake of turbine turb
            point_velocity = wind_speed - deficit_sum

            !print *, "point velocity, deficit_sum, turbI, p: ", point_velocity, deficit_sum, turbI, p

            ! put sample point height in global reference frame
            point_z = LocalRotorPointZ + turbineZ(turbI)

            !print *, "point_z, turbI, p: ", point_z, turbI, p
            ! adjust sample point velocity for shear
            call wind_shear_func(point_z, point_velocity, z_ref, z_0, shear_exp, point_velocity_with_shear)
            !print *, "v, vs, x, turb, turbI, p: ", point_velocity, point_velocity_with_shear, x, turb, turbI, p
            ! add sample point velocity to turbine velocity to be averaged later
            wtVelocity(turbI) = wtVelocity(turbI) + point_velocity_with_shear

        end do

        ! final velocity calculation for turbine turbI (average equally across all points)
        rpts = REAL(nRotorPoints, dp)
!         print *, rpts, nRotorPoints, wtVelocity(turbI), wtVelocity(turbI)/rpts, wtVelocity(turbI)/nRotorPoints
!         STOP 1
        wtVelocity(turbI) = wtVelocity(turbI)/rpts
!         print *, wtVelocity(turbI)
        if (use_ct_curve) then
          ! print *, "wtVelocity(turbI): ", wtVelocity(turbI)
            call interpolation(nCtPoints, interp_type, ct_curve_wind_speed, ct_curve_ct, &
                              & wtVelocity(turbI), Ct_local(turbI), 0.0_dp, 0.0_dp, .false.)
          ! print *, "Ct_local(turbI): ", Ct_local(turbI)
        end if
        ! print *, Ct_local(turbI)
    end do

   !!  print TIturbs values to a file
!     if (print_ti) then
!         open(unit=2, file="TIturbs_tmp.txt")
!         do, turb=1, nTurbines
!             write(2,*) TIturbs(turb)
!         end do
!         close(2)
!     end if

    !print *, "TIturbs: ", TIturbs
    !print *, wtVelocity

    !! make sure turbine inflow velocity is non-negative
!             if (wtVelocity(turbI) .lt. 0.0_dp) then
!                 wtVelocity(turbI) = 0.0_dp
!             end if
    !print *, "fortran"

end subroutine porteagel_analyze


! calculates the onset of far-wake conditions
subroutine x0_func(rotor_diameter, yaw, Ct, alpha, TI, beta, x0)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: rotor_diameter, yaw, Ct, alpha, TI, beta

    ! out
    real(dp), intent(out) :: x0

    intrinsic cos, sqrt


    ! determine the onset location of far wake
    x0 = rotor_diameter * (cos(yaw) * (1.0_dp + sqrt(1.0_dp - Ct)) / &
                                (sqrt(2.0_dp) * (alpha * TI + beta * &
                                                & (1.0_dp - sqrt(1.0_dp - Ct)))))

end subroutine x0_func


! calculates the wake angle at the onset of far wake conditions
subroutine theta_c_0_func(yaw, Ct, theta_c_0)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: yaw, Ct

    ! out
    real(dp), intent(out) :: theta_c_0

    intrinsic cos, sqrt

    ! determine the initial wake angle at the onset of far wake
    theta_c_0 = 0.3_dp * yaw * (1.0_dp - sqrt(1.0_dp - Ct * cos(yaw))) / cos(yaw)

end subroutine theta_c_0_func


! calculates the horizontal spread of the wake at a given distance from the onset of far
! wake condition
subroutine sigmay_func(ky, deltax0, rotor_diameter, yaw, sigmay)

    implicit none

    ! define precision to be the standard for a double precision on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: ky, deltax0, rotor_diameter, yaw

    ! out
    real(dp), intent(out) :: sigmay

    intrinsic cos, sqrt

    ! horizontal spread
    sigmay = rotor_diameter * (ky * deltax0 / rotor_diameter + cos(yaw) / sqrt(8.0_dp))

end subroutine sigmay_func


! calculates the vertical spread of the wake at a given distance from the onset of far
! wake condition
subroutine sigmaz_func(kz, deltax0, rotor_diameter, sigmaz)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: kz, deltax0, rotor_diameter

    ! out
    real(dp), intent(out) :: sigmaz

    ! load necessary intrinsic functions
    intrinsic sqrt

    ! vertical spread
    sigmaz = rotor_diameter * (kz * deltax0 / rotor_diameter + 1.0_dp / sqrt(8.0_dp))

end subroutine sigmaz_func


! calculates the horizontal distance from the wake center to the hub of the turbine making
! the wake
subroutine wake_offset_func(rotor_diameter, theta_c_0, x0, yaw, ky, kz, Ct, sigmay, &
                            & sigmaz, wake_offset)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: rotor_diameter, theta_c_0, x0, yaw, ky, kz, Ct, sigmay
    real(dp), intent(in) :: sigmaz

    ! out
    real(dp), intent(out) :: wake_offset

    intrinsic cos, sqrt, log

    ! horizontal cross-wind wake displacement from hub
    wake_offset = rotor_diameter * (                                           &
                  theta_c_0 * x0 / rotor_diameter +                            &
                  (theta_c_0 / 14.7_dp) * sqrt(cos(yaw) / (ky * kz * Ct)) *    &
                  (2.9_dp + 1.3_dp * sqrt(1.0_dp - Ct) - Ct) *                 &
                  log(                                                         &
                    ((1.6_dp + sqrt(Ct)) *                                     &
                     (1.6_dp * sqrt(8.0_dp * sigmay * sigmaz /                 &
                                    (cos(yaw) * rotor_diameter ** 2))          &
                      - sqrt(Ct))) /                                           &
                    ((1.6_dp - sqrt(Ct)) *                                     &
                     (1.6_dp * sqrt(8.0_dp * sigmay * sigmaz /                 &
                                    (cos(yaw) * rotor_diameter ** 2))          &
                      + sqrt(Ct)))                                             &
                  )                                                            &
    )
end subroutine wake_offset_func


! calculates the velocity difference between hub velocity and free stream for a given wake
! for use in the far wake region
subroutine deltav_func(deltay, deltaz, Ct, yaw, sigmay, sigmaz, &
                      & rotor_diameter_ust, version, k, deltax, wec_factor, deltav)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: deltay, deltaz, Ct, yaw, sigmay
    real(dp), intent(in) :: sigmaz, rotor_diameter_ust, wec_factor
    real(dp), intent(in) :: k, deltax    ! only for 2014 version
    integer, intent(in) :: version

    ! local
    real(dp) :: beta_2014, epsilon_2014 ! only for 2014 version

    ! out
    real(dp), intent(out) :: deltav

    ! load intrinsic functions
    intrinsic cos, sqrt, exp

    !print *, "rotor_diameter in deltav entry", rotor_diameter_ust
!     print *, 'wake model version in deltav: ', version

    if (version == 2014) then
        !print *, "in 2014 version"
        beta_2014 = 0.5_dp*(1.0_dp + sqrt(1.0_dp - Ct))/sqrt(1.0_dp - Ct)
        epsilon_2014 = 0.2_dp*sqrt(beta_2014)

       ! print *, "beta = ", beta_2014, "epsilon = ", epsilon_2014
       ! print *, "k, deltax: ", k, deltax
       ! print *, "term: ", Ct                                                   &
!                            / (8.0_dp * (k*deltax/rotor_diameter_ust+epsilon_2014)**2)
        deltav = (                                                                       &
            (1.0_dp - sqrt(1.0_dp - Ct                                                   &
                           / (8.0_dp * ((k*deltax/rotor_diameter_ust)+epsilon_2014)**2)))* &
            exp((-1.0_dp/(2.0_dp*((k*deltax/rotor_diameter_ust) + epsilon_2014)**2))*      &
            ((deltaz/(wec_factor*rotor_diameter_ust))**2 + (deltay/(wec_factor*rotor_diameter_ust))**2))           &
        )
       ! print *, "deltav 2014 = ", deltav
    else if (version == 2016) then
        ! velocity difference in the wake at each sample point
        deltav = (                                                                    &
            (1.0_dp - sqrt(1.0_dp - Ct *                                                         &
                           cos(yaw) / (8.0_dp * sigmay * sigmaz / (rotor_diameter_ust ** 2)))) *     &
            exp(-0.5_dp * (deltay / (wec_factor*sigmay)) ** 2) * exp(-0.5_dp * (deltaz / (wec_factor*sigmaz)) ** 2)&
        )
    else
        print *, "Invalid Bastankhah and Porte Agel model version. Must be 2014 or 2016. ", version, " was given."
        stop 1
    end if

    !print *, "rotor_diameter in deltav exit", rotor_diameter_ust

end subroutine deltav_func


! calculates the velocity difference between hub velocity and free stream for a given wake
! for use in the near wake region only
subroutine deltav_near_wake_lin_func(deltay, deltaz, Ct, yaw,  &
                                 & sigmay, sigmaz, rotor_diameter_ust, x, &
                                 & discontinuity_point, sigmay0, sigmaz0, version, k, &
                                 & deltax0_dp, wec_factor, deltav)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: deltay, deltaz, Ct, yaw, sigmay
    real(dp), intent(in) :: sigmaz, rotor_diameter_ust, wec_factor
    real(dp), intent(in) :: x, discontinuity_point, sigmay0, sigmaz0
    real(dp), intent(in) :: k, deltax0_dp    ! only for 2014 version
    integer, intent(in) :: version

    ! local
    real(dp) :: deltav0m, deltavr
    real(dp) :: beta_2014, epsilon_2014 ! only for 2014 version

    ! out
    real(dp), intent(out) :: deltav

    ! load intrinsic functions
    intrinsic cos, sqrt, exp

   !  print *, 'wake model version in deltav near wake: ', version
!     print *, "inputs: ", deltay, deltaz, Ct, yaw,  &
!                                  & sigmay, sigmaz, rotor_diameter_ust, x, &
!                                  & discontinuity_point, sigmay0, sigmaz0, version, k, &
!                                  & deltax0_dp, wec_factor
    if (version == 2014) then
        if (yaw > 0.0_dp) then
            print *, "model version 2014 may only be used when yaw=0"
            stop 1
        end if
        beta_2014 = 0.5_dp*(1.0_dp + sqrt(1.0_dp - Ct))/sqrt(1.0_dp - Ct)
        epsilon_2014 = 0.2_dp*sqrt(beta_2014)

        ! magnitude term of gaussian at x0
        deltav0m = (1.0_dp - sqrt(1.0_dp - Ct                                            &
                           / (8.0_dp * (k*deltax0_dp/rotor_diameter_ust+epsilon_2014)**2)))

        ! initialize the gaussian magnitude term at the rotor for the linear interpolation
        deltavr = deltav0m

        ! linearized gaussian magnitude term for near wake
        deltav = (                                                                       &
             (((deltav0m - deltavr)/discontinuity_point) * x + deltavr) *                &
            exp((-1.0_dp/(2.0_dp*(k*deltax0_dp/rotor_diameter_ust + epsilon_2014)**2))*      &
            ((deltaz/(wec_factor*rotor_diameter_ust))**2 + (deltay/(wec_factor*rotor_diameter_ust))**2))           &
        )
    else if (version == 2016) then

        ! magnitude term of gaussian at x0
        deltav0m = (                                         &
                    (1.0_dp - sqrt(1.0_dp - Ct *                          &
                    cos(yaw) / (8.0_dp * sigmay0 * sigmaz0 /              &
                                                (rotor_diameter_ust ** 2)))))
        ! initialize the gaussian magnitude term at the rotor for the linear interpolation
        deltavr = deltav0m

        ! linearized gaussian magnitude term for near wake
        deltav = (((deltav0m - deltavr)/discontinuity_point) * x + deltavr) *       &
            exp(-0.5_dp * (deltay / (wec_factor*sigmay)) ** 2) *                                 &
            exp(-0.5_dp * (deltaz / (wec_factor*sigmaz)) ** 2)
    else
        print *, "Invalid Bastankhah and Porte Agel model version. Must be 2014 or 2016. ", version, " was given."
        stop 1
    end if

end subroutine deltav_near_wake_lin_func

! calculates the overlap area between a given wake and a rotor area
subroutine overlap_area_func(turbine_y, turbine_z, rotor_diameter, &
                            wake_center_y, wake_center_z, wake_diameter, &
                            wake_overlap)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: turbine_y, turbine_z, rotor_diameter
    real(dp), intent(in) :: wake_center_y, wake_center_z, wake_diameter

    ! out
    real(dp), intent(out) :: wake_overlap

    ! local
    real(dp), parameter :: pi = 3.141592653589793_dp, tol = 0.000001_dp
    real(dp) :: OVdYd, OVr, OVRR, OVL, OVz, OVz2

    ! load intrinsic functions
    intrinsic acos, sqrt

!     print *, turbine_y, turbine_z, rotor_diameter, &
!                             wake_center_y, wake_center_z, wake_diameter, &
!                             wake_overlap

   ! distance between wake center and rotor center
    if ((wake_center_z > (turbine_z + tol)) .or. (wake_center_z < (turbine_z - tol))) then
        OVdYd = sqrt((wake_center_y-turbine_y)**2_dp + (wake_center_z - turbine_z)**2_dp)
    else if (wake_center_y > (turbine_y + tol)) then! potential source of gradient issues, abs() did not cause a problem in FLORIS
        OVdYd = wake_center_y - turbine_y
    else if (turbine_y > (wake_center_y + tol)) then
        OVdYd = turbine_y - wake_center_y
    else
        OVdYd = 0.0_dp
    end if

    !print *, "OVdYd: ", OVdYd
    ! find rotor radius
    OVr = rotor_diameter/2.0_dp
    !print *, "OVr: ", OVr

    ! find wake radius
    OVRR = wake_diameter/2.0_dp
    !print *, "OVRR: ", OVRR

    ! determine if there is overlap
    if (OVdYd < (OVr+OVRR)) then ! if the rotor overlaps the wake zone

        ! check that turbine and wake centers are not perfectly aligned
        if (OVdYd > (0.0_dp + tol)) then

            ! check if the rotor is wholly contained in the wake
            if ((OVdYd + OVr) < OVRR + tol) then
                wake_overlap = pi*OVr*OVr
!                 print *, "1"
            ! check if the wake is wholly contained in the rotor swept area
            else if ((OVdYd + OVRR) < OVr + tol) then
                wake_overlap = pi*OVRR*OVRR
!                 print *, "2"
            else

                ! calculate the distance from the wake center to the chord connecting the lens
                ! cusps
                OVL = (-OVr*OVr+OVRR*OVRR+OVdYd*OVdYd)/(2.0_dp*OVdYd)

                OVz = sqrt(OVRR*OVRR-OVL*OVL)
                OVz2 = sqrt(OVr*OVr-(OVdYd-OVL)*(OVdYd-OVL))

                wake_overlap = OVRR*OVRR*acos(OVL/OVRR) + OVr*OVr*acos((OVdYd-OVL)/OVr) - &
                               & OVL*OVz - (OVdYd-OVL)*OVz2
!                 print *, OVRR, OVr, OVdYd, OVL, OVz, OVz2
!                 print *, "3"
            end if

        ! perfect overlap case where the wake is larger than the rotor
        else if (OVRR > OVr) then
            wake_overlap = pi*OVr*OVr
!             print *, "4"

        ! perfect overlap case where the rotor is larger than the wake
        else
            wake_overlap = pi*OVRR*OVRR
!             print *, "5"
        end if

    ! case with no overlap
    else
        wake_overlap = 0.0_dp
    end if

!     print *, "wake overlap in func: ", wake_overlap/(pi*OVr**2)
!     print *, "wake overlap in func: ", wake_overlap/(pi*OVRR**2)

    if ((wake_overlap/(pi*OVr*OVr) > (1.0_dp + tol)) .or. (wake_overlap/(pi*OVRR*OVRR) > (1.0_dp + tol))) then
        print *, "wake overlap in func: ", wake_overlap/(pi*OVr*OVr)
        print *, "wake overlap in func: ", wake_overlap/(pi*OVRR*OVRR)
        STOP 1
    end if

end subroutine overlap_area_func

! combines wakes using various methods
subroutine wake_combination_func(wind_speed, turb_inflow, deltav,                  &
                                 wake_combination_method, deficit_sum)

    ! combines wakes to calculate velocity at a given turbine
    ! wind_speed                = Free stream velocity
    ! turb_inflow               = Effective velocity as seen by the upstream rotor
    ! deltav                    = Velocity deficit percentage for current turbine pair
    ! wake_combination_method   = Use for selecting which method to use for wake combo
    ! deficit_sum (in)          = Combined deficits prior to including the current deltav
    ! deficit_sum (out)         = Combined deficits after to including the current deltav

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: wind_speed, turb_inflow, deltav
    integer, intent(in) :: wake_combination_method

    ! out
    real(dp), intent(inout) :: deficit_sum

    ! intrinsic functions
    intrinsic sqrt

    ! freestream linear superposition (Lissaman 1979)
    if (wake_combination_method == 0) then
        deficit_sum = deficit_sum + wind_speed*deltav

    ! local velocity linear superposition (Niayifar and Porte Agel 2015, 2016)
    else if (wake_combination_method == 1) then
        deficit_sum = deficit_sum + turb_inflow*deltav
        !print *, "here"

    ! sum of squares freestream superposition (Katic et al. 1986)
    else if (wake_combination_method == 2) then
        deficit_sum = sqrt(deficit_sum**2 + (wind_speed*deltav)**2)

    ! sum of squares local velocity superposition (Voutsinas 1990)
    else if (wake_combination_method == 3) then
        deficit_sum = sqrt(deficit_sum**2 + (turb_inflow*deltav)**2)

    ! wake combination method error
    else
        print *, "Invalid wake combination method. Must be one of [0,1,2,3]."
        stop 1
    end if

end subroutine wake_combination_func

! combines wakes using various methods
subroutine added_ti_func(TI, Ct_ust, x, k_star_ust, rotor_diameter_ust, rotor_diameter_dst, &
                        & deltay, wake_height, turbine_height, sm_smoothing, TI_ust, &
                        & TI_calculation_method, TI_area_ratio_in, TI_dst_in, TI_area_ratio, TI_dst)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: Ct_ust, x, k_star_ust, rotor_diameter_ust, rotor_diameter_dst
    real(dp), intent(in) :: deltay, wake_height, turbine_height, sm_smoothing
    real(dp), intent(in) :: TI_ust, TI, TI_area_ratio_in, TI_dst_in
    integer, intent(in) :: TI_calculation_method

    ! local
    real(dp) :: axial_induction_ust, beta, epsilon, sigma, wake_diameter, wake_overlap
    real(dp) :: TI_added, TI_tmp, rotor_area_dst, TI_area_ratio_tmp
    real(dp), parameter :: pi = 3.141592653589793_dp

    ! out
    real(dp), intent(out) :: TI_dst, TI_area_ratio

    ! intrinsic functions
    intrinsic sqrt

    ! initialize output variables
    TI_area_ratio = TI_area_ratio_in
    TI_dst = TI_dst_in

    ! initialize wake overlap to zero
    wake_overlap = 0.0_dp

    !print *, "TI_dst in: ", TI_dst
    ! Niayifar and Porte Agel 2015, 2016 (adjusted by Annoni and Thomas for SOWFA match
    ! and optimization)
    if (TI_calculation_method == 1) then

        ! calculate axial induction based on the Ct value
        call ct_to_axial_ind_func(Ct_ust, axial_induction_ust)

        ! calculate BPA spread parameters Bastankhah and Porte Agel 2014
        beta = 0.5_dp*((1.0_dp + sqrt(1.0_dp - Ct_ust))/sqrt(1.0_dp - Ct_ust))
        epsilon = 0.2_dp*sqrt(beta)
        !print *, "epsilon = ", epsilon
        ! calculate wake spread for TI calcs
        sigma = k_star_ust*x + rotor_diameter_ust*epsilon
        wake_diameter = 4.0_dp*sigma
        !print *, "sigma = ", sigma
        ! calculate wake overlap ratio
        call overlap_area_func(deltay, turbine_height, rotor_diameter_dst, &
                            0.0_dp, wake_height, wake_diameter, &
                            wake_overlap)
        !print *, "wake_overlap = ", wake_overlap
        ! Calculate the turbulence added to the inflow of the downstream turbine by the
        ! wake of the upstream turbine
        TI_added = 0.73_dp*(axial_induction_ust**0.8325_dp)*(TI_ust**0.0325_dp)* &
                    ((x/rotor_diameter_ust)**(-0.32_dp))
        !print *, "TI_added = ", TI_added
        rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
        ! Calculate the total turbulence intensity at the downstream turbine
        !sum_of_squares = TI_dst**2 + (TI_added*wake_overlap)**2
        ! print *, "sum of squares = ", sum_of_squares
!         TI_dst = sqrt(sum_of_squares)
!         !print *, "TI_dst = ", TI_dst
        TI_dst = sqrt(TI_dst_in**2.0_dp + (TI_added*wake_overlap/rotor_area_dst)**2.0_dp)


    ! Niayifar and Porte Agel 2015, 2016
    else if (TI_calculation_method == 2) then

        ! calculate axial induction based on the Ct value
        call ct_to_axial_ind_func(Ct_ust, axial_induction_ust)

        ! calculate BPA spread parameters Bastankhah and Porte Agel 2014
        beta = 0.5_dp*((1.0_dp + sqrt(1.0_dp - Ct_ust))/sqrt(1.0_dp - Ct_ust))
        epsilon = 0.2_dp*sqrt(beta)

        ! calculate wake spread for TI calcs
        sigma = k_star_ust*x + rotor_diameter_ust*epsilon
        wake_diameter = 4.0_dp*sigma

        ! calculate wake overlap ratio
        call overlap_area_func(deltay, turbine_height, rotor_diameter_dst, &
                            0.0_dp, wake_height, wake_diameter, &
                            wake_overlap)

        ! Calculate the turbulence added to the inflow of the downstream turbine by the
        ! wake of the upstream turbine
        TI_added = 0.73_dp*(axial_induction_ust**0.8325_dp)*(TI_ust**0.0325_dp)* &
                    ((x/rotor_diameter_ust)**(-0.32_dp))

        ! Calculate the total turbulence intensity at the downstream turbine based on
        ! current upstream turbine
        rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
        TI_tmp = sqrt(TI**2.0_dp + (TI_added*(wake_overlap/rotor_area_dst))**2.0_dp)

        ! Check if this is the max and use it if it is
        if (TI_tmp > TI_dst_in) then
!            print *, "TI_tmp > TI_dst"
           TI_dst = TI_tmp
        end if

    ! Niayifar and Porte Agel 2015, 2016 with smooth max
     else if (TI_calculation_method == 3) then

        ! calculate axial induction based on the Ct value
        call ct_to_axial_ind_func(Ct_ust, axial_induction_ust)

        ! calculate BPA spread parameters Bastankhah and Porte Agel 2014
        beta = 0.5_dp*((1.0_dp + sqrt(1.0_dp - Ct_ust))/sqrt(1.0_dp - Ct_ust))
        epsilon = 0.2_dp*sqrt(beta)

        ! calculate wake spread for TI calcs
        sigma = k_star_ust*x + rotor_diameter_ust*epsilon
        wake_diameter = 4.0_dp*sigma

!         print *, "sigma, k_star_ust, x, rotor_diameter_ust, epsilon ", sigma, k_star_ust, x, rotor_diameter_ust, epsilon

        ! print *, "deltay, turbine_height, rotor_diameter_dst, wake_height, wake_diameter", &
!                 & deltay, turbine_height, rotor_diameter_dst, &
!                             wake_height, wake_diameter

        ! calculate wake overlap ratio
        call overlap_area_func(deltay, turbine_height, rotor_diameter_dst, &
                            0.0_dp, wake_height, wake_diameter, &
                            wake_overlap)

        ! Calculate the turbulence added to the inflow of the downstream turbine by the
        ! wake of the upstream turbine
        TI_added = 0.73_dp*(axial_induction_ust**0.8325_dp)*(TI_ust**0.0325_dp)* &
                    ((x/rotor_diameter_ust)**(-0.32_dp))

        ! Calculate the total turbulence intensity at the downstream turbine based on
        ! current upstream turbine
        rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
        TI_tmp = sqrt(TI**2.0_dp + (TI_added*(wake_overlap/rotor_area_dst))**2.0_dp)

        !print *, "TI, TI_added, wake_overlap, rotor_area_dst: ", TI, TI_added, wake_overlap, rotor_area_dst

        ! Check if this is the max and use it if it is
        !if (TI_tmp > TI_dst) then
        !    TI_dst = TI_tmp
        !end if
!         print *, "before: ", TI_dst, TI_tmp
!         TI_dst_in = TI_dst
        call smooth_max(sm_smoothing, TI_dst_in, TI_tmp, TI_dst)
!         print *, "after:: ", TI_dst, TI_tmp

    ! Niayifar and Porte Agel 2015, 2016 using max on area TI ratio
    else if (TI_calculation_method == 4) then

        ! calculate axial induction based on the Ct value
        call ct_to_axial_ind_func(Ct_ust, axial_induction_ust)

        ! calculate BPA spread parameters Bastankhah and Porte Agel 2014
        beta = 0.5_dp*((1.0_dp + sqrt(1.0_dp - Ct_ust))/sqrt(1.0_dp - Ct_ust))
        epsilon = 0.2_dp*sqrt(beta)

        ! calculate wake spread for TI calcs
        sigma = k_star_ust*x + rotor_diameter_ust*epsilon
        wake_diameter = 4.0_dp*sigma

        ! calculate wake overlap ratio
        call overlap_area_func(deltay, turbine_height, rotor_diameter_dst, &
                            0.0_dp, wake_height, wake_diameter, &
                            wake_overlap)

        ! Calculate the turbulence added to the inflow of the downstream turbine by the
        ! wake of the upstream turbine
        TI_added = 0.73_dp*(axial_induction_ust**0.8325_dp)*(TI_ust**0.0325_dp)* &
                    ((x/rotor_diameter_ust)**(-0.32_dp))

        ! Calculate the total turbulence intensity at the downstream turbine based on
        ! current upstream turbine
        rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
        TI_area_ratio_tmp = TI_added*(wake_overlap/rotor_area_dst)

        ! Check if this is the max and use it if it is
        if (TI_area_ratio_tmp > TI_area_ratio_in) then
!            print *, "ti_area_ratio_tmp > ti_area_ratio"
           !TI_dst = TI_tmp
           TI_area_ratio = TI_area_ratio_tmp
           TI_dst = sqrt(TI**2.0_dp + (TI_area_ratio)**2.0_dp)

        end if

    ! Niayifar and Porte Agel 2015, 2016 using smooth max on area TI ratio
    else if (TI_calculation_method == 5) then

        ! calculate axial induction based on the Ct value
        call ct_to_axial_ind_func(Ct_ust, axial_induction_ust)

        ! calculate BPA spread parameters Bastankhah and Porte Agel 2014
        beta = 0.5_dp*((1.0_dp + sqrt(1.0_dp - Ct_ust))/sqrt(1.0_dp - Ct_ust))
        epsilon = 0.2_dp*sqrt(beta)

        ! calculate wake spread for TI calcs
        sigma = k_star_ust*x + rotor_diameter_ust*epsilon
        wake_diameter = 4.0_dp*sigma

        ! calculate wake overlap ratio
        call overlap_area_func(deltay, turbine_height, rotor_diameter_dst, &
                            0.0_dp, wake_height, wake_diameter, &
                             wake_overlap)
        ! only include turbines with area overlap in the softmax
        if (wake_overlap > 0.0_dp) then

            ! Calculate the turbulence added to the inflow of the downstream turbine by the
            ! wake of the upstream turbine
            TI_added = 0.73_dp*(axial_induction_ust**0.8325_dp)*(TI_ust**0.0325_dp)* &
                        ((x/rotor_diameter_ust)**(-0.32_dp))


            rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
            TI_area_ratio_tmp = TI_added*(wake_overlap/rotor_area_dst)
            !TI_tmp = sqrt(TI**2.0_dp + (TI_added*(wake_overlap/rotor_area_dst))**2.0_dp)

            ! Run through the smooth max to get an approximation of the true max TI area ratio
            call smooth_max(sm_smoothing, TI_area_ratio_in, TI_area_ratio_tmp, TI_area_ratio)

            ! Calculate the total turbulence intensity at the downstream turbine based on
            ! the result of the smooth max function
            TI_dst = sqrt(TI**2.0_dp + TI_area_ratio**2.0_dp)

        end if

    ! wake combination method error
    else
        print *, "Invalid added TI calculation method. Must be one of [0,1,2,3,4,5]."
        stop 1
    end if

    !print *, "ratio: ", wake_overlap/rotor_area_dst
    !print *, "Dr, Dw: ", rotor_diameter_dst, wake_diameter
    !print *, "Ar, Aol: ", rotor_area_dst, wake_overlap

end subroutine added_ti_func

! compute wake spread parameter based on local turbulence intensity
subroutine k_star_func(TI_ust, k_star_ust)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: TI_ust

    ! out
    real(dp), intent(out) :: k_star_ust

    ! calculate wake spread parameter from Niayifar and Porte Agel (2015, 2016)
    k_star_ust = 0.3837*TI_ust+0.003678

end subroutine k_star_func

! calculate axial induction from Ct
subroutine ct_to_axial_ind_func(CT, axial_induction)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: CT

    ! out
    real(dp), intent(out) :: axial_induction

    ! initialize axial induction to zero
    axial_induction = 0.0_dp

    ! calculate axial induction
    if (CT > 0.96) then  ! Glauert condition
        axial_induction = 0.143_dp + sqrt(0.0203_dp-0.6427_dp*(0.889_dp - CT))
    else
        axial_induction = 0.5_dp*(1.0_dp-sqrt(1.0_dp-CT))
    end if

end subroutine ct_to_axial_ind_func

! adjust wind speed for wind shear
subroutine wind_shear_func(point_z, u_ref, z_ref, z_0, shear_exp, adjusted_wind_speed)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: point_z, u_ref, z_ref, z_0, shear_exp

    ! out
    real(dp), intent(out) :: adjusted_wind_speed

    ! initialize adjusted wind speed to zero
    adjusted_wind_speed = 0.0_dp

    ! check that the point of interest is above ground level
    if (point_z >= z_0) then
        ! adjusted wind speed for wind shear if point is above ground
        adjusted_wind_speed = u_ref*((point_z-z_0)/(z_ref-z_0))**shear_exp
    else
        ! if the point of interest is below ground, set the wind speed to 0.0
        adjusted_wind_speed = 0.0_dp
    end if

end subroutine wind_shear_func


! calculate the point where the Bastankhah and Porte Agel wake model becomes undefined
subroutine discontinuity_point_func(x0, rotor_diameter, ky, kz, yaw, Ct, discontinuity_point)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x0, rotor_diameter, ky, kz, yaw, Ct

    ! local
    real(dp) :: a, b, c

    ! out
    real(dp), intent(out) :: discontinuity_point

    intrinsic cos, sqrt

    ! for clarity, break out the terms in the equation
    a = ky + kz*cos(yaw)
    b = 4.0_dp * ky * kz * cos(yaw)*(Ct - 1.0_dp)
    c = 2.0_dp * sqrt(8.0_dp) * ky * kz

    ! distance from rotor to the last point where the wake model is undefined
    discontinuity_point = x0 + rotor_diameter * (a - sqrt(a**2 - b))/c

end subroutine discontinuity_point_func

subroutine smooth_max(s, x, y, g)

    ! based on John D. Cook's writings at
    ! (1) https://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! and
    ! (2) https://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/

    ! s controls the level of smoothing used in the smooth max
    ! x and y are the values to be compared

    ! g is the result

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: s, x, y

    ! local
    real(dp) :: max_val, min_val

    ! out
    real(dp), intent(out) :: g

    intrinsic log, exp, max, min

    ! LogSumExponential Method - used this in the past
    ! g = (x*exp(s*x)+y*exp(s*y))/(exp(s*x)+exp(s*y))

    ! non-overflowing version of Smooth Max function (see ref 2 above)
    max_val = max(x, y)
    min_val = min(x, y)
    g = (log(1.0_dp + exp(s*(min_val - max_val)))+s*max_val)/s

end subroutine smooth_max

subroutine interpolation(nPoints, interp_type, x, y, xval, yval, dy0in, dy1in, usedyin)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nPoints, interp_type
    real(dp), dimension(nPoints), intent(in) :: x, y
    real(dp), intent(in) :: xval
    real(dp), intent(in):: dy0in, dy1in
    logical :: usedyin

    ! local
    integer :: idx
    real(dp) :: x0, x1, y0, dy0, y1, dy1

    ! out
    real(dp), intent(out) :: yval

!     print *, "in interpolation"

    ! if ((xval < x(1)) .or. (xval > x(nPoints))) then
!         print *, "interpolation point is out of bounds"
! !         STOP 1
!     end if

    if (usedyin .and. (interp_type == 1)) then
        print *, "end point derivatives may not be specified for linear interpolation"
        STOP 1
    end if

    if (xval < x(1)) then
        yval = y(1)
    else if (xval > x(nPoints)) then
        yval = y(nPoints)

    else
        idx = 1

        do while ((xval > x(idx)) .and. (idx <= nPoints))
            idx = idx + 1
        end do

        idx = idx - 1

        x0 = x(idx)
        x1 = x((idx + 1))
        y0 = y(idx)
        y1 = y((idx + 1))

        ! Hermite cubic piecewise spline interpolation
        if (interp_type == 0) then

            ! approximate derivative at left end of interval
            if (idx == 1) then

                if (usedyin) then
                    dy0 = dy0in
                else
                    dy0 = 0.0_dp
                endif

            else
                dy0 = (y(idx+1) - y(idx-1))/(x(idx+1) - x(idx-1))
            end if

            ! approximate derivative at the right end of interval
            if (idx >= nPoints-1) then

                if(usedyin)then
                    dy1 = dy1in
                else
                    dy1 = 0.0_dp
                endif
            else
                dy1 = (y(idx+2) - y(idx))/(x(idx+2) - x(idx))
            end if

            ! call Hermite spline routine
            call Hermite_Spline(xval, x0, x1, y0, dy0, y1, dy1, yval)

        ! linear interpolation
        else if (interp_type == 1) then
            yval = (xval-x0)*(y1-y0)/(x1-x0) + y0

        end if
    end if

!     print *, "yval = ", yval

end subroutine interpolation























!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (master-db54337a6) - 29 Jul 2019 10:54
!
!  Differentiation of calcaep in forward (tangent) mode:
!   variations   of useful results: aep
!   with respect to varying inputs: turbinex turbiney
!   RW status of diff variables: turbinex:in turbiney:in aep:out
SUBROUTINE CALCAEP_DV(nturbines, ndirections, nrotorpoints, nctpoints, &
& turbinex, turbinexd, turbiney, turbineyd, turbinez, rotordiameter, ct&
& , yawdeg, winddirections, windspeeds, windfrequencies, shearexp, &
& rated_ws, rated_power, cut_in_speed, cut_out_speed, zref, z0, ky, kz, &
& alpha, beta, ti, wec_factor, rotorpointsy, rotorpointsz, &
& ct_curve_wind_speed, ct_curve_ct, sm_smoothing, &
& wake_combination_method, ti_calculation_method, wake_model_version, &
& interp_type, calc_k_star, print_ti, use_ct_curve, aep, aepd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines, ndirections, nrotorpoints, nctpoints
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney, &
& turbinez, rotordiameter, ct, yawdeg
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbinexd, &
& turbineyd
  REAL(dp), DIMENSION(ndirections), INTENT(IN) :: winddirections, &
& windspeeds, windfrequencies
  REAL(dp), INTENT(IN) :: shearexp, rated_ws, rated_power, cut_in_speed&
& , cut_out_speed, zref, z0
  REAL(dp), INTENT(IN) :: ky, kz, alpha, beta, ti, wec_factor
  REAL(dp), DIMENSION(nrotorpoints), INTENT(IN) :: rotorpointsy, &
& rotorpointsz
  REAL(dp), DIMENSION(nctpoints), INTENT(IN) :: ct_curve_wind_speed, &
& ct_curve_ct
  REAL(dp), INTENT(IN) :: sm_smoothing
  INTEGER, INTENT(IN) :: wake_combination_method, ti_calculation_method&
& , wake_model_version, interp_type
  LOGICAL, INTENT(IN) :: calc_k_star, print_ti, use_ct_curve
! out
  REAL(dp), INTENT(OUT) :: aep
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: aepd
! local
  REAL(dp), DIMENSION(ndirections) :: dir_powers
  REAL(dp), DIMENSION(nbdirs, ndirections) :: dir_powersd
  REAL(dp), DIMENSION(nturbines) :: turbinexw, turbineyw, vinf_floris, &
& wtvelocity, loss
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbinexwd, turbineywd, &
& wtvelocityd
  REAL(dp) :: hrs_per_year, pwrdir, vinf
  REAL(dp), DIMENSION(nbdirs) :: pwrdird
  INTEGER :: n, i
  INTEGER, DIMENSION(nturbines) :: sorted_x_idx
  INTRINSIC SUM
  REAL(dp), DIMENSION(ndirections) :: arg1
  REAL(dp), DIMENSION(nbdirs, ndirections) :: arg1d
  INTEGER :: nd
  INTEGER :: nbdirs
  dir_powersd(:, :) = 0.0_8
  wtvelocityd(:, :) = 0.0_8
  DO n=1,ndirections
    CALL WINDFRAME_DV(nturbines, winddirections(n), turbinex, turbinexd&
&               , turbiney, turbineyd, turbinexw, turbinexwd, turbineyw&
&               , turbineywd, nbdirs)
    CALL POWWIND(nturbines, windspeeds(n), turbinez, shearexp, zref, z0&
&          , vinf_floris)
    vinf = vinf_floris(1)
    CALL SORT_TURBS(nturbines, turbinexw, sorted_x_idx)
! call GaussianWake(nTurbines, turbineXw, turbineYw, rotorDiameter(1), relaxationFactor, loss)
! wtVelocity = Vinf*(1.0_dp-loss)
    IF (vinf .GT. cut_in_speed) THEN
      CALL PORTEAGEL_ANALYZE_DV(nturbines, nrotorpoints, nctpoints, &
&                         turbinexw, turbinexwd, sorted_x_idx, turbineyw&
&                         , turbineywd, turbinez, rotordiameter, ct, &
&                         vinf, yawdeg, ky, kz, alpha, beta, ti, &
&                         rotorpointsy, rotorpointsz, zref, z0, shearexp&
&                         , wake_combination_method, &
&                         ti_calculation_method, calc_k_star, wec_factor&
&                         , print_ti, wake_model_version, interp_type, &
&                         use_ct_curve, ct_curve_wind_speed, ct_curve_ct&
&                         , sm_smoothing, wtvelocity, wtvelocityd, &
&                         nbdirs)
    ELSE
      DO i=1,nturbines
        DO nd=1,nbdirs
          wtvelocityd(nd, i) = 0.0_8
        END DO
        wtvelocity(i) = 0.0
      END DO
    END IF
! print *, wtVelocity
    CALL DIRPOWER_DV(nturbines, wtvelocity, wtvelocityd, rated_ws, &
&              rated_power, cut_in_speed, cut_out_speed, pwrdir, pwrdird&
&              , nbdirs)
    DO nd=1,nbdirs
      dir_powersd(nd, n) = pwrdird(nd)
    END DO
    dir_powers(n) = pwrdir
  END DO
  hrs_per_year = 365.*24.
  DO nd=1,nbdirs
    arg1d(nd, :) = windfrequencies*dir_powersd(nd, :)
    aepd(nd) = hrs_per_year*SUM(arg1d(nd, :))
  END DO
  arg1(:) = windfrequencies*dir_powers
  aep = hrs_per_year*SUM(arg1(:))
END SUBROUTINE CALCAEP_DV

!  Differentiation of windframe in forward (tangent) mode:
!   variations   of useful results: turbinexw turbineyw
!   with respect to varying inputs: turbinex turbiney
!global functions
SUBROUTINE WINDFRAME_DV(nturbines, wind_direction, turbinex, turbinexd, &
& turbiney, turbineyd, turbinexw, turbinexwd, turbineyw, turbineywd, &
& nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: wind_direction
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinex, turbiney
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbinexd, &
& turbineyd
! out
  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: turbinexw, turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: turbinexwd, &
& turbineywd
! local
  REAL(dp) :: winddirectiondeg, winddirectionrad
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  INTRINSIC COS
  INTRINSIC SIN
  INTEGER :: nd
  INTEGER :: nbdirs
  winddirectiondeg = 270. - wind_direction
  IF (winddirectiondeg .LT. 0.) winddirectiondeg = winddirectiondeg + &
&     360.
  winddirectionrad = pi*winddirectiondeg/180.0
  DO nd=1,nbdirs
    turbinexwd(nd, :) = COS(-winddirectionrad)*turbinexd(nd, :) - SIN(-&
&     winddirectionrad)*turbineyd(nd, :)
    turbineywd(nd, :) = SIN(-winddirectionrad)*turbinexd(nd, :) + COS(-&
&     winddirectionrad)*turbineyd(nd, :)
  END DO
  turbinexw = turbinex*COS(-winddirectionrad) - turbiney*SIN(-&
&   winddirectionrad)
  turbineyw = turbinex*SIN(-winddirectionrad) + turbiney*COS(-&
&   winddirectionrad)
END SUBROUTINE WINDFRAME_DV

!  Differentiation of dirpower in forward (tangent) mode:
!   variations   of useful results: dir_power
!   with respect to varying inputs: wtvelocity
SUBROUTINE DIRPOWER_DV(nturbines, wtvelocity, wtvelocityd, rated_ws, &
& rated_power, cut_in_speed, cut_out_speed, dir_power, dir_powerd, &
& nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: rated_ws, rated_power, cut_in_speed, &
& cut_out_speed
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: wtvelocity
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: wtvelocityd
! out
  REAL(dp), INTENT(OUT) :: dir_power
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: dir_powerd
! local
  REAL(dp), DIMENSION(nturbines) :: wtpower
  REAL(dp), DIMENSION(nbdirs, nturbines) :: wtpowerd
  REAL(dp) :: buffer, x0, x1, y0, y1, dy0, dy1
  INTEGER :: n
  INTRINSIC SUM
  INTEGER :: nd
  INTEGER :: nbdirs
  buffer = 0.1
  wtpowerd(:, :) = 0.0_8
  DO n=1,nturbines
! If we're below cut-in
    IF (wtvelocity(n) .LT. cut_in_speed - buffer) THEN
      DO nd=1,nbdirs
        wtpowerd(nd, n) = 0.0_8
      END DO
      wtpower(n) = 0.
! If we're at the spline of cut-in
    ELSE IF (wtvelocity(n) .GT. cut_in_speed - buffer .AND. wtvelocity(n&
&       ) .LT. cut_in_speed + buffer) THEN
      x0 = cut_in_speed - buffer
      x1 = cut_in_speed + buffer
      y0 = 0.
      y1 = rated_power*((cut_in_speed+buffer)/rated_ws)**3
      dy0 = 0.
      dy1 = 3.*rated_power*(cut_in_speed+buffer)**2/rated_ws**3
      CALL HERMITE_SPLINE_DV(wtvelocity(n), wtvelocityd(:, n), x0, x1, &
&                      y0, dy0, y1, dy1, wtpower(n), wtpowerd(:, n), &
&                      nbdirs)
! If we're between cut-in and rated
    ELSE IF (wtvelocity(n) .GT. cut_in_speed + buffer .AND. wtvelocity(n&
&       ) .LT. rated_ws - buffer) THEN
      DO nd=1,nbdirs
        wtpowerd(nd, n) = rated_power*3*wtvelocity(n)**2*wtvelocityd(nd&
&         , n)/rated_ws**3
      END DO
      wtpower(n) = rated_power*(wtvelocity(n)/rated_ws)**3
! If we're at the spline of rated
    ELSE IF (wtvelocity(n) .GT. rated_ws - buffer .AND. wtvelocity(n) &
&       .LT. rated_ws + buffer) THEN
      x0 = rated_ws - buffer
      x1 = rated_ws + buffer
      y0 = rated_power*((rated_ws-buffer)/rated_ws)**3
      y1 = rated_power
      dy0 = 3.*rated_power*(rated_ws-buffer)**2/rated_ws**3
      dy1 = 0.
      CALL HERMITE_SPLINE_DV(wtvelocity(n), wtvelocityd(:, n), x0, x1, &
&                      y0, dy0, y1, dy1, wtpower(n), wtpowerd(:, n), &
&                      nbdirs)
! If we're between rated and cut-out
    ELSE IF (wtvelocity(n) .GT. rated_ws + buffer .AND. wtvelocity(n) &
&       .LT. cut_out_speed - buffer) THEN
      DO nd=1,nbdirs
        wtpowerd(nd, n) = 0.0_8
      END DO
      wtpower(n) = rated_power
! If we're at the spline of cut-out
    ELSE IF (wtvelocity(n) .GT. cut_out_speed - buffer .AND. wtvelocity(&
&       n) .LT. cut_out_speed + buffer) THEN
      x0 = cut_out_speed - buffer
      x1 = cut_out_speed + buffer
      y0 = rated_power
      y1 = 0.
      dy0 = 0.
      dy1 = 0.
      CALL HERMITE_SPLINE_DV(wtvelocity(n), wtvelocityd(:, n), x0, x1, &
&                      y0, dy0, y1, dy1, wtpower(n), wtpowerd(:, n), &
&                      nbdirs)
! If we're above cut-out
    ELSE IF (wtvelocity(n) .GT. cut_out_speed + buffer) THEN
      DO nd=1,nbdirs
        wtpowerd(nd, n) = 0.0_8
      END DO
      wtpower(n) = 0.
    END IF
  END DO
  DO nd=1,nbdirs
    dir_powerd(nd) = SUM(wtpowerd(nd, :))
  END DO
  dir_power = SUM(wtpower)
END SUBROUTINE DIRPOWER_DV

!  Differentiation of porteagel_analyze in forward (tangent) mode:
!   variations   of useful results: wtvelocity
!   with respect to varying inputs: turbinexw turbineyw
! implementation of the Bastankhah and Porte Agel (BPA) wake model for analysis
SUBROUTINE PORTEAGEL_ANALYZE_DV(nturbines, nrotorpoints, nctpoints, &
& turbinexw, turbinexwd, sorted_x_idx, turbineyw, turbineywd, turbinez, &
& rotordiameter, ct, wind_speed, yawdeg, ky, kz, alpha, beta, ti, &
& rotorpointsy, rotorpointsz, z_ref, z_0, shear_exp, &
& wake_combination_method, ti_calculation_method, calc_k_star, &
& wec_factor, print_ti, wake_model_version, interp_type, use_ct_curve, &
& ct_curve_wind_speed, ct_curve_ct, sm_smoothing, wtvelocity, &
& wtvelocityd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
! print *, Ct_local(turbI)
!!  print TIturbs values to a file
!     if (print_ti) then
!         open(unit=2, file="TIturbs_tmp.txt")
!         do, turb=1, nTurbines
!             write(2,*) TIturbs(turb)
!         end do
!         close(2)
!     end if
!print *, "TIturbs: ", TIturbs
!print *, wtVelocity
!! make sure turbine inflow velocity is non-negative
!             if (wtVelocity(turbI) .lt. 0.0_dp) then
!                 wtVelocity(turbI) = 0.0_dp
!             end if
!print *, "fortran"
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines, nrotorpoints, nctpoints
  INTEGER, INTENT(IN) :: wake_combination_method, ti_calculation_method&
& , wake_model_version, interp_type
  LOGICAL, INTENT(IN) :: calc_k_star, print_ti, use_ct_curve
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw, turbineyw, &
& turbinez
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbinexwd, &
& turbineywd
  INTEGER, DIMENSION(nturbines), INTENT(IN) :: sorted_x_idx
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: rotordiameter, yawdeg
  REAL(dp), DIMENSION(nturbines) :: ct
  REAL(dp), INTENT(IN) :: ky, kz, alpha, beta, ti, wind_speed, z_ref, &
& z_0, shear_exp, wec_factor
  REAL(dp), DIMENSION(nrotorpoints), INTENT(IN) :: rotorpointsy, &
& rotorpointsz
  REAL(dp), DIMENSION(nctpoints), INTENT(IN) :: ct_curve_wind_speed, &
& ct_curve_ct
  REAL(dp), INTENT(IN) :: sm_smoothing
! local (General)
  REAL(dp), DIMENSION(nturbines) :: yaw, titurbs, ct_local, ky_local, &
& kz_local
  REAL(dp), DIMENSION(nbdirs, nturbines) :: titurbsd, ct_locald, &
& ky_locald, kz_locald
  REAL(dp) :: x0, deltax0, deltay, theta_c_0, sigmay, sigmaz, &
& wake_offset, k_star
  REAL(dp), DIMENSION(nbdirs) :: x0d, deltax0d, deltayd, theta_c_0d, &
& sigmayd, sigmazd, wake_offsetd, k_stard
  REAL(dp) :: x, deltav, deltaz, sigmay_dp, sigmaz_dp, deltax0_dp, &
& deficit_sum
  REAL(dp), DIMENSION(nbdirs) :: xd, deltavd, sigmay_dpd, sigmaz_dpd&
& , deltax0_dpd, deficit_sumd
  REAL(dp) :: tol, discontinuity_point, ti_area_ratio
  REAL(dp), DIMENSION(nbdirs) :: discontinuity_pointd, ti_area_ratiod
  REAL(dp) :: ti_area_ratio_tmp, ti_dst_tmp, ti_ust_tmp, rpts
  REAL(dp), DIMENSION(nbdirs) :: ti_area_ratio_tmpd, ti_dst_tmpd, &
& ti_ust_tmpd
  REAL(dp) :: localrotorpointy, localrotorpointz, point_velocity, &
& point_z, point_velocity_with_shear
  REAL(dp), DIMENSION(nbdirs) :: point_velocityd, &
& point_velocity_with_sheard
  INTEGER :: u, d, turb, turbi, p
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! model out
  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: wtvelocity
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: wtvelocityd
  INTRINSIC SIN, COS, ATAN, MAX, SQRT, LOG
  INTRINSIC REAL
  INTEGER :: nd
  INTEGER :: nbdirs
! bastankhah and porte agel 2016 define yaw to be positive clockwise, this is reversed
  yaw = -(yawdeg*pi/180.0_dp)
! set tolerance for location checks
  tol = 0.1_dp
! initialize wind turbine velocities to 0.0
  wtvelocity = 0.0_dp
! initialize TI of all turbines to free-stream value
!print *, "start TIturbs: ", TIturbs
  titurbs = ti
! initialize the local wake factors
  ky_local(:) = ky
  kz_local(:) = kz
  ct_local(:) = ct
  wtvelocityd(:, :) = 0.0_8
  titurbsd(:, :) = 0.0_8
  kz_locald(:, :) = 0.0_8
  ct_locald(:, :) = 0.0_8
  ky_locald(:, :) = 0.0_8
!print *, 'wake model version: ', wake_model_version
!print *, "ky_local: ", ky_local
!print *, "kz_local: ", kz_local
!print *, "TIturbs init: ", TIturbs
  DO d=1,nturbines
! get index of downstream turbine
    turbi = sorted_x_idx(d) + 1
    DO p=1,nrotorpoints
! initialize the TI_area_ratio to 0.0 for each turbine
      ti_area_ratio = 0.0_dp
! initialize deficit summation term to zero
      deficit_sum = 0.0_dp
! scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
      localrotorpointy = rotorpointsy(p)*0.5_dp*rotordiameter(turbi)
      localrotorpointz = rotorpointsz(p)*0.5_dp*rotordiameter(turbi)
      deficit_sumd(:) = 0.0_8
      ti_area_ratiod(:) = 0.0_8
!             !print *, "rotorDiameter after local rotor points", rotorDiameter
!             !print *, "local rotor points Y,Z: ", LocalRotorPointY, LocalRotorPointZ
! at turbineX-locations
      DO u=1,nturbines
! get index of upstream turbine
        turb = sorted_x_idx(u) + 1
! skip this loop if turb = turbI (turbines impact on itself)
        IF (turb .NE. turbi) THEN
          DO nd=1,nbdirs
! downstream distance between upstream turbine and point
            xd(nd) = turbinexwd(nd, turbi) - turbinexwd(nd, turb)
          END DO
          x = turbinexw(turbi) - turbinexw(turb) + localrotorpointy*SIN(&
&           yaw(turbi))
! set this iterations velocity deficit to 0
          deltav = 0.0_dp
! check turbine relative locations
          IF (x .GT. 0.0_dp + tol) THEN
!print *, "rotorDiameter before x0 ", rotorDiameter
! determine the onset location of far wake
            CALL X0_FUNC_DV(rotordiameter(turb), yaw(turb), ct_local(&
&                     turb), ct_locald(:, turb), alpha, titurbs(turb), &
&                     titurbsd(:, turb), beta, x0, x0d, nbdirs)
            DO nd=1,nbdirs
!                     call x0_func(rotorDiameter(turb), yaw(turb), Ct(turb), alpha, &
!                                 & TI, beta, x0)
! downstream distance from far wake onset to downstream turbine
              deltax0d(nd) = xd(nd) - x0d(nd)
            END DO
            deltax0 = x - x0
! calculate wake spreading parameter at each turbine if desired
            IF (calc_k_star .EQV. .true.) THEN
              CALL K_STAR_FUNC_DV(titurbs(turb), titurbsd(:, turb), &
&                           k_star, k_stard, nbdirs)
              DO nd=1,nbdirs
                ky_locald(nd, turb) = k_stard(nd)
                kz_locald(nd, turb) = k_stard(nd)
              END DO
              ky_local(turb) = k_star
              kz_local(turb) = k_star
            END IF
!print *, "ky_local ", ky_local
!print *, "deltax0 ", deltax0
!print *, "turbineZ ", turbineZ
!print *, "rotorDiameter after x0 ", rotorDiameter
!print *, "Ct ", Ct
!print *, "yaw ", yaw
! determine the initial wake angle at the onset of far wake
            CALL THETA_C_0_FUNC_DV(yaw(turb), ct_local(turb), ct_locald(&
&                            :, turb), theta_c_0, theta_c_0d, nbdirs)
!print *, "theta_c_0 ", theta_c_0
! horizontal spread
            CALL SIGMAY_FUNC_DV(ky_local(turb), ky_locald(:, turb), &
&                         deltax0, deltax0d, rotordiameter(turb), yaw(&
&                         turb), sigmay, sigmayd, nbdirs)
!print *, "sigmay ", sigmay
!print *, "rotorDiameter after sigmay", rotorDiameter
! vertical spread
            CALL SIGMAZ_FUNC_DV(kz_local(turb), kz_locald(:, turb), &
&                         deltax0, deltax0d, rotordiameter(turb), sigmaz&
&                         , sigmazd, nbdirs)
!print *, "sigmaz ", sigmaz
!print *, "rotorDiameter after sigmaz ", rotorDiameter
! horizontal cross-wind wake displacement from hub
            CALL WAKE_OFFSET_FUNC_DV(rotordiameter(turb), theta_c_0, &
&                              theta_c_0d, x0, x0d, yaw(turb), ky_local(&
&                              turb), ky_locald(:, turb), kz_local(turb)&
&                              , kz_locald(:, turb), ct_local(turb), &
&                              ct_locald(:, turb), sigmay, sigmayd, &
&                              sigmaz, sigmazd, wake_offset, &
&                              wake_offsetd, nbdirs)
            DO nd=1,nbdirs
!print *, "wake_offset ", wake_offset
! cross wind distance from downstream point location to wake center
              deltayd(nd) = turbineywd(nd, turbi) - turbineywd(nd, turb)&
&               - wake_offsetd(nd)
            END DO
            deltay = localrotorpointy*COS(yaw(turbi)) + turbineyw(turbi)&
&             - (turbineyw(turb)+wake_offset)
! cross wind distance from hub height to height of point of interest
            deltaz = localrotorpointz + turbinez(turbi) - turbinez(turb)
!print *, "dx, dy, dz: ", x, deltay, deltaz
!print *, "local y,z : ", LocalRotorPointY, LocalRotorPointZ, turb, turbI, p
!print *, deltaz, deltay
! far wake region
! find the final point where the original model is undefined
            CALL DISCONTINUITY_POINT_FUNC_DV(x0, x0d, rotordiameter(turb&
&                                      ), ky_local(turb), ky_locald(:, &
&                                      turb), kz_local(turb), kz_locald(&
&                                      :, turb), yaw(turb), ct_local(&
&                                      turb), ct_locald(:, turb), &
&                                      discontinuity_point, &
&                                      discontinuity_pointd, nbdirs)
!print *, "discontinuity point is: ", discontinuity_point
            IF (x .GT. discontinuity_point) THEN
!print *, x
! velocity difference in the wake
              CALL DELTAV_FUNC_DV(deltay, deltayd, deltaz, ct_local(turb&
&                           ), ct_locald(:, turb), yaw(turb), sigmay, &
&                           sigmayd, sigmaz, sigmazd, rotordiameter(turb&
&                           ), wake_model_version, kz_local(turb), &
&                           kz_locald(:, turb), x, xd, wec_factor, &
&                           deltav, deltavd, nbdirs)
!print *, "rotorDiameter after far deltav ", rotorDiameter
! near wake region (linearized)
            ELSE
              DO nd=1,nbdirs
! determine distance from discontinuity point to far wake onset
                deltax0_dpd(nd) = discontinuity_pointd(nd) - x0d(nd)
              END DO
              deltax0_dp = discontinuity_point - x0
! horizontal spread at far wake onset
              CALL SIGMAY_FUNC_DV(ky_local(turb), ky_locald(:, turb), &
&                           deltax0_dp, deltax0_dpd, rotordiameter(turb)&
&                           , yaw(turb), sigmay_dp, sigmay_dpd, nbdirs)
! vertical spread at far wake onset
              CALL SIGMAZ_FUNC_DV(kz_local(turb), kz_locald(:, turb), &
&                           deltax0_dp, deltax0_dpd, rotordiameter(turb)&
&                           , sigmaz_dp, sigmaz_dpd, nbdirs)
!  print *, "inputs in parent: ", deltay, deltaz, Ct(turb), yaw(turb), sigmay_dp, sigmaz_dp, &
!                                          & rotorDiameter(turb), x, discontinuity_point, sigmay_dp, sigmaz_dp, &
!                                          & wake_model_version, kz_local, x0, &
!                                          & wec_factor
! velocity deficit in the nearwake (linear model)
              CALL DELTAV_NEAR_WAKE_LIN_FUNC_DV(deltay, deltayd, deltaz&
&                                         , ct_local(turb), ct_locald(:&
&                                         , turb), yaw(turb), sigmay_dp&
&                                         , sigmay_dpd, sigmaz_dp, &
&                                         sigmaz_dpd, rotordiameter(turb&
&                                         ), x, xd, discontinuity_point&
&                                         , discontinuity_pointd, &
&                                         sigmay_dp, sigmay_dpd, &
&                                         sigmaz_dp, sigmaz_dpd, &
&                                         wake_model_version, kz_local(&
&                                         turb), kz_locald(:, turb), x0&
&                                         , x0d, wec_factor, deltav, &
&                                         deltavd, nbdirs)
!print *, "rotorDiameter after deltav near ", rotorDiameter
            END IF
! combine deficits according to selected method wake combination method
            CALL WAKE_COMBINATION_FUNC_DV(wind_speed, wtvelocity(turb), &
&                                   wtvelocityd(:, turb), deltav, &
&                                   deltavd, wake_combination_method, &
&                                   deficit_sum, deficit_sumd, nbdirs)
            IF (x .GT. 0.0_dp .AND. ti_calculation_method .GT. 0) THEN
              DO nd=1,nbdirs
!print *, "turbI, turb: ", turbI, turb
! calculate TI value at each turbine
!                         print *, "turb, turbI: ", turb, turbI
! save ti_area_ratio and ti_dst to new memory locations to avoid
! aliasing during differentiation
                ti_area_ratio_tmpd(nd) = ti_area_ratiod(nd)
                ti_dst_tmpd(nd) = titurbsd(nd, turbi)
                ti_ust_tmpd(nd) = titurbsd(nd, turb)
              END DO
              ti_area_ratio_tmp = ti_area_ratio
              ti_dst_tmp = titurbs(turbi)
              ti_ust_tmp = titurbs(turb)
              CALL ADDED_TI_FUNC_DV(ti, ct_local(turb), ct_locald(:, &
&                             turb), x, xd, ky_local(turb), ky_locald(:&
&                             , turb), rotordiameter(turb), &
&                             rotordiameter(turbi), deltay, deltayd, &
&                             turbinez(turb), turbinez(turbi), &
&                             sm_smoothing, ti_ust_tmp, ti_ust_tmpd, &
&                             ti_calculation_method, ti_area_ratio_tmp, &
&                             ti_area_ratio_tmpd, ti_dst_tmp, &
&                             ti_dst_tmpd, ti_area_ratio, ti_area_ratiod&
&                             , titurbs(turbi), titurbsd(:, turbi), &
&                             nbdirs)
!print *, "rotorDiameter after TI calcs", rotorDiameter
            END IF
          END IF
        END IF
      END DO
      DO nd=1,nbdirs
!                     !print *, "deficit_sum, turbI, p, turb: ", deficit_sum, turbI, p, turb
! print *, deficit_sum
! find velocity at point p due to the wake of turbine turb
        point_velocityd(nd) = -deficit_sumd(nd)
      END DO
      point_velocity = wind_speed - deficit_sum
!print *, "point velocity, deficit_sum, turbI, p: ", point_velocity, deficit_sum, turbI, p
! put sample point height in global reference frame
      point_z = localrotorpointz + turbinez(turbi)
!print *, "point_z, turbI, p: ", point_z, turbI, p
! adjust sample point velocity for shear
      CALL WIND_SHEAR_FUNC_DV(point_z, point_velocity, point_velocityd, &
&                       z_ref, z_0, shear_exp, point_velocity_with_shear&
&                       , point_velocity_with_sheard, nbdirs)
      DO nd=1,nbdirs
!print *, "v, vs, x, turb, turbI, p: ", point_velocity, point_velocity_with_shear, x, turb, turbI, p
! add sample point velocity to turbine velocity to be averaged later
        wtvelocityd(nd, turbi) = wtvelocityd(nd, turbi) + &
&         point_velocity_with_sheard(nd)
      END DO
      wtvelocity(turbi) = wtvelocity(turbi) + point_velocity_with_shear
    END DO
! final velocity calculation for turbine turbI (average equally across all points)
    rpts = REAL(nrotorpoints, dp)
    DO nd=1,nbdirs
!         print *, rpts, nRotorPoints, wtVelocity(turbI), wtVelocity(turbI)/rpts, wtVelocity(turbI)/nRotorPoints
!         STOP 1
      wtvelocityd(nd, turbi) = wtvelocityd(nd, turbi)/rpts
    END DO
    wtvelocity(turbi) = wtvelocity(turbi)/rpts
!         print *, wtVelocity(turbI)
    IF (use_ct_curve) THEN
! print *, "wtVelocity(turbI): ", wtVelocity(turbI)
      CALL INTERPOLATION_DV(nctpoints, interp_type, ct_curve_wind_speed&
&                     , ct_curve_ct, wtvelocity(turbi), wtvelocityd(:, &
&                     turbi), ct_local(turbi), ct_locald(:, turbi), &
&                     0.0_dp, 0.0_dp, .false., nbdirs)
! print *, "Ct_local(turbI): ", Ct_local(turbI)
    END IF
  END DO
END SUBROUTINE PORTEAGEL_ANALYZE_DV

!  Differentiation of x0_func in forward (tangent) mode:
!   variations   of useful results: x0
!   with respect to varying inputs: ti ct
! calculates the onset of far-wake conditions
SUBROUTINE X0_FUNC_DV(rotor_diameter, yaw, ct, ctd, alpha, ti, tid, beta&
& , x0, x0d, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: rotor_diameter, yaw, ct, alpha, ti, beta
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ctd, tid
! out
  REAL(dp), INTENT(OUT) :: x0
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: x0d
  INTRINSIC COS, SQRT
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  REAL(dp) :: result3
  REAL(dp), DIMENSION(nbdirs) :: result3d
  INTEGER :: nd
  INTEGER :: nbdirs
  result1 = SQRT(1.0_dp - ct)
  result2 = SQRT(2.0_dp)
  result3 = SQRT(1.0_dp - ct)
  DO nd=1,nbdirs
! determine the onset location of far wake
    result1d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
    result3d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
    x0d(nd) = rotor_diameter*(COS(yaw)*result1d(nd)*result2*(alpha*ti+&
&     beta*(1.0_dp-result3))-COS(yaw)*(1.0_dp+result1)*result2*(alpha*&
&     tid(nd)-beta*result3d(nd)))/(result2**2*(alpha*ti+beta*(1.0_dp-&
&     result3))**2)
  END DO
  x0 = rotor_diameter*(COS(yaw)*(1.0_dp+result1)/(result2*(alpha*ti+beta&
&   *(1.0_dp-result3))))
END SUBROUTINE X0_FUNC_DV

!  Differentiation of theta_c_0_func in forward (tangent) mode:
!   variations   of useful results: theta_c_0
!   with respect to varying inputs: ct
! calculates the wake angle at the onset of far wake conditions
SUBROUTINE THETA_C_0_FUNC_DV(yaw, ct, ctd, theta_c_0, theta_c_0d, nbdirs&
&)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: yaw, ct
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ctd
! out
  REAL(dp), INTENT(OUT) :: theta_c_0
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: theta_c_0d
  INTRINSIC COS, SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  INTEGER :: nd
  INTEGER :: nbdirs
  arg1 = 1.0_dp - ct*COS(yaw)
  DO nd=1,nbdirs
! determine the initial wake angle at the onset of far wake
    arg1d(nd) = -(COS(yaw)*ctd(nd))
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0_8
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    theta_c_0d(nd) = (-(0.3_dp*yaw*result1d(nd)))/COS(yaw)
  END DO
  result1 = SQRT(arg1)
  theta_c_0 = 0.3_dp*yaw*(1.0_dp-result1)/COS(yaw)
END SUBROUTINE THETA_C_0_FUNC_DV

!  Differentiation of sigmay_func in forward (tangent) mode:
!   variations   of useful results: sigmay
!   with respect to varying inputs: ky deltax0
! calculates the horizontal spread of the wake at a given distance from the onset of far
! wake condition
SUBROUTINE SIGMAY_FUNC_DV(ky, kyd, deltax0, deltax0d, rotor_diameter, &
& yaw, sigmay, sigmayd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: ky, deltax0, rotor_diameter, yaw
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: kyd, deltax0d
! out
  REAL(dp), INTENT(OUT) :: sigmay
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: sigmayd
  INTRINSIC COS, SQRT
  REAL(dp) :: result1
  INTEGER :: nd
  INTEGER :: nbdirs
! horizontal spread
  result1 = SQRT(8.0_dp)
  DO nd=1,nbdirs
    sigmayd(nd) = kyd(nd)*deltax0 + ky*deltax0d(nd)
  END DO
  sigmay = rotor_diameter*(ky*deltax0/rotor_diameter+COS(yaw)/result1)
END SUBROUTINE SIGMAY_FUNC_DV

!  Differentiation of sigmaz_func in forward (tangent) mode:
!   variations   of useful results: sigmaz
!   with respect to varying inputs: kz deltax0
! calculates the vertical spread of the wake at a given distance from the onset of far
! wake condition
SUBROUTINE SIGMAZ_FUNC_DV(kz, kzd, deltax0, deltax0d, rotor_diameter, &
& sigmaz, sigmazd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: kz, deltax0, rotor_diameter
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: kzd, deltax0d
! out
  REAL(dp), INTENT(OUT) :: sigmaz
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: sigmazd
! load necessary intrinsic functions
  INTRINSIC SQRT
  REAL(dp) :: result1
  INTEGER :: nd
  INTEGER :: nbdirs
! vertical spread
  result1 = SQRT(8.0_dp)
  DO nd=1,nbdirs
    sigmazd(nd) = kzd(nd)*deltax0 + kz*deltax0d(nd)
  END DO
  sigmaz = rotor_diameter*(kz*deltax0/rotor_diameter+1.0_dp/result1)
END SUBROUTINE SIGMAZ_FUNC_DV

!  Differentiation of wake_offset_func in forward (tangent) mode:
!   variations   of useful results: wake_offset
!   with respect to varying inputs: theta_c_0 sigmay sigmaz ky
!                kz x0 ct
! calculates the horizontal distance from the wake center to the hub of the turbine making
! the wake
SUBROUTINE WAKE_OFFSET_FUNC_DV(rotor_diameter, theta_c_0, theta_c_0d, x0&
& , x0d, yaw, ky, kyd, kz, kzd, ct, ctd, sigmay, sigmayd, sigmaz, &
& sigmazd, wake_offset, wake_offsetd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: rotor_diameter, theta_c_0, x0, yaw, ky, kz, ct&
& , sigmay
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: theta_c_0d, x0d, kyd, &
& kzd, ctd, sigmayd
  REAL(dp), INTENT(IN) :: sigmaz
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: sigmazd
! out
  REAL(dp), INTENT(OUT) :: wake_offset
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: wake_offsetd
  INTRINSIC COS, SQRT, LOG
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  REAL(dp), DIMENSION(nbdirs) :: result2d
  REAL(dp) :: result3
  REAL(dp), DIMENSION(nbdirs) :: result3d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  REAL(dp) :: result4
  REAL(dp), DIMENSION(nbdirs) :: result4d
  REAL(dp) :: result5
  REAL(dp), DIMENSION(nbdirs) :: result5d
  REAL(dp) :: result6
  REAL(dp), DIMENSION(nbdirs) :: result6d
  REAL(dp) :: arg3
  REAL(dp), DIMENSION(nbdirs) :: arg3d
  REAL(dp) :: result7
  REAL(dp), DIMENSION(nbdirs) :: result7d
  REAL(dp) :: result8
  REAL(dp), DIMENSION(nbdirs) :: result8d
  REAL(dp) :: arg4
  REAL(dp), DIMENSION(nbdirs) :: arg4d
  INTEGER :: nd
  INTEGER :: nbdirs
  arg1 = COS(yaw)/(ky*kz*ct)
  result1 = SQRT(arg1)
  result2 = SQRT(1.0_dp - ct)
  result3 = SQRT(ct)
  arg2 = 8.0_dp*sigmay*sigmaz/(COS(yaw)*rotor_diameter**2)
  result4 = SQRT(arg2)
  result5 = SQRT(ct)
  result6 = SQRT(ct)
  arg3 = 8.0_dp*sigmay*sigmaz/(COS(yaw)*rotor_diameter**2)
  result7 = SQRT(arg3)
  result8 = SQRT(ct)
  arg4 = (1.6_dp+result3)*(1.6_dp*result4-result5)/((1.6_dp-result6)*(&
&   1.6_dp*result7+result8))
  DO nd=1,nbdirs
! horizontal cross-wind wake displacement from hub
    arg1d(nd) = -(COS(yaw)*((kyd(nd)*kz+ky*kzd(nd))*ct+ky*kz*ctd(nd))/(&
&     ky*kz*ct)**2)
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0_8
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    result2d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
    IF (ct .EQ. 0.0) THEN
      result3d(nd) = 0.0_8
    ELSE
      result3d(nd) = ctd(nd)/(2.0*SQRT(ct))
    END IF
    arg2d(nd) = 8.0_dp*(sigmayd(nd)*sigmaz+sigmay*sigmazd(nd))/(COS(yaw)&
&     *rotor_diameter**2)
    IF (arg2 .EQ. 0.0) THEN
      result4d(nd) = 0.0_8
    ELSE
      result4d(nd) = arg2d(nd)/(2.0*SQRT(arg2))
    END IF
    IF (ct .EQ. 0.0) THEN
      result5d(nd) = 0.0_8
    ELSE
      result5d(nd) = ctd(nd)/(2.0*SQRT(ct))
    END IF
    IF (ct .EQ. 0.0) THEN
      result6d(nd) = 0.0_8
    ELSE
      result6d(nd) = ctd(nd)/(2.0*SQRT(ct))
    END IF
    arg3d(nd) = 8.0_dp*(sigmayd(nd)*sigmaz+sigmay*sigmazd(nd))/(COS(yaw)&
&     *rotor_diameter**2)
    IF (arg3 .EQ. 0.0) THEN
      result7d(nd) = 0.0_8
    ELSE
      result7d(nd) = arg3d(nd)/(2.0*SQRT(arg3))
    END IF
    IF (ct .EQ. 0.0) THEN
      result8d(nd) = 0.0_8
    ELSE
      result8d(nd) = ctd(nd)/(2.0*SQRT(ct))
    END IF
    arg4d(nd) = ((result3d(nd)*(1.6_dp*result4-result5)+(1.6_dp+result3)&
&     *(1.6_dp*result4d(nd)-result5d(nd)))*(1.6_dp-result6)*(1.6_dp*&
&     result7+result8)-(1.6_dp+result3)*(1.6_dp*result4-result5)*((&
&     1.6_dp-result6)*(1.6_dp*result7d(nd)+result8d(nd))-result6d(nd)*(&
&     1.6_dp*result7+result8)))/((1.6_dp-result6)*(1.6_dp*result7+&
&     result8))**2
    wake_offsetd(nd) = rotor_diameter*((theta_c_0d(nd)*x0+theta_c_0*x0d(&
&     nd))/rotor_diameter+((theta_c_0d(nd)*result1/14.7_dp+theta_c_0*&
&     result1d(nd)/14.7_dp)*(2.9_dp+1.3_dp*result2-ct)+theta_c_0*result1&
&     *(1.3_dp*result2d(nd)-ctd(nd))/14.7_dp)*LOG(arg4)+theta_c_0*&
&     result1*(2.9_dp+1.3_dp*result2-ct)*arg4d(nd)/(14.7_dp*arg4))
  END DO
  wake_offset = rotor_diameter*(theta_c_0*x0/rotor_diameter+theta_c_0/&
&   14.7_dp*result1*(2.9_dp+1.3_dp*result2-ct)*LOG(arg4))
END SUBROUTINE WAKE_OFFSET_FUNC_DV

!  Differentiation of deltav_func in forward (tangent) mode:
!   variations   of useful results: deltav
!   with respect to varying inputs: k sigmay sigmaz deltax deltay
!                ct
! calculates the velocity difference between hub velocity and free stream for a given wake
! for use in the far wake region
SUBROUTINE DELTAV_FUNC_DV(deltay, deltayd, deltaz, ct, ctd, yaw, sigmay&
& , sigmayd, sigmaz, sigmazd, rotor_diameter_ust, version, k, kd, deltax&
& , deltaxd, wec_factor, deltav, deltavd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!print *, "rotor_diameter in deltav exit", rotor_diameter_ust
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: deltay, deltaz, ct, yaw, sigmay
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: deltayd, ctd, sigmayd
  REAL(dp), INTENT(IN) :: sigmaz, rotor_diameter_ust, wec_factor
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: sigmazd
! only for 2014 version
  REAL(dp), INTENT(IN) :: k, deltax
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: kd, deltaxd
  INTEGER, INTENT(IN) :: version
! local
! only for 2014 version
  REAL(dp) :: beta_2014, epsilon_2014
  REAL(dp), DIMENSION(nbdirs) :: beta_2014d, epsilon_2014d
! out
  REAL(dp), INTENT(OUT) :: deltav
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: deltavd
! load intrinsic functions
  INTRINSIC COS, SQRT, EXP
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  REAL(dp), DIMENSION(nbdirs) :: result2d
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  REAL(dp) :: arg3
  REAL(dp), DIMENSION(nbdirs) :: arg3d
  INTEGER :: nd
  INTEGER :: nbdirs
!print *, "rotor_diameter in deltav entry", rotor_diameter_ust
!     print *, 'wake model version in deltav: ', version
  IF (version .EQ. 2014) THEN
    result1 = SQRT(1.0_dp - ct)
    result2 = SQRT(1.0_dp - ct)
    beta_2014 = 0.5_dp*(1.0_dp+result1)/result2
    DO nd=1,nbdirs
!print *, "in 2014 version"
      result1d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
      result2d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
      beta_2014d(nd) = (0.5_dp*result1d(nd)*result2-0.5_dp*(1.0_dp+&
&       result1)*result2d(nd))/result2**2
      IF (beta_2014 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = beta_2014d(nd)/(2.0*SQRT(beta_2014))
      END IF
      epsilon_2014d(nd) = 0.2_dp*result1d(nd)
    END DO
    result1 = SQRT(beta_2014)
    epsilon_2014 = 0.2_dp*result1
    arg1 = 1.0_dp - ct/(8.0_dp*(k*deltax/rotor_diameter_ust+epsilon_2014&
&     )**2)
    result1 = SQRT(arg1)
    arg2 = (-(1.0_dp/(2.0_dp*(k*deltax/rotor_diameter_ust+epsilon_2014)&
&     **2)))*((deltaz/(wec_factor*rotor_diameter_ust))**2+(deltay/(&
&     wec_factor*rotor_diameter_ust))**2)
    DO nd=1,nbdirs
! print *, "beta = ", beta_2014, "epsilon = ", epsilon_2014
! print *, "k, deltax: ", k, deltax
! print *, "term: ", Ct                                                   &
!                            / (8.0_dp * (k*deltax/rotor_diameter_ust+epsilon_2014)**2)
      arg1d(nd) = -((ctd(nd)*8.0_dp*(k*deltax/rotor_diameter_ust+&
&       epsilon_2014)**2-ct*8.0_dp*2*(k*deltax/rotor_diameter_ust+&
&       epsilon_2014)*((kd(nd)*deltax+k*deltaxd(nd))/rotor_diameter_ust+&
&       epsilon_2014d(nd)))/(8.0_dp*(k*deltax/rotor_diameter_ust+&
&       epsilon_2014)**2)**2)
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      arg2d(nd) = 2*((kd(nd)*deltax+k*deltaxd(nd))/rotor_diameter_ust+&
&       epsilon_2014d(nd))*((deltaz/(wec_factor*rotor_diameter_ust))**2+&
&       (deltay/(wec_factor*rotor_diameter_ust))**2)/(2.0_dp*(k*deltax/&
&       rotor_diameter_ust+epsilon_2014)**3) - 2*deltay*deltayd(nd)/(&
&       2.0_dp*(k*deltax/rotor_diameter_ust+epsilon_2014)**2*wec_factor&
&       **2*rotor_diameter_ust**2)
      deltavd(nd) = (1.0_dp-result1)*arg2d(nd)*EXP(arg2) - result1d(nd)*&
&       EXP(arg2)
    END DO
    deltav = (1.0_dp-result1)*EXP(arg2)
! print *, "deltav 2014 = ", deltav
  ELSE IF (version .EQ. 2016) THEN
    arg1 = 1.0_dp - ct*COS(yaw)/(8.0_dp*sigmay*sigmaz/rotor_diameter_ust&
&     **2)
    result1 = SQRT(arg1)
    arg2 = -(0.5_dp*(deltay/(wec_factor*sigmay))**2)
    arg3 = -(0.5_dp*(deltaz/(wec_factor*sigmaz))**2)
    DO nd=1,nbdirs
! velocity difference in the wake at each sample point
      arg1d(nd) = -((COS(yaw)*ctd(nd)*8.0_dp*sigmay*sigmaz/&
&       rotor_diameter_ust**2-ct*COS(yaw)*8.0_dp*(sigmayd(nd)*sigmaz+&
&       sigmay*sigmazd(nd))/rotor_diameter_ust**2)/(8.0_dp*sigmay*sigmaz&
&       /rotor_diameter_ust**2)**2)
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      arg2d(nd) = -(0.5_dp*2*deltay*(deltayd(nd)*wec_factor*sigmay-&
&       deltay*wec_factor*sigmayd(nd))/(wec_factor**3*sigmay**3))
      arg3d(nd) = 0.5_dp*2*deltaz**2*sigmazd(nd)/(wec_factor**2*sigmaz**&
&       3)
      deltavd(nd) = ((1.0_dp-result1)*arg2d(nd)*EXP(arg2)-result1d(nd)*&
&       EXP(arg2))*EXP(arg3) + (1.0_dp-result1)*EXP(arg2)*arg3d(nd)*EXP(&
&       arg3)
    END DO
    deltav = (1.0_dp-result1)*EXP(arg2)*EXP(arg3)
  ELSE
    PRINT*, 'Invalid Bastankhah and Porte Agel model version. Must be 20&
&14 or 2016. ', version, ' was given.'
    STOP
  END IF
END SUBROUTINE DELTAV_FUNC_DV

!  Differentiation of deltav_near_wake_lin_func in forward (tangent) mode:
!   variations   of useful results: deltav
!   with respect to varying inputs: k discontinuity_point x sigmay
!                sigmaz deltax0_dp deltay sigmay0 ct sigmaz0
! calculates the velocity difference between hub velocity and free stream for a given wake
! for use in the near wake region only
SUBROUTINE DELTAV_NEAR_WAKE_LIN_FUNC_DV(deltay, deltayd, deltaz, ct, ctd&
& , yaw, sigmay, sigmayd, sigmaz, sigmazd, rotor_diameter_ust, x, xd, &
& discontinuity_point, discontinuity_pointd, sigmay0, sigmay0d, sigmaz0&
& , sigmaz0d, version, k, kd, deltax0_dp, deltax0_dpd, wec_factor, &
& deltav, deltavd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: deltay, deltaz, ct, yaw, sigmay
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: deltayd, ctd, sigmayd
  REAL(dp), INTENT(IN) :: sigmaz, rotor_diameter_ust, wec_factor
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: sigmazd
  REAL(dp), INTENT(IN) :: x, discontinuity_point, sigmay0, sigmaz0
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: xd, discontinuity_pointd&
& , sigmay0d, sigmaz0d
! only for 2014 version
  REAL(dp), INTENT(IN) :: k, deltax0_dp
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: kd, deltax0_dpd
  INTEGER, INTENT(IN) :: version
! local
  REAL(dp) :: deltav0m, deltavr
  REAL(dp), DIMENSION(nbdirs) :: deltav0md, deltavrd
! only for 2014 version
  REAL(dp) :: beta_2014, epsilon_2014
  REAL(dp), DIMENSION(nbdirs) :: beta_2014d, epsilon_2014d
! out
  REAL(dp), INTENT(OUT) :: deltav
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: deltavd
! load intrinsic functions
  INTRINSIC COS, SQRT, EXP
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  REAL(dp), DIMENSION(nbdirs) :: result2d
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  INTEGER :: nd
  INTEGER :: nbdirs
!  print *, 'wake model version in deltav near wake: ', version
!     print *, "inputs: ", deltay, deltaz, Ct, yaw,  &
!                                  & sigmay, sigmaz, rotor_diameter_ust, x, &
!                                  & discontinuity_point, sigmay0, sigmaz0, version, k, &
!                                  & deltax0_dp, wec_factor
  IF (version .EQ. 2014) THEN
    IF (yaw .GT. 0.0_dp) THEN
      PRINT*, 'model version 2014 may only be used when yaw=0'
      STOP
    ELSE
      result1 = SQRT(1.0_dp - ct)
      result2 = SQRT(1.0_dp - ct)
      beta_2014 = 0.5_dp*(1.0_dp+result1)/result2
      DO nd=1,nbdirs
        result1d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
        result2d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
        beta_2014d(nd) = (0.5_dp*result1d(nd)*result2-0.5_dp*(1.0_dp+&
&         result1)*result2d(nd))/result2**2
        IF (beta_2014 .EQ. 0.0) THEN
          result1d(nd) = 0.0_8
        ELSE
          result1d(nd) = beta_2014d(nd)/(2.0*SQRT(beta_2014))
        END IF
        epsilon_2014d(nd) = 0.2_dp*result1d(nd)
      END DO
      result1 = SQRT(beta_2014)
      epsilon_2014 = 0.2_dp*result1
      arg1 = 1.0_dp - ct/(8.0_dp*(k*deltax0_dp/rotor_diameter_ust+&
&       epsilon_2014)**2)
      DO nd=1,nbdirs
! magnitude term of gaussian at x0
        arg1d(nd) = -((ctd(nd)*8.0_dp*(k*deltax0_dp/rotor_diameter_ust+&
&         epsilon_2014)**2-ct*8.0_dp*2*(k*deltax0_dp/rotor_diameter_ust+&
&         epsilon_2014)*((kd(nd)*deltax0_dp+k*deltax0_dpd(nd))/&
&         rotor_diameter_ust+epsilon_2014d(nd)))/(8.0_dp*(k*deltax0_dp/&
&         rotor_diameter_ust+epsilon_2014)**2)**2)
        IF (arg1 .EQ. 0.0) THEN
          result1d(nd) = 0.0_8
        ELSE
          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
        END IF
        deltav0md(nd) = -result1d(nd)
! initialize the gaussian magnitude term at the rotor for the linear interpolation
        deltavrd(nd) = deltav0md(nd)
! linearized gaussian magnitude term for near wake
        arg1d(nd) = 2*((kd(nd)*deltax0_dp+k*deltax0_dpd(nd))/&
&         rotor_diameter_ust+epsilon_2014d(nd))*((deltaz/(wec_factor*&
&         rotor_diameter_ust))**2+(deltay/(wec_factor*rotor_diameter_ust&
&         ))**2)/(2.0_dp*(k*deltax0_dp/rotor_diameter_ust+epsilon_2014)&
&         **3) - 2*deltay*deltayd(nd)/(2.0_dp*(k*deltax0_dp/&
&         rotor_diameter_ust+epsilon_2014)**2*wec_factor**2*&
&         rotor_diameter_ust**2)
      END DO
      result1 = SQRT(arg1)
      deltav0m = 1.0_dp - result1
      deltavr = deltav0m
      arg1 = (-(1.0_dp/(2.0_dp*(k*deltax0_dp/rotor_diameter_ust+&
&       epsilon_2014)**2)))*((deltaz/(wec_factor*rotor_diameter_ust))**2&
&       +(deltay/(wec_factor*rotor_diameter_ust))**2)
      DO nd=1,nbdirs
        deltavd(nd) = (((deltav0md(nd)-deltavrd(nd))*discontinuity_point&
&         -(deltav0m-deltavr)*discontinuity_pointd(nd))*x/&
&         discontinuity_point**2+(deltav0m-deltavr)*xd(nd)/&
&         discontinuity_point+deltavrd(nd))*EXP(arg1) + ((deltav0m-&
&         deltavr)/discontinuity_point*x+deltavr)*arg1d(nd)*EXP(arg1)
      END DO
      deltav = ((deltav0m-deltavr)/discontinuity_point*x+deltavr)*EXP(&
&       arg1)
    END IF
  ELSE IF (version .EQ. 2016) THEN
    arg1 = 1.0_dp - ct*COS(yaw)/(8.0_dp*sigmay0*sigmaz0/&
&     rotor_diameter_ust**2)
    DO nd=1,nbdirs
! magnitude term of gaussian at x0
      arg1d(nd) = -((COS(yaw)*ctd(nd)*8.0_dp*sigmay0*sigmaz0/&
&       rotor_diameter_ust**2-ct*COS(yaw)*8.0_dp*(sigmay0d(nd)*sigmaz0+&
&       sigmay0*sigmaz0d(nd))/rotor_diameter_ust**2)/(8.0_dp*sigmay0*&
&       sigmaz0/rotor_diameter_ust**2)**2)
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      deltav0md(nd) = -result1d(nd)
! initialize the gaussian magnitude term at the rotor for the linear interpolation
      deltavrd(nd) = deltav0md(nd)
! linearized gaussian magnitude term for near wake
      arg1d(nd) = -(0.5_dp*2*deltay*(deltayd(nd)*wec_factor*sigmay-&
&       deltay*wec_factor*sigmayd(nd))/(wec_factor**3*sigmay**3))
      arg2d(nd) = 0.5_dp*2*deltaz**2*sigmazd(nd)/(wec_factor**2*sigmaz**&
&       3)
    END DO
    result1 = SQRT(arg1)
    deltav0m = 1.0_dp - result1
    deltavr = deltav0m
    arg1 = -(0.5_dp*(deltay/(wec_factor*sigmay))**2)
    arg2 = -(0.5_dp*(deltaz/(wec_factor*sigmaz))**2)
    DO nd=1,nbdirs
      deltavd(nd) = ((((deltav0md(nd)-deltavrd(nd))*discontinuity_point-&
&       (deltav0m-deltavr)*discontinuity_pointd(nd))*x/&
&       discontinuity_point**2+(deltav0m-deltavr)*xd(nd)/&
&       discontinuity_point+deltavrd(nd))*EXP(arg1)+((deltav0m-deltavr)/&
&       discontinuity_point*x+deltavr)*arg1d(nd)*EXP(arg1))*EXP(arg2) + &
&       ((deltav0m-deltavr)/discontinuity_point*x+deltavr)*EXP(arg1)*&
&       arg2d(nd)*EXP(arg2)
    END DO
    deltav = ((deltav0m-deltavr)/discontinuity_point*x+deltavr)*EXP(arg1&
&     )*EXP(arg2)
  ELSE
    PRINT*, 'Invalid Bastankhah and Porte Agel model version. Must be 20&
&14 or 2016. ', version, ' was given.'
    STOP
  END IF
END SUBROUTINE DELTAV_NEAR_WAKE_LIN_FUNC_DV

!  Differentiation of wake_combination_func in forward (tangent) mode:
!   variations   of useful results: deficit_sum
!   with respect to varying inputs: turb_inflow deficit_sum deltav
! combines wakes using various methods
SUBROUTINE WAKE_COMBINATION_FUNC_DV(wind_speed, turb_inflow, &
& turb_inflowd, deltav, deltavd, wake_combination_method, deficit_sum, &
& deficit_sumd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: wind_speed, turb_inflow, deltav
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: turb_inflowd, deltavd
  INTEGER, INTENT(IN) :: wake_combination_method
! out
  REAL(dp), INTENT(INOUT) :: deficit_sum
  REAL(dp), DIMENSION(nbdirs), INTENT(INOUT) :: deficit_sumd
! intrinsic functions
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  INTEGER :: nd
  INTEGER :: nbdirs
! freestream linear superposition (Lissaman 1979)
  IF (wake_combination_method .EQ. 0) THEN
    DO nd=1,nbdirs
      deficit_sumd(nd) = deficit_sumd(nd) + wind_speed*deltavd(nd)
    END DO
    deficit_sum = deficit_sum + wind_speed*deltav
! local velocity linear superposition (Niayifar and Porte Agel 2015, 2016)
  ELSE IF (wake_combination_method .EQ. 1) THEN
    DO nd=1,nbdirs
      deficit_sumd(nd) = deficit_sumd(nd) + turb_inflowd(nd)*deltav + &
&       turb_inflow*deltavd(nd)
    END DO
    deficit_sum = deficit_sum + turb_inflow*deltav
!print *, "here"
! sum of squares freestream superposition (Katic et al. 1986)
  ELSE IF (wake_combination_method .EQ. 2) THEN
    arg1 = deficit_sum**2 + (wind_speed*deltav)**2
    DO nd=1,nbdirs
      arg1d(nd) = 2*deficit_sum*deficit_sumd(nd) + 2*wind_speed**2*&
&       deltav*deltavd(nd)
      IF (arg1 .EQ. 0.0) THEN
        deficit_sumd(nd) = 0.0_8
      ELSE
        deficit_sumd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    deficit_sum = SQRT(arg1)
! sum of squares local velocity superposition (Voutsinas 1990)
  ELSE IF (wake_combination_method .EQ. 3) THEN
    arg1 = deficit_sum**2 + (turb_inflow*deltav)**2
    DO nd=1,nbdirs
      arg1d(nd) = 2*deficit_sum*deficit_sumd(nd) + 2*turb_inflow*deltav*&
&       (turb_inflowd(nd)*deltav+turb_inflow*deltavd(nd))
      IF (arg1 .EQ. 0.0) THEN
        deficit_sumd(nd) = 0.0_8
      ELSE
        deficit_sumd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    deficit_sum = SQRT(arg1)
! wake combination method error
  ELSE
    PRINT*, 'Invalid wake combination method. Must be one of [0,1,2,3].'
    STOP
  END IF
END SUBROUTINE WAKE_COMBINATION_FUNC_DV

!  Differentiation of added_ti_func in forward (tangent) mode:
!   variations   of useful results: ti_dst ti_area_ratio
!   with respect to varying inputs: k_star_ust x ti_dst_in ct_ust
!                deltay ti_area_ratio_in ti_ust
! combines wakes using various methods
SUBROUTINE ADDED_TI_FUNC_DV(ti, ct_ust, ct_ustd, x, xd, k_star_ust, &
& k_star_ustd, rotor_diameter_ust, rotor_diameter_dst, deltay, deltayd, &
& wake_height, turbine_height, sm_smoothing, ti_ust, ti_ustd, &
& ti_calculation_method, ti_area_ratio_in, ti_area_ratio_ind, ti_dst_in&
& , ti_dst_ind, ti_area_ratio, ti_area_ratiod, ti_dst, ti_dstd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!print *, "ratio: ", wake_overlap/rotor_area_dst
!print *, "Dr, Dw: ", rotor_diameter_dst, wake_diameter
!print *, "Ar, Aol: ", rotor_area_dst, wake_overlap
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: ct_ust, x, k_star_ust, rotor_diameter_ust, &
& rotor_diameter_dst
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ct_ustd, xd, k_star_ustd
  REAL(dp), INTENT(IN) :: deltay, wake_height, turbine_height, &
& sm_smoothing
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: deltayd
  REAL(dp), INTENT(IN) :: ti_ust, ti, ti_area_ratio_in, ti_dst_in
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ti_ustd, &
& ti_area_ratio_ind, ti_dst_ind
  INTEGER, INTENT(IN) :: ti_calculation_method
! local
  REAL(dp) :: axial_induction_ust, beta, epsilon, sigma, wake_diameter, &
& wake_overlap
  REAL(dp), DIMENSION(nbdirs) :: axial_induction_ustd, betad, &
& epsilond, sigmad, wake_diameterd, wake_overlapd
  REAL(dp) :: ti_added, ti_tmp, rotor_area_dst, ti_area_ratio_tmp
  REAL(dp), DIMENSION(nbdirs) :: ti_addedd, ti_tmpd, &
& ti_area_ratio_tmpd
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp
! out
  REAL(dp), INTENT(OUT) :: ti_dst, ti_area_ratio
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: ti_dstd, ti_area_ratiod
! intrinsic functions
  INTRINSIC SQRT
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: result2
  REAL(dp), DIMENSION(nbdirs) :: result2d
  REAL(dp) :: pwx1
  REAL(dp), DIMENSION(nbdirs) :: pwx1d
  REAL(dp) :: pwr1
  REAL(dp), DIMENSION(nbdirs) :: pwr1d
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  INTEGER :: nd
  INTEGER :: nbdirs
  DO nd=1,nbdirs
! initialize output variables
    ti_area_ratiod(nd) = ti_area_ratio_ind(nd)
    ti_dstd(nd) = ti_dst_ind(nd)
  END DO
  ti_area_ratio = ti_area_ratio_in
  ti_dst = ti_dst_in
! initialize wake overlap to zero
  wake_overlap = 0.0_dp
!print *, "TI_dst in: ", TI_dst
! Niayifar and Porte Agel 2015, 2016 (adjusted by Annoni and Thomas for SOWFA match
! and optimization)
  IF (ti_calculation_method .EQ. 1) THEN
! calculate axial induction based on the Ct value
    CALL CT_TO_AXIAL_IND_FUNC_DV(ct_ust, ct_ustd, axial_induction_ust, &
&                          axial_induction_ustd, nbdirs)
    result1 = SQRT(1.0_dp - ct_ust)
    result2 = SQRT(1.0_dp - ct_ust)
    beta = 0.5_dp*((1.0_dp+result1)/result2)
    pwx1 = x/rotor_diameter_ust
    pwr1 = pwx1**(-0.32_dp)
    DO nd=1,nbdirs
! calculate BPA spread parameters Bastankhah and Porte Agel 2014
      result1d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      result2d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      betad(nd) = 0.5_dp*(result1d(nd)*result2-(1.0_dp+result1)*result2d&
&       (nd))/result2**2
      IF (beta .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = betad(nd)/(2.0*SQRT(beta))
      END IF
      epsilond(nd) = 0.2_dp*result1d(nd)
!print *, "epsilon = ", epsilon
! calculate wake spread for TI calcs
      sigmad(nd) = k_star_ustd(nd)*x + k_star_ust*xd(nd) + &
&       rotor_diameter_ust*epsilond(nd)
      wake_diameterd(nd) = 4.0_dp*sigmad(nd)
!print *, "wake_overlap = ", wake_overlap
! Calculate the turbulence added to the inflow of the downstream turbine by the
! wake of the upstream turbine
      pwx1d(nd) = xd(nd)/rotor_diameter_ust
      IF (pwx1 .GT. 0.0) THEN
        pwr1d(nd) = -(0.32_dp*pwx1**(-1.32)*pwx1d(nd))
      ELSE
        pwr1d(nd) = 0.0
      END IF
      ti_addedd(nd) = 0.73_dp*((0.8325_dp*axial_induction_ust**(-0.1675)&
&       *axial_induction_ustd(nd)*pwr1+axial_induction_ust**0.8325_dp*&
&       pwr1d(nd))*ti_ust**0.0325_dp+axial_induction_ust**0.8325_dp*pwr1&
&       *0.0325_dp*ti_ust**(-0.9675)*ti_ustd(nd))
    END DO
    result1 = SQRT(beta)
    epsilon = 0.2_dp*result1
    sigma = k_star_ust*x + rotor_diameter_ust*epsilon
    wake_diameter = 4.0_dp*sigma
!print *, "sigma = ", sigma
! calculate wake overlap ratio
    CALL OVERLAP_AREA_FUNC_DV(deltay, deltayd, turbine_height, &
&                       rotor_diameter_dst, 0.0_dp, wake_height, &
&                       wake_diameter, wake_diameterd, wake_overlap, &
&                       wake_overlapd, nbdirs)
    ti_added = 0.73_dp*axial_induction_ust**0.8325_dp*ti_ust**0.0325_dp*&
&     pwr1
!print *, "TI_added = ", TI_added
    rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
    arg1 = ti_dst_in**2.0_dp + (ti_added*wake_overlap/rotor_area_dst)**&
&     2.0_dp
    DO nd=1,nbdirs
! Calculate the total turbulence intensity at the downstream turbine
!sum_of_squares = TI_dst**2 + (TI_added*wake_overlap)**2
! print *, "sum of squares = ", sum_of_squares
!         TI_dst = sqrt(sum_of_squares)
!         !print *, "TI_dst = ", TI_dst
      arg1d(nd) = 2.0_dp*ti_dst_in*ti_dst_ind(nd) + 2.0_dp*ti_added*&
&       wake_overlap*(ti_addedd(nd)*wake_overlap+ti_added*wake_overlapd(&
&       nd))/rotor_area_dst**2
      IF (arg1 .EQ. 0.0) THEN
        ti_dstd(nd) = 0.0_8
      ELSE
        ti_dstd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    ti_dst = SQRT(arg1)
! Niayifar and Porte Agel 2015, 2016
  ELSE IF (ti_calculation_method .EQ. 2) THEN
! calculate axial induction based on the Ct value
    CALL CT_TO_AXIAL_IND_FUNC_DV(ct_ust, ct_ustd, axial_induction_ust, &
&                          axial_induction_ustd, nbdirs)
    result1 = SQRT(1.0_dp - ct_ust)
    result2 = SQRT(1.0_dp - ct_ust)
    beta = 0.5_dp*((1.0_dp+result1)/result2)
    pwx1 = x/rotor_diameter_ust
    pwr1 = pwx1**(-0.32_dp)
    DO nd=1,nbdirs
! calculate BPA spread parameters Bastankhah and Porte Agel 2014
      result1d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      result2d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      betad(nd) = 0.5_dp*(result1d(nd)*result2-(1.0_dp+result1)*result2d&
&       (nd))/result2**2
      IF (beta .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = betad(nd)/(2.0*SQRT(beta))
      END IF
      epsilond(nd) = 0.2_dp*result1d(nd)
! calculate wake spread for TI calcs
      sigmad(nd) = k_star_ustd(nd)*x + k_star_ust*xd(nd) + &
&       rotor_diameter_ust*epsilond(nd)
      wake_diameterd(nd) = 4.0_dp*sigmad(nd)
! Calculate the turbulence added to the inflow of the downstream turbine by the
! wake of the upstream turbine
      pwx1d(nd) = xd(nd)/rotor_diameter_ust
      IF (pwx1 .GT. 0.0) THEN
        pwr1d(nd) = -(0.32_dp*pwx1**(-1.32)*pwx1d(nd))
      ELSE
        pwr1d(nd) = 0.0
      END IF
      ti_addedd(nd) = 0.73_dp*((0.8325_dp*axial_induction_ust**(-0.1675)&
&       *axial_induction_ustd(nd)*pwr1+axial_induction_ust**0.8325_dp*&
&       pwr1d(nd))*ti_ust**0.0325_dp+axial_induction_ust**0.8325_dp*pwr1&
&       *0.0325_dp*ti_ust**(-0.9675)*ti_ustd(nd))
    END DO
    result1 = SQRT(beta)
    epsilon = 0.2_dp*result1
    sigma = k_star_ust*x + rotor_diameter_ust*epsilon
    wake_diameter = 4.0_dp*sigma
! calculate wake overlap ratio
    CALL OVERLAP_AREA_FUNC_DV(deltay, deltayd, turbine_height, &
&                       rotor_diameter_dst, 0.0_dp, wake_height, &
&                       wake_diameter, wake_diameterd, wake_overlap, &
&                       wake_overlapd, nbdirs)
    ti_added = 0.73_dp*axial_induction_ust**0.8325_dp*ti_ust**0.0325_dp*&
&     pwr1
! Calculate the total turbulence intensity at the downstream turbine based on
! current upstream turbine
    rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
    arg1 = ti**2.0_dp + (ti_added*(wake_overlap/rotor_area_dst))**2.0_dp
    DO nd=1,nbdirs
      arg1d(nd) = 2.0_dp*ti_added*wake_overlap*(ti_addedd(nd)*&
&       wake_overlap/rotor_area_dst+ti_added*wake_overlapd(nd)/&
&       rotor_area_dst)/rotor_area_dst
      IF (arg1 .EQ. 0.0) THEN
        ti_tmpd(nd) = 0.0_8
      ELSE
        ti_tmpd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    ti_tmp = SQRT(arg1)
! Check if this is the max and use it if it is
    IF (ti_tmp .GT. ti_dst_in) THEN
      DO nd=1,nbdirs
!            print *, "TI_tmp > TI_dst"
        ti_dstd(nd) = ti_tmpd(nd)
      END DO
      ti_dst = ti_tmp
    END IF
  ELSE IF (ti_calculation_method .EQ. 3) THEN
! Niayifar and Porte Agel 2015, 2016 with smooth max
! calculate axial induction based on the Ct value
    CALL CT_TO_AXIAL_IND_FUNC_DV(ct_ust, ct_ustd, axial_induction_ust, &
&                          axial_induction_ustd, nbdirs)
    result1 = SQRT(1.0_dp - ct_ust)
    result2 = SQRT(1.0_dp - ct_ust)
    beta = 0.5_dp*((1.0_dp+result1)/result2)
    pwx1 = x/rotor_diameter_ust
    pwr1 = pwx1**(-0.32_dp)
    DO nd=1,nbdirs
! calculate BPA spread parameters Bastankhah and Porte Agel 2014
      result1d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      result2d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      betad(nd) = 0.5_dp*(result1d(nd)*result2-(1.0_dp+result1)*result2d&
&       (nd))/result2**2
      IF (beta .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = betad(nd)/(2.0*SQRT(beta))
      END IF
      epsilond(nd) = 0.2_dp*result1d(nd)
! calculate wake spread for TI calcs
      sigmad(nd) = k_star_ustd(nd)*x + k_star_ust*xd(nd) + &
&       rotor_diameter_ust*epsilond(nd)
      wake_diameterd(nd) = 4.0_dp*sigmad(nd)
! Calculate the turbulence added to the inflow of the downstream turbine by the
! wake of the upstream turbine
      pwx1d(nd) = xd(nd)/rotor_diameter_ust
      IF (pwx1 .GT. 0.0) THEN
        pwr1d(nd) = -(0.32_dp*pwx1**(-1.32)*pwx1d(nd))
      ELSE
        pwr1d(nd) = 0.0
      END IF
      ti_addedd(nd) = 0.73_dp*((0.8325_dp*axial_induction_ust**(-0.1675)&
&       *axial_induction_ustd(nd)*pwr1+axial_induction_ust**0.8325_dp*&
&       pwr1d(nd))*ti_ust**0.0325_dp+axial_induction_ust**0.8325_dp*pwr1&
&       *0.0325_dp*ti_ust**(-0.9675)*ti_ustd(nd))
    END DO
    result1 = SQRT(beta)
    epsilon = 0.2_dp*result1
    sigma = k_star_ust*x + rotor_diameter_ust*epsilon
    wake_diameter = 4.0_dp*sigma
!         print *, "sigma, k_star_ust, x, rotor_diameter_ust, epsilon ", sigma, k_star_ust, x, rotor_diameter_ust, epsilon
! print *, "deltay, turbine_height, rotor_diameter_dst, wake_height, wake_diameter", &
!                 & deltay, turbine_height, rotor_diameter_dst, &
!                             wake_height, wake_diameter
! calculate wake overlap ratio
    CALL OVERLAP_AREA_FUNC_DV(deltay, deltayd, turbine_height, &
&                       rotor_diameter_dst, 0.0_dp, wake_height, &
&                       wake_diameter, wake_diameterd, wake_overlap, &
&                       wake_overlapd, nbdirs)
    ti_added = 0.73_dp*axial_induction_ust**0.8325_dp*ti_ust**0.0325_dp*&
&     pwr1
! Calculate the total turbulence intensity at the downstream turbine based on
! current upstream turbine
    rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
    arg1 = ti**2.0_dp + (ti_added*(wake_overlap/rotor_area_dst))**2.0_dp
    DO nd=1,nbdirs
      arg1d(nd) = 2.0_dp*ti_added*wake_overlap*(ti_addedd(nd)*&
&       wake_overlap/rotor_area_dst+ti_added*wake_overlapd(nd)/&
&       rotor_area_dst)/rotor_area_dst
      IF (arg1 .EQ. 0.0) THEN
        ti_tmpd(nd) = 0.0_8
      ELSE
        ti_tmpd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    ti_tmp = SQRT(arg1)
!print *, "TI, TI_added, wake_overlap, rotor_area_dst: ", TI, TI_added, wake_overlap, rotor_area_dst
! Check if this is the max and use it if it is
!if (TI_tmp > TI_dst) then
!    TI_dst = TI_tmp
!end if
!         print *, "before: ", TI_dst, TI_tmp
!         TI_dst_in = TI_dst
    CALL SMOOTH_MAX_DV(sm_smoothing, ti_dst_in, ti_dst_ind, ti_tmp, &
&                ti_tmpd, ti_dst, ti_dstd, nbdirs)
!         print *, "after:: ", TI_dst, TI_tmp
! Niayifar and Porte Agel 2015, 2016 using max on area TI ratio
  ELSE IF (ti_calculation_method .EQ. 4) THEN
! calculate axial induction based on the Ct value
    CALL CT_TO_AXIAL_IND_FUNC_DV(ct_ust, ct_ustd, axial_induction_ust, &
&                          axial_induction_ustd, nbdirs)
    result1 = SQRT(1.0_dp - ct_ust)
    result2 = SQRT(1.0_dp - ct_ust)
    beta = 0.5_dp*((1.0_dp+result1)/result2)
    pwx1 = x/rotor_diameter_ust
    pwr1 = pwx1**(-0.32_dp)
    DO nd=1,nbdirs
! calculate BPA spread parameters Bastankhah and Porte Agel 2014
      result1d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      result2d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      betad(nd) = 0.5_dp*(result1d(nd)*result2-(1.0_dp+result1)*result2d&
&       (nd))/result2**2
      IF (beta .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = betad(nd)/(2.0*SQRT(beta))
      END IF
      epsilond(nd) = 0.2_dp*result1d(nd)
! calculate wake spread for TI calcs
      sigmad(nd) = k_star_ustd(nd)*x + k_star_ust*xd(nd) + &
&       rotor_diameter_ust*epsilond(nd)
      wake_diameterd(nd) = 4.0_dp*sigmad(nd)
! Calculate the turbulence added to the inflow of the downstream turbine by the
! wake of the upstream turbine
      pwx1d(nd) = xd(nd)/rotor_diameter_ust
      IF (pwx1 .GT. 0.0) THEN
        pwr1d(nd) = -(0.32_dp*pwx1**(-1.32)*pwx1d(nd))
      ELSE
        pwr1d(nd) = 0.0
      END IF
      ti_addedd(nd) = 0.73_dp*((0.8325_dp*axial_induction_ust**(-0.1675)&
&       *axial_induction_ustd(nd)*pwr1+axial_induction_ust**0.8325_dp*&
&       pwr1d(nd))*ti_ust**0.0325_dp+axial_induction_ust**0.8325_dp*pwr1&
&       *0.0325_dp*ti_ust**(-0.9675)*ti_ustd(nd))
    END DO
    result1 = SQRT(beta)
    epsilon = 0.2_dp*result1
    sigma = k_star_ust*x + rotor_diameter_ust*epsilon
    wake_diameter = 4.0_dp*sigma
! calculate wake overlap ratio
    CALL OVERLAP_AREA_FUNC_DV(deltay, deltayd, turbine_height, &
&                       rotor_diameter_dst, 0.0_dp, wake_height, &
&                       wake_diameter, wake_diameterd, wake_overlap, &
&                       wake_overlapd, nbdirs)
    ti_added = 0.73_dp*axial_induction_ust**0.8325_dp*ti_ust**0.0325_dp*&
&     pwr1
! Calculate the total turbulence intensity at the downstream turbine based on
! current upstream turbine
    rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
    DO nd=1,nbdirs
      ti_area_ratio_tmpd(nd) = ti_addedd(nd)*wake_overlap/rotor_area_dst&
&       + ti_added*wake_overlapd(nd)/rotor_area_dst
    END DO
    ti_area_ratio_tmp = ti_added*(wake_overlap/rotor_area_dst)
! Check if this is the max and use it if it is
    IF (ti_area_ratio_tmp .GT. ti_area_ratio_in) THEN
      ti_area_ratio = ti_area_ratio_tmp
      arg1 = ti**2.0_dp + ti_area_ratio**2.0_dp
      DO nd=1,nbdirs
!            print *, "ti_area_ratio_tmp > ti_area_ratio"
!TI_dst = TI_tmp
        ti_area_ratiod(nd) = ti_area_ratio_tmpd(nd)
        arg1d(nd) = 2.0_dp*ti_area_ratio*ti_area_ratiod(nd)
        IF (arg1 .EQ. 0.0) THEN
          ti_dstd(nd) = 0.0_8
        ELSE
          ti_dstd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
        END IF
      END DO
      ti_dst = SQRT(arg1)
    END IF
  ELSE IF (ti_calculation_method .EQ. 5) THEN
! Niayifar and Porte Agel 2015, 2016 using smooth max on area TI ratio
! calculate axial induction based on the Ct value
    CALL CT_TO_AXIAL_IND_FUNC_DV(ct_ust, ct_ustd, axial_induction_ust, &
&                          axial_induction_ustd, nbdirs)
    result1 = SQRT(1.0_dp - ct_ust)
    result2 = SQRT(1.0_dp - ct_ust)
    beta = 0.5_dp*((1.0_dp+result1)/result2)
    DO nd=1,nbdirs
! calculate BPA spread parameters Bastankhah and Porte Agel 2014
      result1d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      result2d(nd) = -(ct_ustd(nd)/(2.0*SQRT(1.0_dp-ct_ust)))
      betad(nd) = 0.5_dp*(result1d(nd)*result2-(1.0_dp+result1)*result2d&
&       (nd))/result2**2
      IF (beta .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = betad(nd)/(2.0*SQRT(beta))
      END IF
      epsilond(nd) = 0.2_dp*result1d(nd)
! calculate wake spread for TI calcs
      sigmad(nd) = k_star_ustd(nd)*x + k_star_ust*xd(nd) + &
&       rotor_diameter_ust*epsilond(nd)
      wake_diameterd(nd) = 4.0_dp*sigmad(nd)
    END DO
    result1 = SQRT(beta)
    epsilon = 0.2_dp*result1
    sigma = k_star_ust*x + rotor_diameter_ust*epsilon
    wake_diameter = 4.0_dp*sigma
! calculate wake overlap ratio
    CALL OVERLAP_AREA_FUNC_DV(deltay, deltayd, turbine_height, &
&                       rotor_diameter_dst, 0.0_dp, wake_height, &
&                       wake_diameter, wake_diameterd, wake_overlap, &
&                       wake_overlapd, nbdirs)
! only include turbines with area overlap in the softmax
    IF (wake_overlap .GT. 0.0_dp) THEN
      pwx1 = x/rotor_diameter_ust
      pwr1 = pwx1**(-0.32_dp)
      ti_added = 0.73_dp*axial_induction_ust**0.8325_dp*ti_ust**&
&       0.0325_dp*pwr1
      rotor_area_dst = 0.25_dp*pi*rotor_diameter_dst**2_dp
      DO nd=1,nbdirs
! Calculate the turbulence added to the inflow of the downstream turbine by the
! wake of the upstream turbine
        pwx1d(nd) = xd(nd)/rotor_diameter_ust
        IF (pwx1 .GT. 0.0) THEN
          pwr1d(nd) = -(0.32_dp*pwx1**(-1.32)*pwx1d(nd))
        ELSE
          pwr1d(nd) = 0.0
        END IF
        ti_addedd(nd) = 0.73_dp*((0.8325_dp*axial_induction_ust**(&
&         -0.1675)*axial_induction_ustd(nd)*pwr1+axial_induction_ust**&
&         0.8325_dp*pwr1d(nd))*ti_ust**0.0325_dp+axial_induction_ust**&
&         0.8325_dp*pwr1*0.0325_dp*ti_ust**(-0.9675)*ti_ustd(nd))
        ti_area_ratio_tmpd(nd) = ti_addedd(nd)*wake_overlap/&
&         rotor_area_dst + ti_added*wake_overlapd(nd)/rotor_area_dst
      END DO
      ti_area_ratio_tmp = ti_added*(wake_overlap/rotor_area_dst)
!TI_tmp = sqrt(TI**2.0_dp + (TI_added*(wake_overlap/rotor_area_dst))**2.0_dp)
! Run through the smooth max to get an approximation of the true max TI area ratio
      CALL SMOOTH_MAX_DV(sm_smoothing, ti_area_ratio_in, &
&                  ti_area_ratio_ind, ti_area_ratio_tmp, &
&                  ti_area_ratio_tmpd, ti_area_ratio, ti_area_ratiod, &
&                  nbdirs)
      arg1 = ti**2.0_dp + ti_area_ratio**2.0_dp
      DO nd=1,nbdirs
! Calculate the total turbulence intensity at the downstream turbine based on
! the result of the smooth max function
        arg1d(nd) = 2.0_dp*ti_area_ratio*ti_area_ratiod(nd)
        IF (arg1 .EQ. 0.0) THEN
          ti_dstd(nd) = 0.0_8
        ELSE
          ti_dstd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
        END IF
      END DO
      ti_dst = SQRT(arg1)
    END IF
  ELSE
! wake combination method error
    PRINT*, &
&   'Invalid added TI calculation method. Must be one of [0,1,2,3,4,5].'
    STOP
  END IF
END SUBROUTINE ADDED_TI_FUNC_DV

!  Differentiation of overlap_area_func in forward (tangent) mode:
!   variations   of useful results: wake_overlap
!   with respect to varying inputs: wake_diameter turbine_y
! calculates the overlap area between a given wake and a rotor area
SUBROUTINE OVERLAP_AREA_FUNC_DV(turbine_y, turbine_yd, turbine_z, &
& rotor_diameter, wake_center_y, wake_center_z, wake_diameter, &
& wake_diameterd, wake_overlap, wake_overlapd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: turbine_y, turbine_z, rotor_diameter
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: turbine_yd
  REAL(dp), INTENT(IN) :: wake_center_y, wake_center_z, wake_diameter
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: wake_diameterd
! out
  REAL(dp), INTENT(OUT) :: wake_overlap
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: wake_overlapd
! local
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, tol=0.000001_dp
  REAL(dp) :: ovdyd, ovr, ovrr, ovl, ovz, ovz2
  REAL(dp), DIMENSION(nbdirs) :: ovdydd, ovrrd, ovld, ovzd, ovz2d
! load intrinsic functions
  INTRINSIC ACOS, SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  REAL(dp) :: result2
  REAL(dp), DIMENSION(nbdirs) :: result2d
  INTEGER :: nd
  INTEGER :: nbdirs
!     print *, turbine_y, turbine_z, rotor_diameter, &
!                             wake_center_y, wake_center_z, wake_diameter, &
!                             wake_overlap
! distance between wake center and rotor center
  IF (wake_center_z .GT. turbine_z + tol .OR. wake_center_z .LT. &
&     turbine_z - tol) THEN
    arg1 = (wake_center_y-turbine_y)**2_dp + (wake_center_z-turbine_z)**&
&     2_dp
    DO nd=1,nbdirs
      arg1d(nd) = -(2_dp*(wake_center_y-turbine_y)*turbine_yd(nd))
      IF (arg1 .EQ. 0.0) THEN
        ovdydd(nd) = 0.0_8
      ELSE
        ovdydd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
    END DO
    ovdyd = SQRT(arg1)
  ELSE IF (wake_center_y .GT. turbine_y + tol) THEN
    DO nd=1,nbdirs
! potential source of gradient issues, abs() did not cause a problem in FLORIS
      ovdydd(nd) = -turbine_yd(nd)
    END DO
    ovdyd = wake_center_y - turbine_y
  ELSE IF (turbine_y .GT. wake_center_y + tol) THEN
    DO nd=1,nbdirs
      ovdydd(nd) = turbine_yd(nd)
    END DO
    ovdyd = turbine_y - wake_center_y
  ELSE
    ovdyd = 0.0_dp
    ovdydd(:) = 0.0_8
  END IF
!print *, "OVdYd: ", OVdYd
! find rotor radius
  ovr = rotor_diameter/2.0_dp
  DO nd=1,nbdirs
!print *, "OVr: ", OVr
! find wake radius
    ovrrd(nd) = wake_diameterd(nd)/2.0_dp
  END DO
  ovrr = wake_diameter/2.0_dp
!print *, "OVRR: ", OVRR
! determine if there is overlap
  IF (ovdyd .LT. ovr + ovrr) THEN
! if the rotor overlaps the wake zone
! check that turbine and wake centers are not perfectly aligned
    IF (ovdyd .GT. 0.0_dp + tol) THEN
! check if the rotor is wholly contained in the wake
      IF (ovdyd + ovr .LT. ovrr + tol) THEN
        wake_overlap = pi*ovr*ovr
!                 print *, "1"
! check if the wake is wholly contained in the rotor swept area
        wake_overlapd(:) = 0.0_8
      ELSE IF (ovdyd + ovrr .LT. ovr + tol) THEN
        DO nd=1,nbdirs
          wake_overlapd(nd) = pi*(ovrrd(nd)*ovrr+ovrr*ovrrd(nd))
        END DO
        wake_overlap = pi*ovrr*ovrr
!                 print *, "2"
      ELSE
        ovl = (-(ovr*ovr)+ovrr*ovrr+ovdyd*ovdyd)/(2.0_dp*ovdyd)
        arg1 = ovrr*ovrr - ovl*ovl
        arg2 = (ovdyd-ovl)/ovr
        DO nd=1,nbdirs
! calculate the distance from the wake center to the chord connecting the lens
! cusps
          ovld(nd) = ((ovrrd(nd)*ovrr+ovrr*ovrrd(nd)+ovdydd(nd)*ovdyd+&
&           ovdyd*ovdydd(nd))*2.0_dp*ovdyd-(-(ovr*ovr)+ovrr*ovrr+ovdyd*&
&           ovdyd)*2.0_dp*ovdydd(nd))/(2.0_dp*ovdyd)**2
          arg1d(nd) = ovrrd(nd)*ovrr + ovrr*ovrrd(nd) - ovld(nd)*ovl - &
&           ovl*ovld(nd)
          IF (arg1 .EQ. 0.0) THEN
            ovzd(nd) = 0.0_8
          ELSE
            ovzd(nd) = arg1d(nd)/(2.0*SQRT(arg1))
          END IF
          arg1d(nd) = -((ovdydd(nd)-ovld(nd))*(ovdyd-ovl)) - (ovdyd-ovl)&
&           *(ovdydd(nd)-ovld(nd))
          arg2d(nd) = (ovdydd(nd)-ovld(nd))/ovr
          IF (arg2 .EQ. 1.0 .OR. arg2 .EQ. (-1.0)) THEN
            result2d(nd) = 0.0_8
          ELSE
            result2d(nd) = -(arg2d(nd)/SQRT(1.0-arg2**2))
          END IF
        END DO
        ovz = SQRT(arg1)
        arg1 = ovr*ovr - (ovdyd-ovl)*(ovdyd-ovl)
        DO nd=1,nbdirs
          IF (arg1 .EQ. 0.0) THEN
            ovz2d(nd) = 0.0_8
          ELSE
            ovz2d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
          END IF
          arg1d(nd) = (ovld(nd)*ovrr-ovl*ovrrd(nd))/ovrr**2
        END DO
        ovz2 = SQRT(arg1)
        arg1 = ovl/ovrr
        result1 = ACOS(arg1)
        DO nd=1,nbdirs
          IF (arg1 .EQ. 1.0 .OR. arg1 .EQ. (-1.0)) THEN
            result1d(nd) = 0.0_8
          ELSE
            result1d(nd) = -(arg1d(nd)/SQRT(1.0-arg1**2))
          END IF
          wake_overlapd(nd) = (ovrrd(nd)*ovrr+ovrr*ovrrd(nd))*result1 + &
&           ovrr**2*result1d(nd) + ovr**2*result2d(nd) - ovld(nd)*ovz - &
&           ovl*ovzd(nd) - (ovdydd(nd)-ovld(nd))*ovz2 - (ovdyd-ovl)*&
&           ovz2d(nd)
        END DO
        result2 = ACOS(arg2)
        wake_overlap = ovrr*ovrr*result1 + ovr*ovr*result2 - ovl*ovz - (&
&         ovdyd-ovl)*ovz2
!                 print *, OVRR, OVr, OVdYd, OVL, OVz, OVz2
!                 print *, "3"
      END IF
    ELSE IF (ovrr .GT. ovr) THEN
! perfect overlap case where the wake is larger than the rotor
      wake_overlap = pi*ovr*ovr
!             print *, "4"
! perfect overlap case where the rotor is larger than the wake
      wake_overlapd(:) = 0.0_8
    ELSE
      DO nd=1,nbdirs
        wake_overlapd(nd) = pi*(ovrrd(nd)*ovrr+ovrr*ovrrd(nd))
      END DO
      wake_overlap = pi*ovrr*ovrr
!             print *, "5"
    END IF
  ELSE
! case with no overlap
    wake_overlap = 0.0_dp
    wake_overlapd(:) = 0.0_8
  END IF
!     print *, "wake overlap in func: ", wake_overlap/(pi*OVr**2)
!     print *, "wake overlap in func: ", wake_overlap/(pi*OVRR**2)
  IF (wake_overlap/(pi*ovr*ovr) .GT. 1.0_dp + tol .OR. wake_overlap/(pi*&
&     ovrr*ovrr) .GT. 1.0_dp + tol) THEN
    PRINT*, 'wake overlap in func: ', wake_overlap/(pi*ovr*ovr)
    PRINT*, 'wake overlap in func: ', wake_overlap/(pi*ovrr*ovrr)
    STOP
  END IF
END SUBROUTINE OVERLAP_AREA_FUNC_DV

!  Differentiation of k_star_func in forward (tangent) mode:
!   variations   of useful results: k_star_ust
!   with respect to varying inputs: ti_ust
! compute wake spread parameter based on local turbulence intensity
SUBROUTINE K_STAR_FUNC_DV(ti_ust, ti_ustd, k_star_ust, k_star_ustd, &
& nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: ti_ust
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ti_ustd
! out
  REAL(dp), INTENT(OUT) :: k_star_ust
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: k_star_ustd
  INTEGER :: nd
  INTEGER :: nbdirs
  DO nd=1,nbdirs
! calculate wake spread parameter from Niayifar and Porte Agel (2015, 2016)
    k_star_ustd(nd) = 0.3837*ti_ustd(nd)
  END DO
  k_star_ust = 0.3837*ti_ust + 0.003678
END SUBROUTINE K_STAR_FUNC_DV

!  Differentiation of ct_to_axial_ind_func in forward (tangent) mode:
!   variations   of useful results: axial_induction
!   with respect to varying inputs: ct
! calculate axial induction from Ct
SUBROUTINE CT_TO_AXIAL_IND_FUNC_DV(ct, ctd, axial_induction, &
& axial_inductiond, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: ct
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: ctd
! out
  REAL(dp), INTENT(OUT) :: axial_induction
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: axial_inductiond
  INTRINSIC SQRT
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  INTEGER :: nd
  INTEGER :: nbdirs
! initialize axial induction to zero
  axial_induction = 0.0_dp
! calculate axial induction
  IF (ct .GT. 0.96) THEN
    arg1 = 0.0203_dp - 0.6427_dp*(0.889_dp-ct)
    DO nd=1,nbdirs
! Glauert condition
      arg1d(nd) = 0.6427_dp*ctd(nd)
      IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.0_8
      ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
      axial_inductiond(nd) = result1d(nd)
    END DO
    result1 = SQRT(arg1)
    axial_induction = 0.143_dp + result1
  ELSE
    DO nd=1,nbdirs
      result1d(nd) = -(ctd(nd)/(2.0*SQRT(1.0_dp-ct)))
      axial_inductiond(nd) = -(0.5_dp*result1d(nd))
    END DO
    result1 = SQRT(1.0_dp - ct)
    axial_induction = 0.5_dp*(1.0_dp-result1)
  END IF
END SUBROUTINE CT_TO_AXIAL_IND_FUNC_DV

!  Differentiation of wind_shear_func in forward (tangent) mode:
!   variations   of useful results: adjusted_wind_speed
!   with respect to varying inputs: u_ref
! adjust wind speed for wind shear
SUBROUTINE WIND_SHEAR_FUNC_DV(point_z, u_ref, u_refd, z_ref, z_0, &
& shear_exp, adjusted_wind_speed, adjusted_wind_speedd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: point_z, u_ref, z_ref, z_0, shear_exp
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: u_refd
! out
  REAL(dp), INTENT(OUT) :: adjusted_wind_speed
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: adjusted_wind_speedd
  REAL(dp) :: pwx1
  REAL(dp) :: pwr1
  INTEGER :: nd
  INTEGER :: nbdirs
! initialize adjusted wind speed to zero
  adjusted_wind_speed = 0.0_dp
! check that the point of interest is above ground level
  IF (point_z .GE. z_0) THEN
! adjusted wind speed for wind shear if point is above ground
    pwx1 = (point_z-z_0)/(z_ref-z_0)
    pwr1 = pwx1**shear_exp
    DO nd=1,nbdirs
      adjusted_wind_speedd(nd) = pwr1*u_refd(nd)
    END DO
    adjusted_wind_speed = u_ref*pwr1
  ELSE
! if the point of interest is below ground, set the wind speed to 0.0
    adjusted_wind_speed = 0.0_dp
    adjusted_wind_speedd(:) = 0.0_8
  END IF
END SUBROUTINE WIND_SHEAR_FUNC_DV

!  Differentiation of discontinuity_point_func in forward (tangent) mode:
!   variations   of useful results: discontinuity_point
!   with respect to varying inputs: ky kz x0 ct
! calculate the point where the Bastankhah and Porte Agel wake model becomes undefined
SUBROUTINE DISCONTINUITY_POINT_FUNC_DV(x0, x0d, rotor_diameter, ky, kyd&
& , kz, kzd, yaw, ct, ctd, discontinuity_point, discontinuity_pointd, &
& nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x0, rotor_diameter, ky, kz, yaw, ct
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: x0d, kyd, kzd, ctd
! local
  REAL(dp) :: a, b, c
  REAL(dp), DIMENSION(nbdirs) :: ad, bd, cd
! out
  REAL(dp), INTENT(OUT) :: discontinuity_point
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: discontinuity_pointd
  INTRINSIC COS, SQRT
  REAL(dp) :: result1
  REAL(dp), DIMENSION(nbdirs) :: result1d
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  INTEGER :: nd
  INTEGER :: nbdirs
  a = ky + kz*COS(yaw)
  b = 4.0_dp*ky*kz*COS(yaw)*(ct-1.0_dp)
  result1 = SQRT(8.0_dp)
  arg1 = a**2 - b
  DO nd=1,nbdirs
! for clarity, break out the terms in the equation
    ad(nd) = kyd(nd) + COS(yaw)*kzd(nd)
    bd(nd) = 4.0_dp*COS(yaw)*((kyd(nd)*kz+ky*kzd(nd))*(ct-1.0_dp)+ky*kz*&
&     ctd(nd))
    cd(nd) = 2.0_dp*result1*(kyd(nd)*kz+ky*kzd(nd))
! distance from rotor to the last point where the wake model is undefined
    arg1d(nd) = 2*a*ad(nd) - bd(nd)
    IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.0_8
    ELSE
      result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
  END DO
  c = 2.0_dp*result1*ky*kz
  result1 = SQRT(arg1)
  DO nd=1,nbdirs
    discontinuity_pointd(nd) = x0d(nd) + (rotor_diameter*(ad(nd)-&
&     result1d(nd))*c-rotor_diameter*(a-result1)*cd(nd))/c**2
  END DO
  discontinuity_point = x0 + rotor_diameter*(a-result1)/c
END SUBROUTINE DISCONTINUITY_POINT_FUNC_DV

!  Differentiation of smooth_max in forward (tangent) mode:
!   variations   of useful results: g
!   with respect to varying inputs: x y
SUBROUTINE SMOOTH_MAX_DV(s, x, xd, y, yd, g, gd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: s, x, y
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: xd, yd
! local
  REAL(dp) :: max_val, min_val
  REAL(dp), DIMENSION(nbdirs) :: max_vald, min_vald
! out
  REAL(dp), INTENT(OUT) :: g
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: gd
  INTRINSIC LOG, EXP, MAX, MIN
  REAL(dp) :: arg1
  REAL(dp), DIMENSION(nbdirs) :: arg1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  INTEGER :: nd
  INTEGER :: nbdirs
  IF (x .LT. y) THEN
    DO nd=1,nbdirs
      max_vald(nd) = yd(nd)
    END DO
    max_val = y
  ELSE
    DO nd=1,nbdirs
      max_vald(nd) = xd(nd)
    END DO
    max_val = x
  END IF
  IF (x .GT. y) THEN
    DO nd=1,nbdirs
      min_vald(nd) = yd(nd)
    END DO
    min_val = y
  ELSE
    DO nd=1,nbdirs
      min_vald(nd) = xd(nd)
    END DO
    min_val = x
  END IF
  arg1 = s*(min_val-max_val)
  arg2 = 1.0_dp + EXP(arg1)
  DO nd=1,nbdirs
    arg1d(nd) = s*(min_vald(nd)-max_vald(nd))
    arg2d(nd) = arg1d(nd)*EXP(arg1)
    gd(nd) = (arg2d(nd)/arg2+s*max_vald(nd))/s
  END DO
  g = (LOG(arg2)+s*max_val)/s
END SUBROUTINE SMOOTH_MAX_DV

!  Differentiation of interpolation in forward (tangent) mode:
!   variations   of useful results: yval
!   with respect to varying inputs: yval xval
SUBROUTINE INTERPOLATION_DV(npoints, interp_type, x, y, xval, xvald, &
& yval, yvald, dy0in, dy1in, usedyin, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!     print *, "yval = ", yval
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: npoints, interp_type
  REAL(dp), DIMENSION(npoints), INTENT(IN) :: x, y
  REAL(dp), INTENT(IN) :: xval
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: xvald
  REAL(dp), INTENT(IN) :: dy0in, dy1in
  LOGICAL :: usedyin
! local
  INTEGER :: idx
  REAL(dp) :: x0, x1, y0, dy0, y1, dy1
! out
  REAL(dp), INTENT(OUT) :: yval
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: yvald
  INTEGER :: nd
  INTEGER :: nbdirs
!     print *, "in interpolation"
! if ((xval < x(1)) .or. (xval > x(nPoints))) then
!         print *, "interpolation point is out of bounds"
! !         STOP 1
!     end if
  IF (usedyin .AND. interp_type .EQ. 1) THEN
    PRINT*, &
&  'end point derivatives may not be specified for linear interpolation'
    STOP
  ELSE IF (xval .LT. x(1)) THEN
    yval = y(1)
    yvald(:) = 0.0_8
  ELSE IF (xval .GT. x(npoints)) THEN
    yval = y(npoints)
    yvald(:) = 0.0_8
  ELSE
    idx = 1
    DO WHILE (xval .GT. x(idx) .AND. idx .LE. npoints)
      idx = idx + 1
    END DO
    idx = idx - 1
    x0 = x(idx)
    x1 = x(idx+1)
    y0 = y(idx)
    y1 = y(idx+1)
! Hermite cubic piecewise spline interpolation
    IF (interp_type .EQ. 0) THEN
! approximate derivative at left end of interval
      IF (idx .EQ. 1) THEN
        IF (usedyin) THEN
          dy0 = dy0in
        ELSE
          dy0 = 0.0_dp
        END IF
      ELSE
        dy0 = (y(idx+1)-y(idx-1))/(x(idx+1)-x(idx-1))
      END IF
! approximate derivative at the right end of interval
      IF (idx .GE. npoints - 1) THEN
        IF (usedyin) THEN
          dy1 = dy1in
        ELSE
          dy1 = 0.0_dp
        END IF
      ELSE
        dy1 = (y(idx+2)-y(idx))/(x(idx+2)-x(idx))
      END IF
! call Hermite spline routine
      CALL HERMITE_SPLINE_DV(xval, xvald, x0, x1, y0, dy0, y1, dy1, yval&
&                      , yvald, nbdirs)
! linear interpolation
    ELSE IF (interp_type .EQ. 1) THEN
      DO nd=1,nbdirs
        yvald(nd) = (y1-y0)*xvald(nd)/(x1-x0)
      END DO
      yval = (xval-x0)*(y1-y0)/(x1-x0) + y0
    END IF
  END IF
END SUBROUTINE INTERPOLATION_DV

!  Differentiation of hermite_spline in forward (tangent) mode:
!   variations   of useful results: y
!   with respect to varying inputs: x
SUBROUTINE HERMITE_SPLINE_DV(x, xd, x0, x1, y0, dy0, y1, dy1, y, yd, &
& nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
!dy_dx = c3*3*x**2 + c2*2*x + c1
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  REAL(dp), INTENT(IN) :: x, x0, x1, y0, dy0, y1, dy1
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: xd
! out
!, dy_dx
  REAL(dp), INTENT(OUT) :: y
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: yd
! local
  REAL(dp) :: c3, c2, c1, c0
  INTEGER :: nd
  INTEGER :: nbdirs
! initialize coefficients for parametric cubic spline
  c3 = 2.0_dp*y1/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) - 2.0_dp*&
&   y0/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + dy0/(x0**2-2.0_dp&
&   *x0*x1+x1**2) + dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  c2 = 3.0_dp*y0*(x0+x1)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) -&
&   dy1*(2.0_dp*x0+x1)/(x0**2-2.0_dp*x0*x1+x1**2) - dy0*(x0+2.0_dp*x1)/(&
&   x0**2-2.0_dp*x0*x1+x1**2) - 3.0_dp*y1*(x0+x1)/(x0**3-3.0_dp*x0**2*x1&
&   +3.0_dp*x0*x1**2-x1**3)
  c1 = dy0*(x1**2+2.0_dp*x0*x1)/(x0**2-2.0_dp*x0*x1+x1**2) + dy1*(x0**2+&
&   2.0_dp*x1*x0)/(x0**2-2.0_dp*x0*x1+x1**2) - 6.0_dp*x0*x1*y0/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3) + 6.0_dp*x0*x1*y1/(x0**3-&
&   3.0_dp*x0**2*x1+3.0_dp*x0*x1**2-x1**3)
  c0 = y0*(-(x1**3)+3.0_dp*x0*x1**2)/(x0**3-3.0_dp*x0**2*x1+3.0_dp*x0*x1&
&   **2-x1**3) - y1*(-(x0**3)+3.0_dp*x1*x0**2)/(x0**3-3.0_dp*x0**2*x1+&
&   3.0_dp*x0*x1**2-x1**3) - x0*x1**2*dy0/(x0**2-2.0_dp*x0*x1+x1**2) - &
&   x0**2*x1*dy1/(x0**2-2.0_dp*x0*x1+x1**2)
  DO nd=1,nbdirs
! Solve for y and dy values at the given point
    yd(nd) = c3*3*x**2*xd(nd) + c2*2*x*xd(nd) + c1*xd(nd)
  END DO
  y = c3*x**3 + c2*x**2 + c1*x + c0
END SUBROUTINE HERMITE_SPLINE_DV
