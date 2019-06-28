
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


subroutine makeGrid_fortran(nRows, nTurbines, dx, dy, shear, rotate, turbs_per_row, x_start, &
                    & y0, turbineX, turbineY)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nRows, nTurbines
    real(dp), intent(in) :: dx, dy, shear, rotate, y0
    integer, dimension(nRows), intent(in) :: turbs_per_row
    real(dp), dimension(nRows), intent(in) ::  x_start

    ! out
    real(dp), dimension(nTurbines), intent(out) :: turbineX, turbineY

    ! local
    integer :: i, j, index
    real(dp) :: rotate_rad
    real(dp), dimension(nTurbines) :: x, y

    index = 1

    do i = 1, nRows
        do j = 1, turbs_per_row(i)
            x(index) = x_start(i) + dx*j + i*shear
            y(index) = y0 + dy*i
            index = index + 1
        end do
    end do

    rotate_rad = (rotate*3.1415926535)/180.0
    turbineX = cos(rotate_rad)*x - sin(rotate_rad)*y
    turbineY = sin(rotate_rad)*x + cos(rotate_rad)*y


end subroutine makeGrid_fortran


subroutine calcAEP_grid(nTurbines, nDirections, nRows, turbineZ, rotorDiameter, windDirections,&
            &dx, dy, shear, rotate, turbs_per_row, x_start, y0,&
            &windSpeeds, windFrequencies, shearExp, relaxationFactor, rated_ws, rated_power,&
            &cut_in_speed, cut_out_speed, zref, z0, AEP)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines, nDirections, nRows
    real(dp), intent(in) :: shearExp, relaxationFactor, rated_ws, rated_power, cut_in_speed,&
              &cut_out_speed, zref, z0, dx, dy, shear, rotate, y0
    real(dp), dimension(nTurbines), intent(in) :: turbineZ, rotorDiameter
    real(dp), dimension(nDirections), intent(in) :: windDirections, windSpeeds, windFrequencies
    integer, dimension(nRows), intent(in) :: turbs_per_row
    real(dp), dimension(nRows), intent(in) ::  x_start

    ! out
    real(dp), intent(out) :: AEP

    ! local
    real(dp), dimension(nTurbines) :: turbineX, turbineY
    real(dp), dimension(nDirections) :: dir_powers
    real(dp), dimension(nTurbines) :: turbineXw, turbineYw, Vinf_floris, wtVelocity, loss
    real(dp) :: hrs_per_year, pwrDir, Vinf
    integer :: n, i

    call makeGrid_fortran(nRows, nTurbines, dx, dy, shear, rotate, turbs_per_row, x_start, &
                        & y0, turbineX, turbineY)

    do n = 1, nDirections
        call WindFrame(nTurbines, windDirections(n), turbineX, turbineY, turbineXw, turbineYw)
        call PowWind(nTurbines, windSpeeds(n), turbineZ, shearExp, zref, z0, Vinf_floris)
        Vinf = Vinf_floris(1)
        call GaussianWake(nTurbines, turbineXw, turbineYw, rotorDiameter(1), relaxationFactor, loss)
        wtVelocity = Vinf*(1.0_dp-loss)
        call DirPower(nTurbines, wtVelocity, rated_ws, rated_power, cut_in_speed, cut_out_speed, pwrDir)
        dir_powers(n) = pwrDir
    end do

    hrs_per_year = 365.*24.
    AEP = hrs_per_year * (sum(windFrequencies * dir_powers))

end subroutine calcAEP_grid


subroutine GaussianWake(nTurbines, turbineXw, turbineYw, turb_diam, relaxationFactor, loss)

    implicit none

    ! define precision to be the standard for a double precision ! on local system
    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: nTurbines
    real(dp), intent(in) :: turb_diam, relaxationFactor
    real(dp), dimension(nTurbines), intent(in) :: turbineXw, turbineYw

    ! out
    real(dp), dimension(nTurbines), intent(out) :: loss

    ! local
    real(dp) :: CT, k, x, y, sigma, exponent, radical
    real(dp), dimension(nTurbines) :: loss_array
    real(dp), parameter :: pi = 3.141592653589793_dp, e = 2.718281828459045_dp, tol = 0.000001_dp
    integer :: i, j

    CT = 4.0*1./3.*(1.0-1./3.)
    k = 0.0324555

    do i = 1, nTurbines
        do j = 1, nTurbines
            x = turbineXw(i) - turbineXw(j)
            y = turbineYw(i) - turbineYw(j)
            if (x > 0.0) then
                sigma = k*x + turb_diam/sqrt(8.0)
                exponent = -0.5 * (y/(relaxationFactor*sigma))**2
                radical = 1. - CT/(8.*sigma**2 / turb_diam**2)
                loss_array(j) = (1.-sqrt(radical)) * e**exponent
            else
                loss_array(j) = 0.0
            end if
        end do
        loss(i) = sqrt(sum(loss_array**2))
    end do

end subroutine GaussianWake





!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7259) - 18 Jan 2019 09:31
!
!  Differentiation of windframe in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
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

!  Differentiation of hermite_spline in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
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

!  Differentiation of dirpower in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
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

!  Differentiation of makegrid_fortran in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: turbinex turbiney
!   with respect to varying inputs: rotate dx dy shear
SUBROUTINE MAKEGRID_FORTRAN_DV(nrows, nturbines, dx, dxd, dy, dyd, shear&
& , sheard, rotate, rotated, turbs_per_row, x_start, y0, turbinex, &
& turbinexd, turbiney, turbineyd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nrows, nturbines
  REAL(dp), INTENT(IN) :: dx, dy, shear, rotate, y0
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: dxd, dyd, sheard, &
& rotated
  INTEGER, DIMENSION(nrows), INTENT(IN) :: turbs_per_row
  REAL(dp), DIMENSION(nrows), INTENT(IN) :: x_start
! out
  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: turbinex, turbiney
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: turbinexd, &
& turbineyd
! local
  INTEGER :: i, j, index
  REAL(dp) :: rotate_rad
  REAL(dp), DIMENSION(nbdirs) :: rotate_radd
  REAL(dp), DIMENSION(nturbines) :: x, y
  REAL(dp), DIMENSION(nbdirs, nturbines) :: xd, yd
  INTRINSIC COS
  INTRINSIC SIN
  INTEGER :: nd
  INTEGER :: nbdirs
  index = 1
  xd(:, :) = 0.0_8
  yd(:, :) = 0.0_8
  DO i=1,nrows
    DO j=1,turbs_per_row(i)
      DO nd=1,nbdirs
        xd(nd, index) = j*dxd(nd) + i*sheard(nd)
        yd(nd, index) = i*dyd(nd)
      END DO
      x(index) = x_start(i) + dx*j + i*shear
      y(index) = y0 + dy*i
      index = index + 1
    END DO
  END DO
  rotate_rad = rotate*3.1415926535/180.0
  DO nd=1,nbdirs
    rotate_radd(nd) = 3.1415926535*rotated(nd)/180.0
    turbinexd(nd, :) = COS(rotate_rad)*xd(nd, :) - rotate_radd(nd)*SIN(&
&     rotate_rad)*x - rotate_radd(nd)*COS(rotate_rad)*y - SIN(rotate_rad&
&     )*yd(nd, :)
    turbineyd(nd, :) = rotate_radd(nd)*COS(rotate_rad)*x + SIN(&
&     rotate_rad)*xd(nd, :) + COS(rotate_rad)*yd(nd, :) - rotate_radd(nd&
&     )*SIN(rotate_rad)*y
  END DO
  turbinex = COS(rotate_rad)*x - SIN(rotate_rad)*y
  turbiney = SIN(rotate_rad)*x + COS(rotate_rad)*y
END SUBROUTINE MAKEGRID_FORTRAN_DV

!  Differentiation of gaussianwake in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: loss
!   with respect to varying inputs: turbinexw turbineyw loss
SUBROUTINE GAUSSIANWAKE_DV(nturbines, turbinexw, turbinexwd, turbineyw, &
& turbineywd, turb_diam, relaxationfactor, loss, lossd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines
  REAL(dp), INTENT(IN) :: turb_diam, relaxationfactor
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinexw, turbineyw
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(IN) :: turbinexwd, &
& turbineywd
! out
  REAL(dp), DIMENSION(nturbines), INTENT(OUT) :: loss
  REAL(dp), DIMENSION(nbdirs, nturbines), INTENT(OUT) :: lossd
! local
  REAL(dp) :: ct, k, x, y, sigma, exponent, radical
  REAL(dp), DIMENSION(nbdirs) :: xd, yd, sigmad, exponentd, radicald
  REAL(dp), DIMENSION(nturbines) :: loss_array
  REAL(dp), DIMENSION(nbdirs, nturbines) :: loss_arrayd
  REAL(dp), PARAMETER :: pi=3.141592653589793_dp, e=2.718281828459045_dp&
& , tol=0.000001_dp
  INTEGER :: i, j
  INTRINSIC SQRT
  INTRINSIC SUM
  REAL :: result1
  REAL(dp) :: result10
  REAL(dp), DIMENSION(nbdirs) :: result10d
  REAL(dp) :: pwr1
  REAL(dp), DIMENSION(nbdirs) :: pwr1d
  REAL(dp), DIMENSION(nturbines) :: arg1
  REAL(dp), DIMENSION(nbdirs, nturbines) :: arg1d
  REAL(dp) :: arg2
  REAL(dp), DIMENSION(nbdirs) :: arg2d
  INTEGER :: nd
  INTEGER :: nbdirs
  ct = 4.0*1./3.*(1.0-1./3.)
  k = 0.0324555
  loss_arrayd(:, :) = 0.0_8
  DO i=1,nturbines
    DO j=1,nturbines
      DO nd=1,nbdirs
        xd(nd) = turbinexwd(nd, i) - turbinexwd(nd, j)
        yd(nd) = turbineywd(nd, i) - turbineywd(nd, j)
      END DO
      x = turbinexw(i) - turbinexw(j)
      y = turbineyw(i) - turbineyw(j)
      IF (x .GT. 0.0) THEN
        result1 = SQRT(8.0)
        sigma = k*x + turb_diam/result1
        exponent = -(0.5*(y/(relaxationfactor*sigma))**2)
        radical = 1. - ct/(8.*sigma**2/turb_diam**2)
        result10 = SQRT(radical)
        pwr1 = e**exponent
        DO nd=1,nbdirs
          sigmad(nd) = k*xd(nd)
          exponentd(nd) = -(0.5*2*y*(yd(nd)*relaxationfactor*sigma-y*&
&           relaxationfactor*sigmad(nd))/(relaxationfactor**3*sigma**3))
          radicald(nd) = ct*8.*2*sigma*sigmad(nd)/turb_diam**2/(8.*sigma&
&           **2/turb_diam**2)**2
          IF (radical .EQ. 0.0) THEN
            result10d(nd) = 0.0_8
          ELSE
            result10d(nd) = radicald(nd)/(2.0*SQRT(radical))
          END IF
          IF (e .GT. 0.0) THEN
            pwr1d(nd) = LOG(e)*e**exponent*exponentd(nd)
          ELSE
            pwr1d(nd) = 0.0
          END IF
          loss_arrayd(nd, j) = (1.-result10)*pwr1d(nd) - result10d(nd)*&
&           pwr1
        END DO
        loss_array(j) = (1.-result10)*pwr1
      ELSE
        DO nd=1,nbdirs
          loss_arrayd(nd, j) = 0.0_8
        END DO
        loss_array(j) = 0.0
      END IF
    END DO
    arg1(:) = loss_array**2
    arg2 = SUM(arg1(:))
    DO nd=1,nbdirs
      arg1d(nd, :) = 2*loss_array*loss_arrayd(nd, :)
      arg2d(nd) = SUM(arg1d(nd, :))
      IF (arg2 .EQ. 0.0) THEN
        lossd(nd, i) = 0.0_8
      ELSE
        lossd(nd, i) = arg2d(nd)/(2.0*SQRT(arg2))
      END IF
    END DO
    loss(i) = SQRT(arg2)
  END DO
END SUBROUTINE GAUSSIANWAKE_DV

!  Differentiation of calcaep_grid in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: aep
!   with respect to varying inputs: rotate dx dy shear
!   RW status of diff variables: rotate:in dx:in dy:in aep:out
!                shear:in
SUBROUTINE CALCAEP_GRID_DV(nturbines, ndirections, nrows, turbinez, &
& rotordiameter, winddirections, dx, dxd, dy, dyd, shear, sheard, rotate&
& , rotated, turbs_per_row, x_start, y0, windspeeds, windfrequencies, &
& shearexp, relaxationfactor, rated_ws, rated_power, cut_in_speed, &
& cut_out_speed, zref, z0, aep, aepd, nbdirs)

!  Hint: nbdirs should be the maximum number of differentiation directions
  IMPLICIT NONE
  INTRINSIC KIND
! define precision to be the standard for a double precision ! on local system
  INTEGER, PARAMETER :: dp=KIND(0.d0)
! in
  INTEGER, INTENT(IN) :: nturbines, ndirections, nrows
  REAL(dp), INTENT(IN) :: shearexp, relaxationfactor, rated_ws, &
& rated_power, cut_in_speed, cut_out_speed, zref, z0, dx, dy, shear, &
& rotate, y0
  REAL(dp), DIMENSION(nbdirs), INTENT(IN) :: dxd, dyd, sheard, &
& rotated
  REAL(dp), DIMENSION(nturbines), INTENT(IN) :: turbinez, rotordiameter
  REAL(dp), DIMENSION(ndirections), INTENT(IN) :: winddirections, &
& windspeeds, windfrequencies
  INTEGER, DIMENSION(nrows), INTENT(IN) :: turbs_per_row
  REAL(dp), DIMENSION(nrows), INTENT(IN) :: x_start
! out
  REAL(dp), INTENT(OUT) :: aep
  REAL(dp), DIMENSION(nbdirs), INTENT(OUT) :: aepd
! local
  REAL(dp), DIMENSION(nturbines) :: turbinex, turbiney
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbinexd, turbineyd
  REAL(dp), DIMENSION(ndirections) :: dir_powers
  REAL(dp), DIMENSION(nbdirs, ndirections) :: dir_powersd
  REAL(dp), DIMENSION(nturbines) :: turbinexw, turbineyw, vinf_floris, &
& wtvelocity, loss
  REAL(dp), DIMENSION(nbdirs, nturbines) :: turbinexwd, turbineywd, &
& wtvelocityd, lossd
  REAL(dp) :: hrs_per_year, pwrdir, vinf
  REAL(dp), DIMENSION(nbdirs) :: pwrdird
  INTEGER :: n, i
  INTRINSIC SUM
  REAL(dp), DIMENSION(ndirections) :: arg1
  REAL(dp), DIMENSION(nbdirs, ndirections) :: arg1d
  INTEGER :: nd
  INTEGER :: nbdirs
  CALL MAKEGRID_FORTRAN_DV(nrows, nturbines, dx, dxd, dy, dyd, shear, &
&                    sheard, rotate, rotated, turbs_per_row, x_start, y0&
&                    , turbinex, turbinexd, turbiney, turbineyd, nbdirs)
  dir_powersd(:, :) = 0.0_8
  lossd(:, :) = 0.0_8
  DO n=1,ndirections
    CALL WINDFRAME_DV(nturbines, winddirections(n), turbinex, turbinexd&
&               , turbiney, turbineyd, turbinexw, turbinexwd, turbineyw&
&               , turbineywd, nbdirs)
    CALL POWWIND(nturbines, windspeeds(n), turbinez, shearexp, zref, z0&
&          , vinf_floris)
    vinf = vinf_floris(1)
    CALL GAUSSIANWAKE_DV(nturbines, turbinexw, turbinexwd, turbineyw, &
&                  turbineywd, rotordiameter(1), relaxationfactor, loss&
&                  , lossd, nbdirs)
    DO nd=1,nbdirs
      wtvelocityd(nd, :) = -(vinf*lossd(nd, :))
    END DO
    wtvelocity = vinf*(1.0_dp-loss)
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
END SUBROUTINE CALCAEP_GRID_DV
