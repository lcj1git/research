

module linesrch

  use constants

  implicit none
  private
  save

  public :: mcsrch_par, mcsrch

! +----------------------------------------------------------------+
! |                                                                |
! | linesearch                                                     |
! |                                                                |
! |                                                                |
! | This module contains subroutines implementing the line-search  |
! | algorithm by More' and Thuente.                                |
! |                                                                |
! | Originally obtained from Jorge Nocedal's LBFGS implementation. |
! |                                                                |
! +----------------------------------------------------------------+

  ! global parameters

  integer :: n, maxfev

  real(kind=dp) :: ftol, gtol, xtol, stpmin, stpmax

  ! other variables (must be saved)

  integer :: infoc
  logical :: brackt, stage1

  real(kind=dp) :: finit, dginit, dgtest
  real(kind=dp) :: dgx, dgy, fx, fy
  real(kind=dp) :: stx, sty, stmin, stmax
  real(kind=dp) :: width, width1


contains


subroutine mcsrch_par (n1, mxfev, ftol1, gtol1, xtol1, stmin, stmax)

  integer, intent(in) :: n1, mxfev

  real(kind=dp), intent(in) :: ftol1, gtol1, xtol1, stmin, stmax

  n      = n1
  maxfev = mxfev

  ftol   = ftol1
  gtol   = gtol1
  xtol   = xtol1
  stpmin = stmin
  stpmax = stmax

  return
end subroutine mcsrch_par



subroutine mcsrch (x, f, g, s, stp, info, nfev, wa)

! +----------------------------------------------------------------+
! |                                                                |
! | mcsrch                                                         |
! |                                                                |
! |                                                                |
! | From Jorge Nocedal : A slight modification of the subroutine   |
! |   mcsrch of More' and Thuente. The changes are to allow        |
! |   reverse communication and do not affect the performance of   |
! |   the routine.                                                 |
! |                                                                |
! |                                                                |
! | The purpose of mcsrch is to find a step which satisfies a      |
! | sufficient decrease condition and a curvature condition.       |
! |                                                                |
! | At each stage the subroutine updates an interval of            |
! | uncertainty with endpoints stx and sty. The interval of        |
! | uncertainty is initially chosen so that it contains a          |
! | minimizer of the modified function                             |
! |                                                                |
! |   f(x+stp*s) - f(x) - ftol*stp*(gradf(x)! s).                  |
! |                                                                |
! | If a step is obtained for which the modified function has a    |
! | non-positive function value and non-negative derivative, then  |
! | the interval of uncertainty is chosen so that it contains a    |
! | minimizer of f(x+stp*s).                                       |
! |                                                                |
! | The algorithm is designed to find a step which satisfies the   |
! | sufficient decrease condition                                  |
! |                                                                |
! |   f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)! s),                 |
! |                                                                |
! | and the curvature condition                                    |
! |                                                                |
! |   abs(gradf(x+stp*s)! s) <= gtol*abs(gradf(x)! s).             |
! |                                                                |
! | If ftol is less than gtol and if, for example, the function is |
! | bounded below, then there is always a step which satisfies     |
! | both conditions. If no step can be found which satisfies both  |
! | conditions, then the algorithm usually stops when rounding     |
! | errors prevent further progress. In this case stp only         |
! | satisfies the sufficient decrease condition.                   |
! |                                                                |
! |                                                                |
! | - n: positive integer input variable set to the number of      |
! |     variables.                                                 |
! |                                                                |
! | - x: array of length n. On input it must contain the base      |
! |     point for the line search. On output it contains x+stp*s.  |
! |                                                                |
! | - f: a variable. On input it must contain the value of f at    |
! |     x. On output it contains the value of f at x+stp*s.        |
! |                                                                |
! | - g: array of length n. On input it must contain the gradient  |
! |     of f at x. On output it contains the gradient of f at      |
! |     x+stp*s.                                                   |
! |                                                                |
! | - s: input array of length n which specifies the search        |
! |     direction.                                                 |
! |                                                                |
! | - stp: a non-negative variable. On input stp contains an       |
! |     initial estimate of a satisfactory step. On output stp     |
! |     contains the final estimate.                               |
! |                                                                |
! | - ftol, gtol: non-negative input variables. Termination        |
! |     occurs when the sufficient decrease condition and the      |
! |     directional derivative condition are satisfied.            |
! |                                                                |
! | - xtol: non-negative input variable. Termination occurs when   |
! |     the relative width of the interval of uncertainty is at    |
! |     most xtol.                                                 |
! |                                                                |
! | - stpmin, stpmax: non-negative input variables which specify   |
! |     lower and upper bounds for the step.                       |
! |                                                                |
! | - maxfev: positive integer input variable. Termination occurs  |
! |     when the number of calls to fcn is at least maxfev by the  |
! |     end of an iteration.                                       |
! |                                                                |
! | - info: integer output variable set as follows:                |
! |                                                                |
! |     info =  0,  improper input parameters                      |
! |                                                                |
! |          = -1,  return is made to compute function and grad    |
! |                                                                |
! |          =  1,  the sufficient decrease condition and the      |
! |                 directional derivative condition hold          |
! |                                                                |
! |          =  2,  relative width of the interval of uncertainty  |
! |                 is at most xtol                                |
! |                                                                |
! |          =  3,  number of calls to fcn has reached maxfev      |
! |                                                                |
! |          =  4,  the step is at the lower bound stpmin          |
! |                                                                |
! |          =  5,  the step is at the upper bound stpmax          |
! |                                                                |
! |          =  6,  rounding errors prevent further progress.      |
! |                 There may not be a step which satisfies the    |
! |                 sufficient decrease and curvature conditions.  |
! |                 Tolerances may be too small.                   |
! |                                                                |
! | - nfev: integer output variable set to number of calls to fcn. |
! |                                                                |
! | - wa: work array of length n.                                  |
! |                                                                |
! | Argonne National Laboratory. MINPACK PROJECT. June 1983        |
! | Jorge J. More', David J. Thuente                               |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  integer, intent(inout) :: nfev
  integer, intent(inout) :: info

  real(kind=dp), intent(in) :: f
  real(kind=dp), intent(inout) :: stp
  real(kind=dp), dimension(n), intent(in) :: g, s
  real(kind=dp), dimension(n), intent(inout) :: x, wa


  ! local variables

  integer :: j

  real(kind=dp) :: ftest1
  real(kind=dp) :: dg, fm, dgm, fxm, dgxm, fym, dgym

  ! some parameters

  real(kind=dp) :: p5, p66, xtrapf
  parameter ( p5 = 0.5e0_dp, p66 = 0.66e0_dp, xtrapf = 4.0e0_dp )


  if ( info == -1 )  go to 100    ! returning...

  info  = 0
  infoc = 1

  ! Check the input parameters for errors.

  if ( n <= 0 .or. stp <= d0 .or. ftol < d0 .or. gtol < d0 &
     & .or. xtol < d0 .or. stpmin < d0 .or. stpmax < stpmin &
     & .or. maxfev <= 0 )  return


  ! Compute the initial gradient in the search direction and check
  ! that s is a descent direction.

  dginit = d0
  do j = 1, n
    dginit = dginit + g(j)*s(j)
  end do

  if ( dginit >= d0 ) then
    ! write (*, 15)
    return
  end if

  ! 15 format (/' The search direction is not a descent direction.')


  ! Initialize local variables.

  brackt = .false.
  stage1 = .true.

  nfev   = 0
  finit  = f
  dgtest = ftol*dginit
  width  = stpmax - stpmin
  width1 = width / p5

  wa(1:n) = x(1:n)


  ! The variables stx, fx, dgx contain the values of the step, function
  ! and directional derivative at the best step.

  ! The variables sty, fy, dgy contain the values of the step, function
  ! and derivative at the other endpoint of the interval of
  ! uncertainty.

  ! The variables stp, f, dg contain the values of the step, function
  ! and derivative at the current step.

  stx = d0
  fx  = finit
  dgx = dginit

  sty = d0
  fy  = finit
  dgy = dginit


  ! Start of iteration.

  10  continue

    ! Set the minimum and maximum steps to correspond to the present
    ! interval of uncertainty.

    if ( brackt ) then
      stmin = min (stx,sty)
      stmax = max (stx,sty)
    else
      stmin = stx
      stmax = stp + xtrapf*(stp-stx)
    end if

    ! Force the step to be within the bounds stpmax and stpmin.

    stp = max (stp,stpmin)
    stp = min (stp,stpmax)

    ! If an unusual termination is to occur then let stp be the lowest
    ! point obtained so far.

    if ( ( brackt .and. ( stp <= stmin .or. stp >= stmax ) ) &
       & .or. nfev >= maxfev-1 .or. infoc == 0 &
       & .or. ( brackt .and. stmax-stmin <= xtol*stmax ) )  stp = stx


    ! Evaluate the function and gradient at stp and compute the
    ! directional derivative.
    ! (* Return to the main program to obtain f and g. *)

    x(1:n) = wa(1:n) + stp*s(1:n)

      info = -1
      return

    100  continue    ! return here

    info = 0
    nfev = nfev + 1

    dg = d0
    do j = 1, n
      dg = dg + g(j)*s(j)
    end do

    ftest1 = finit + stp*dgtest


    ! Test for convergence.

    if ( ( brackt .and. ( stp <= stmin .or. stp >= stmax ) ) &
       & .or. infoc == 0 )  info = 6

    if ( stp == stpmax .and. &
       & f <= ftest1 .and. dg <= dgtest )  info = 5

    if ( stp == stpmin .and. &
       & ( f > ftest1 .or. dg >= dgtest ) )  info = 4

    if ( nfev >= maxfev )  info = 3

    if ( brackt .and. stmax-stmin <= xtol*stmax )  info = 2

    if ( f <= ftest1 .and. abs(dg) <= gtol*(-dginit) )  info = 1

    ! Check for termination.

    if ( info /= 0 )  return


    ! In the first stage we seek a step for which the modified
    ! function has a non-positive value and non-negative derivative.

    if ( stage1 .and. f <= ftest1 .and. &
       & dg >= min(ftol,gtol)*dginit )  stage1 = .false.


    ! A modified function is used to predict the step only if we have
    ! not obtained a step for which the modified function has a
    ! non-positive function value and non-negative derivative, and if
    ! a lower function value has been obtained but the decrease is not
    ! sufficient.

    if ( stage1 .and. f <= fx .and. f > ftest1 ) then

      ! Define the modified function and derivative values.

      fm  = f  - stp*dgtest
      fxm = fx - stx*dgtest
      fym = fy - sty*dgtest

      dgm  = dg  - dgtest
      dgxm = dgx - dgtest
      dgym = dgy - dgtest

      ! Call mcstep to update the interval of uncertainty and to
      ! compute the new step.

      call mcstep (stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, &
           & brackt, stmin, stmax, infoc)

      ! Reset the function and gradient values for f.

      fx = fxm + stx*dgtest
      fy = fym + sty*dgtest

      dgx = dgxm + dgtest
      dgy = dgym + dgtest

    else

      ! Call mcstep to update the interval of uncertainty and to
      ! compute the new step.

      call mcstep (stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, &
           & stmin, stmax, infoc)

    end if

    ! Force a sufficient decrease in the size of the interval of
    ! uncertainty.

    if ( brackt ) then
      if ( abs(sty-stx) >= p66*width1 ) then
        stp = stx + p5*(sty-stx)
      end if

      width1 = width
      width  = abs(sty-stx)
    end if

  go to 10

  ! End of iteration.


  return
end subroutine mcsrch



subroutine mcstep (stx, fx, dx, sty, fy, dy, stp, fp, d1p, brackt, &
     & stpmin, stpmax, info)

! +----------------------------------------------------------------+
! |                                                                |
! | mcstep                                                         |
! |                                                                |
! |                                                                |
! | The purpose of mcstep is to compute a safeguarded step for a   |
! | line search and to update an interval of uncertainty for a     |
! | minimizer of the function.                                     |
! |                                                                |
! | The parameter stx contains the step with the least function    |
! | value. The parameter stp contains the current step. It is      |
! | assumed that the derivative at stx is negative in the          |
! | direction of the step. If brackt is set to true then a         |
! | minimizer has been bracketed in an interval of uncertainty     |
! | with endpoints stx and sty.                                    |
! |                                                                |
! |                                                                |
! | - stx, fx, dx: variables which specify the step, the function, |
! |     and the derivative at the best step obtained so far. The   |
! |     derivative must be negative in the direction of the step.  |
! |     That is, dx and stp-stx must have opposite signs. On       |
! |     output these parameters are updated appropriately.         |
! |                                                                |
! | - sty, fy, dy: variables which specify the step, the function, |
! |     and the derivative at the other endpoint of the interval   |
! |     of uncertainty. On output these parameters are updated     |
! |     appropriately.                                             |
! |                                                                |
! | - stp, fp, dp: variables which specify the step, the function, |
! |     and the derivative at the current step. If brackt is set   |
! |     to true, then on input stp must be between stx and sty.    |
! |     On output stp is set to the new step.                      |
! |                                                                |
! | - brackt: logical variable which specifies if a minimizer has  |
! |     been bracketed. If the minimizer has not been bracketed    |
! |     then on input brack must be set to false. If the minimizer |
! |     is bracketed then on output brackt is set to true.         |
! |                                                                |
! | - stpmin, stpmax: input variables which specify lower and      |
! |     upper bounds for the step.                                 |
! |                                                                |
! | - info: integer output variable set as follows:                |
! |     info = 1 .. 4  =>  the step has been computed according    |
! |                        to one of the cases below               |
! |          = 0       =>  improper input parameters               |
! |                                                                |
! | Argonne National Laboratory. MINPACK PROJECT. June 1983        |
! | Jorge J. More', David J. Thuente                               |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  integer, intent(out) :: info

  real(kind=dp), intent(in) :: stpmin, stpmax

  real(kind=dp), intent(in) :: fp
  real(kind=dp), intent(inout) :: stx, fx, dx
  real(kind=dp), intent(inout) :: sty, fy, dy
  real(kind=dp), intent(inout) :: stp, d1p

  logical, intent(inout) :: brackt


  ! other variables

  real(kind=dp) :: sgnd, theta, gamma, s, p, q, r
  real(kind=dp) :: stpc, stpq, stpf

  logical :: bound


  info = 0

  ! Check the input parameters for errors.

  if ( ( brackt .and. ( stp <= min (stx,sty) &
              &  .or.  stp >= max (stx,sty) ) ) .or. &
      & dx*(stp-stx) >= d0 .or. stpmax < stpmin ) return

  ! Determine if the derivatives have opposite sign.

  sgnd = d1p*(dx/abs(dx))


  ! First case. A higher function value. The minimum is bracketed. If
  ! the cubic step is closer to stx than the quadratic step, the cubic
  ! step is taken, else the average of the cubic and quadratic steps
  ! is taken.

  if ( fp > fx ) then

    info  = 1
    bound = .true.

    theta = d3 * (fx-fp)/(stp-stx) + dx + d1p
    s     = max (abs(theta), abs(dx), abs(d1p))
    gamma = s * sqrt((theta/s)**2 - (dx/s)*(d1p/s))

    if ( stp < stx ) gamma = -gamma

    p = (gamma - dx) + theta
    q = ((gamma - dx) + gamma) + d1p
    r = p / q

    stpc = stx + r*(stp-stx)
    stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/d2)*(stp-stx)

    if ( abs(stpc-stx) < abs(stpq-stx) ) then
      stpf = stpc
    else
      stpf = stpc + (stpq - stpc)/d2
    end if

    brackt = .true.

  ! Second case. A lower function value and derivatives of opposite
  ! sign. The minimum is bracketed. If the cubic step is closer to
  ! stx than the quadratic (secant) step, the cubic step is taken,
  ! else the quadratic step is taken.

  else if ( sgnd < d0 ) then

    info  = 2
    bound = .false.

    theta = d3 * (fx-fp)/(stp-stx) + dx + d1p
    s     = max (abs(theta), abs(dx), abs(d1p))
    gamma = s * sqrt((theta/s)**2 - (dx/s)*(d1p/s))

    if ( stp > stx ) gamma = -gamma

    p = (gamma - d1p) + theta
    q = ((gamma - d1p) + gamma) + dx
    r = p / q

    stpc = stp + r*(stx-stp)
    stpq = stp + (d1p/(d1p-dx))*(stx-stp)

    if ( abs(stpc-stp) > abs(stpq-stp) ) then
      stpf = stpc
    else
      stpf = stpq
    end if

    brackt = .true.

  ! Third case. A lower function value, derivatives of the same sign,
  ! and the magnitude of the derivative decreases. The cubic step is
  ! only used if the cubic tends to infinity in the direction of the
  ! step or if the minimum of the cubic is beyond stp. Otherwise the
  ! cubic step is defined to be either stpmin or stpmax. The quadratic
  ! (secant) step is also computed and if the minimum is bracketed
  ! then the step closest to stx is taken, else the step farthest away
  ! is taken.

  else if ( abs(d1p) < abs(dx) ) then

    info  = 3
    bound = .true.

    theta = d3 * (fx-fp)/(stp-stx) + dx + d1p
    s     = max (abs(theta), abs(dx), abs(d1p))
    gamma = s * sqrt(max (d0, (theta/s)**2 - (dx/s)*(d1p/s)))

      ! the case gamma = 0 only arises if the cubic does not tend to
      ! infinity in the direction of the step

    if ( stp > stx ) gamma = -gamma

    p = (gamma - d1p) + theta
    q = (gamma + (dx-d1p)) + gamma
    r = p / q

    if ( r < d0 .and. gamma /= d0 ) then
      stpc = stp + r*(stx-stp)
    else if ( stp > stx ) then
      stpc = stpmax
    else
      stpc = stpmin
    end if

    stpq = stp + (d1p/(d1p-dx))*(stx-stp)

    if ( brackt ) then
      if ( abs(stp-stpc) < abs(stp-stpq) ) then
        stpf = stpc
      else
        stpf = stpq
      end if

    else
      if ( abs(stp-stpc) > abs(stp-stpq) ) then
        stpf = stpc
      else
        stpf = stpq
      end if
    end if

  ! Fourth case. A lower function value, derivatives of the same sign,
  ! and the magnitude of the derivative does not decrease. If the
  ! minimum is not bracketed, the step is either stpmin or stpmax,
  ! else the cubic step is taken.

  else

    info  = 4
    bound = .false.

    if ( brackt ) then

      theta = d3 * (fp-fy)/(sty-stp) + dy + d1p
      s     = max (abs(theta), abs(dy), abs(d1p))
      gamma = s * sqrt((theta/s)**2 - (dy/s)*(d1p/s))

      if ( stp > sty ) gamma = -gamma

      p = (gamma - d1p) + theta
      q = ((gamma - d1p) + gamma) + dy
      r = p / q

      stpc = stp + r*(sty-stp)
      stpf = stpc

    else if ( stp > stx ) then
      stpf = stpmax
    else
      stpf = stpmin
    end if

  end if


  ! Update the interval of uncertainty. This update does not depend on
  ! the new step or the case analysis above.

  if ( fp > fx ) then
    sty = stp
    fy  = fp
    dy  = d1p

  else
    if ( sgnd < 0 ) then
      sty = stx
      fy  = fx
      dy  = dx
    end if

    stx = stp
    fx  = fp
    dx  = d1p
  end if


  ! Compute the new step and safeguard it.

  stpf = min (stpmax, stpf)
  stpf = max (stpmin, stpf)
  stp  = stpf

  if ( brackt .and. bound ) then
    if ( sty > stx ) then
      stp = min (stx+0.66e0_dp*(sty-stx), stp)
    else
      stp = max (stx+0.66e0_dp*(sty-stx), stp)
    end if
  end if


  return
end subroutine mcstep


end module linesrch


