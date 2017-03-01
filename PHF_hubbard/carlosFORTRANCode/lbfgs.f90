

module lbfgs

  use constants
  use util, only : error_alloc
  use linesrch

  implicit none
  private
  save

  public :: lbfgs_par, lbfgs_setup, shutdown_lbfgs, lbfgs_drv

! +----------------------------------------------------------------+
! |                                                                |
! | lbfgs                                                          |
! |                                                                |
! |                                                                |
! | This module contains subroutines implementing the LBFGS        |
! | algorithm for large-scale optimizations.                       |
! |                                                                |
! | Original implementation by Jorge Nocedal, July 1990.           |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | CAJH:  Small modifications done so far:                        |
! |                                                                |
! | - Changed to fortran 90 and the use of modules.                |
! |                                                                |
! | - The subroutine lbfgs_setup must be called before any         |
! |   calls to lbfgs_drv.                                          |
! |                                                                |
! | - The subroutine shutdown_lbfgs can be used to deallocate      |
! |   memory of working arrays once optimization is complete.      |
! |                                                                |
! | - Global parameters can be modifed with lbfgs_par.             |
! |                                                                |
! | - Modified iprint to be of dimension (3):                      |
! |     iprint(3) > 0  =>  unit mp must be opened                  |
! |                        and closed when printing output         |
! |                                                                |
! | - Changed convergence criterion:                               |
! |                                                                |
! |   It used to converge when                                     |
! |     ||G|| < eps * max(1,||X||),                                |
! |   This has changed to                                          |
! |     ||G|| < eps.                                               |
! |                                                                |
! +----------------------------------------------------------------+


  ! output variables

  integer :: lp, mp

  character(len=20) :: outfile

  ! global variables

  integer :: n, m, maxfev
  logical :: diagco

  integer, dimension(3) :: iprint

  real(kind=dp) :: eps, ftol, gtol, xtol, stpmin, stpmax

  real(kind=dp), dimension(:), allocatable :: wrk_scr, wrk_alp
  real(kind=dp), dimension(:), allocatable :: wrk_rho
  real(kind=dp), dimension(:,:), allocatable :: wrk_spt, wrk_ypt

  ! other variables (must be saved)

  integer :: nfev, info, point, bound, npt
  logical :: finish

  real(kind=dp) :: stp


contains


subroutine lbfgs_par (ftol1, gtol1, xtol1, stmin, stmax, mxfev)

! +----------------------------------------------------------------+
! |                                                                |
! | lbfgs_par  --  CAJH, 11.2012                                   |
! |                                                                |
! |                                                                |
! | Subroutine that may be used to change some global parameters   |
! | in the lbfgs module.                                           |
! |                                                                |
! +----------------------------------------------------------------+

  integer, intent(in) :: mxfev

  real(kind=dp), intent(in) :: ftol1, xtol1, stmin, stmax, gtol1

  ftol   = ftol1
  gtol   = gtol1
  xtol   = xtol1
  stpmin = stmin
  stpmax = stmax

  maxfev = mxfev

  call mcsrch_par (n, maxfev, ftol, gtol, xtol, stpmin, stpmax)


  return
end subroutine lbfgs_par



subroutine lbfgs_setup (iout, outfl, n1, m1, ldiag, iprnt, eps1)

! +----------------------------------------------------------------+
! |                                                                |
! | lbfgs_setup  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Setup some global variables for lbfgs.                         |
! | Allocate memory for scratch arrays used by lbfgs.              |
! | This subroutine must be called before any call to lbfgs_drv.   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  integer, intent(in) :: iout, n1, m1
  integer, dimension(3), intent(in) :: iprnt

  logical, intent(in) :: ldiag
  character(len=*), intent(in) :: outfl

  real(kind=dp), intent(in) :: eps1

  ! other variables

  integer :: istatus


  outfile = outfl

  mp = iout
  lp = 0
  n  = n1
  m  = m1

  eps    = eps1
  diagco = ldiag

  iprint(1:3) = iprnt(1:3)

  ! parameters for line search routine

  ftol   = 1.0e-4_dp
  gtol   = 0.9e0_dp
  xtol   = 1.0e-16_dp
  stpmin = 1.0e-20_dp
  stpmax = 1.0e+20_dp
  maxfev = 20

  call mcsrch_par (n, maxfev, ftol, gtol, xtol, stpmin, stpmax)

  ! allocate working arrays

  allocate (wrk_scr(n), wrk_alp(m), wrk_rho(m), &
          & wrk_spt(n,m), wrk_ypt(n,m), stat=istatus)
  if ( istatus /= 0 ) call error_alloc(1, 'wrk', 'lbfgs_setup')


  return
end subroutine lbfgs_setup



subroutine shutdown_lbfgs

  integer :: istatus


  deallocate (wrk_scr, wrk_alp, wrk_rho, wrk_spt, wrk_ypt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc(2, 'wrk', 'shutdown_lbfgs')

  return
end subroutine shutdown_lbfgs



subroutine lbfgs_drv (x, f, g, diag, iter, nfun, iflag)

! +----------------------------------------------------------------+
! |                                                                |
! | lbfgs_drv                                                      |
! |                                                                |
! |                                                                |
! | == Limited-memory BFGS method for large-scale optimization ==  |
! |    Original implementation by  Jorge Nocedal, July 1990.       |
! |                                                                |
! |                                                                |
! | This subroutine solves the unconstrained minimization problem  |
! |                                                                |
! |   min F(x),   x = (x1,x2,...,xN),                              |
! |                                                                |
! | using the limited-memory BFGS method. The routine is           |
! | especially effective on problems involving a large number of   |
! | variables. In a typical iteration of this method an            |
! | approximation Hk to the inverse of the Hessian is obtained by  |
! | applying M BFGS updates to a diagonal matrix Hk0, using        |
! | information from the previous M steps. The user specifies the  |
! | number M, which determines the amount of storage required by   |
! | the routine. The user may also provide the diagonal matrices   |
! | Hk0 if not satisfied with the default choice. The algorithm    |
! | is described on                                                |
! |                                                                |
! |   D. Liu and J. Nocedal, 'On the limited memory BFGS method    |
! |     for large scale optimization', Mathematical Programming B, |
! |     45, 503-528 (1989).                                        |
! |                                                                |
! | The user is required to calculate the function value F and its |
! | gradient G. In order to allow the user complete control over   |
! | these computations, reverse communication is used. The routine |
! | must be called repeatedly under the control of the parameter   |
! | IFLAG.                                                         |
! |                                                                |
! | The steplength is determined at each iteration by means of the |
! | line search routine mcsrch, which is a slight modification of  |
! | the routine csrch written by More' and Thuente.                |
! |                                                                |
! |                                                                |
! | - n: integer variable that must be set by the user to the      |
! |     number of variables. It is not altered by the routine.     |
! |     Restriction: N > 0.                                        |
! |                                                                |
! | - m: integer variable that must be set by the user to the      |
! |     number of corrections used in the BFGS update. It is not   |
! |     altered by the routine. Values of M less than 3 are not    |
! |     recommended; large values of M will result in excessive    |
! |     computing time. 3 <= M <= 7 is recommended.                |
! |     Restriction: M > 0.                                        |
! |                                                                |
! | - x: double precision array of length N. On initial entry it   |
! |     must be set by the user to the values of the initial       |
! |     estimate of the solution vector. On exit with IFLAG = 0,   |
! |     it contains the values of the variables at the best point  |
! |     found (usually a solution).                                |
! |                                                                |
! | - f: double precision variable. Before initial entry and on    |
! |     a re-entry with IFLAG = 1, it must be set by the user to   |
! |     contain the value of the function F at the point X.        |
! |                                                                |
! | - g: double precision array of length N. Before initial entry  |
! |     and on a re-entry with IFLAG = 1, it must be set by the    |
! |     user to contain the components of the gradient G at the    |
! |     point X.                                                   |
! |                                                                |
! | - diagco: logical variable that must be set to .TRUE. if the   |
! |     user wishes to provide the diagonal matrix Hk0 at each     |
! |     iteration. Otherwise it should be set to .FALSE., in which |
! |     case LBFGS will use a default value described below. If    |
! |     DIAGCO is set to .TRUE. the routine will return at each    |
! |     iteration of the algorithm with IFLAG = 2, and the         |
! |     diagonal matrix Hk0 must be provided in the array DIAG.    |
! |                                                                |
! | - diag: double precision array of length N. If DIAGCO=.TRUE.,  |
! |     then on initial entry or on re-entry with IFLAG = 2, DIAG  |
! |     must be set by the user to contain the values of the       |
! |     diagonal matrix Hk0.                                       |
! |     Restriction: all elements of DIAG must be positive.        |
! |                                                                |
! | - iprint: integer array of length 2 which must be set by the   |
! |     user.                                                      |
! |                                                                |
! |     iprint(1) specifies the frequency of the output:           |
! |       iprint(1) < 0 : no output is generated.                  |
! |                 = 0 : output only at first and last iteration. |
! |                 > 0 : output every iprint(1) iterations.       |
! |                                                                |
! |     iprint(2) specifies the type of output generated:          |
! |       iprint(2) = 0 : iteration count, number of function      |
! |                       evaluations, function value, norm of     |
! |                       the gradient, step length                |
! |                 = 1 : same as iprint(2) = 0, plus vector of    |
! |                       variables and gradient vector at the     |
! |                       initial point                            |
! |                 = 2 : same as iprint(2) = 1, plus vector of    |
! |                       variables                                |
! |                 = 3 : same as iprint(2) = 2, plus gradient     |
! |                       vector                                   |
! |                                                                |
! | - eps: positive double precision variable that must be set     |
! |     by the user, and determines the accuracy with which the    |
! |     solution is to be found. The subroutine terminates when    |
! |                                                                |
! |       ||G|| < eps * max(1,||X||),                              |
! |                                                                |
! |     where ||.|| denotes the Euclidean norm.                    |
! |                                                                |
! | - xtol: double precision variable that must be set by the user |
! |     to an estimate of the machine precision (e.g. 1.0e-16).    |
! |     The line search routine will terminate if the relative     |
! |     width of the interval of uncertainty is less than XTOL.    |
! |                                                                |
! | - w: double precision array of length N*(2M+1)+2M used as      |
! |     workspace for LBFGS. This array must not be altered by     |
! |     the user.                                                  |
! |                                                                |
! | - iflag: integer variable that must be set to 0 on initial     |
! |     entry to the subroutine. A return with IFLAG < 0 indicates |
! |     an error, and IFLAG = 0 indicates that the routine has     |
! |     terminated without detecting errors. On a return with      |
! |     IFLAG = 1, the user must evaluate the function F and       |
! |     gradient G. On a return with IFLAG = 2, the user must      |
! |     provide the diagonal matrix Hk0.                           |
! |                                                                |
! |     The following negative values of IFLAG, detecting an       |
! |     error, are possible:                                       |
! |                                                                |
! |       iflag = -1,  The line search routine mcsrch failed.      |
! |                    The parameter INFO provides more detailed   |
! |                    information (see documentation of mcsrch).  |
! |                                                                |
! |         info = 0,  improper input parameters                   |
! |              = 2,  relative width of the interval of           |
! |                    uncertainty is at most xtol                 |
! |              = 3,  more than 20 function evaluations were      |
! |                    required at the present iteration           |
! |              = 4,  the step is too small                       |
! |              = 5,  the step is too large                       |
! |              = 6,  rounding errors prevent further progress    |
! |                    There may not be a step which satisfies the |
! |                    sufficient decrease and curvature           |
! |                    conditions. Tolerances may be too small.    |
! |                                                                |
! |       iflag = -2,  The i-th diagonal element of the diagonal   |
! |                    inverse Hessian approximation, given in     |
! |                    DIAG, is not positive.                      |
! |                                                                |
! |       iflag = -3,  Improper input parameters for LBFGS (N or   |
! |                    M are not positive).                        |
! |                                                                |
! |                                                                |
! | ON SOME REFERENCE VARIABLES:                                   |
! |                                                                |
! | - mp: integer variable with default value 6. It is used as the |
! |     unit number for the printing of the monitoring info        |
! |     controlled by IPRINT.                                      |
! |                                                                |
! | - lp: integer variable with default value 6. It is used as the |
! |     unit number for the printing of error messages. This       |
! |     printing may be suppressed by setting LP to a non-positive |
! |     value.                                                     |
! |                                                                |
! | - gtol: double precision variable with default value 0.9,      |
! |     which controls the accuracy of the line search routine     |
! |     mcsrch. If the function and gradient evaluations are       |
! |     inexpensive with respect to the cost of the iteration      |
! |     (which is sometimes the case when solving very large       |
! |     problems), it may be advantageous to set GTOL to a small   |
! |     value. A typical small value is 0.1.                       |
! |     Restriction: gtol should be greater than 1.0e-4.           |
! |                                                                |
! | - stpmin, stpmax: non-negative double precision variables      |
! |     which specify lower and upper bounds for the step in the   |
! |     line search. Their default values are 1.0e-20 and 1.0e+20, |
! |     respectively. These values need not be modified unless the |
! |     exponents are too large for the machine being used, or     |
! |     unless the problem is extremely badly scaled (in which     |
! |     case the exponents should be increased).                   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  integer, intent(inout) :: nfun, iter
  integer, intent(inout) :: iflag

  real(kind=dp), intent(in) :: f

  real(kind=dp), dimension(n), intent(inout) :: x, g, diag


  ! local variables

  integer :: i, cp

  real(kind=dp) :: sq, yr, ys, yy, beta
  real(kind=dp) :: gnorm  ! xnorm
  real(kind=dp) :: stp1

  ! external functions

  real(kind=dp), external :: ddot


  ! format statements

  1001 format (/' iflag = -1', /' Line search failed. Error return of ', &
              & 'line search: info = ', I2, '.')
  1002 format (/' iflag = -2', /' The ', I4, '-th diagonal element ', &
              & 'of the inverse Hessian approximation is not positive.')
  1003 format (/' iflag = -3', /' Improper input parameters: ', &
              & 'n or m are not positive.')
  1004 format (/' gtol is <= 1.0d-4.', /' gtol reset to 0.9d0.')


  if ( iflag == 1 ) then         ! returning with iflag = 1
    go to 10
  else if ( iflag == 2 ) then    ! returning with iflag = 2
    go to 20
  end if


  ! Error checking.

  if ( iflag /= 0 ) return

  if ( n <= 0 .or. m <= 0 ) then
    iflag = -3
    if (lp > 0)  write (lp, 1003)
    return
  end if

  if ( gtol <= 1.0e-4_dp ) then
    if (lp > 0)  write (lp, 1004)
    gtol = 0.9e0_dp
  end if

  ! Initialize.

  iter   = 0
  nfun   = 1
  point  = 0
  finish = .false.

  if ( diagco ) then
    do i = 1, n
      if ( diag(i) <= d0 ) then
        iflag = -2
        if (lp > 0)  write (lp, 1002)  i
        return
      end if
    end do

  else
    diag(1:n) = d1
  end if 


  ! Working arrays:

  !   wrk_scr(1:n)  -  scratch array
  !   wrk_rho(1:m)  -  stores scalars rho
  !   wrk_alp(1:m)  -  stores numbers alpha used in the formula that
  !                    computes H . G
  !   wrk_stp(1:n,1:m)  -  stores the last M search steps
  !   wrk_ytp(1:n,1:m)  -  stores the last M gradient differences

  ! The search steps and gradient differences are stored in a circular
  ! order controlled by the parameter point.

  wrk_spt(1:n,1) = -g(1:n) * diag(1:n)

  gnorm = sqrt (ddot (n, g(1:n), 1, g(1:n), 1))
  stp1  = d1 / gnorm

    if ( gnorm <= eps )  finish = .true.

  if ( iprint(1) >= 0 ) then
    call lb_prnt (gnorm, x(1:n), f, g(1:n), stp1, iter, nfun)
  end if
    
    if ( finish ) then
      iflag = 0
      return
    end if


  ! Start of iteration loop.

  100  continue

    iter  = iter + 1
    info  = 0
    bound = iter - 1
    
    if ( iter == 1 )  go to 11
    if ( iter > m )   bound = m
    
    ys = ddot (n, wrk_ypt(1:n,npt+1), 1, wrk_spt(1:n,npt+1), 1)
    
      cp = point
      if ( point == 0 ) cp = m
    
      wrk_rho(cp) = d1 / ys
      wrk_scr(1:n) = -g(1:n)

    
    if ( .not. diagco ) then
      yy = ddot (n, wrk_ypt(1:n,npt+1), 1, wrk_ypt(1:n,npt+1), 1)
      diag(1:n) = ys / yy
    
    else
      iflag = 2   ! need the inverse Hessian
      return
    end if
    
    20  continue    ! return here with iflag = 2
    
    ! If returning with inverse Hessian, make sure that all elements are
    ! positive.
    
    if ( diagco ) then
      do i = 1, n
        if ( diag(i) <= d0 ) then
          iflag = -2
          if (lp > 0)  write (lp, 1002)  i
          return
        end if
      end do
    end if
    
    ! Compute -H*G using the formula given in:
    
    !   J. Nocedal, 'Updating quasi-Newton matrices with limited
    !     storage', Mathematics of Computation, Vol. 24, No. 151,
    !     pp. 773-782 (1980).
    
    !cp = point
    !if ( point == 0 ) cp = m
    
    !wrk_rho(cp) = d1 / ys
    !wrk_scr(1:n) = -g(1:n)
    
    cp = point
    
    do i = 1, bound
      cp = cp-1
      if ( cp == -1 ) cp = m-1
    
      sq = ddot (n, wrk_spt(1:n,cp+1), 1, wrk_scr(1:n), 1)
    
      wrk_alp(cp+1) = wrk_rho(cp+1)*sq

      call daxpy (n, -wrk_alp(cp+1), wrk_ypt(1:n,cp+1), 1, wrk_scr(1:n), 1)
    end do
    
    wrk_scr(1:n) = diag(1:n) * wrk_scr(1:n)
    
    do i = 1, bound
      yr = ddot (n, wrk_ypt(1:n,cp+1), 1, wrk_scr(1:n), 1)
    
      beta = yr * wrk_rho(cp+1)
      beta = wrk_alp(cp+1) - beta
    
      call daxpy (n, beta, wrk_spt(1:n,cp+1), 1, wrk_scr(1:n), 1)
    
      cp = cp + 1
      if ( cp == m ) cp = 0
    end do
    
    ! Store the new search direction.
    
    wrk_spt(1:n,point+1) = wrk_scr(1:n)
    
    
    ! Obtain the one-dimensional minimizer of the function by using the
    ! line search routine mcsrch.
    
    11  continue
    
    nfev = 0
    stp  = d1
    if ( iter == 1 )  stp = stp1
    
    wrk_scr(1:n) = g(1:n)
    
    
    10  continue    ! return here with iflag = 1
    
    call mcsrch (x(1:n), f, g(1:n), wrk_spt(1:n,point+1), stp, &
         & info, nfev, diag(1:n))
    
    if ( info == -1 ) then
      iflag = 1
      return
    end if
    
    if ( info /=  1 ) then
      iflag = -1
      if (lp > 0)  write (lp, 1001)  info
      return
    end if
    
    nfun = nfun + nfev
    
    
    ! Compute the new step and gradient change.

    npt = point
    
    wrk_spt(1:n,npt+1) = stp * wrk_spt(1:n,npt+1)
    wrk_ypt(1:n,npt+1) = g(1:n)-wrk_scr(1:n)
    
    point = point + 1
    if ( point == m ) point = 0
    
    
    ! Termination test.
    
    gnorm = sqrt (ddot (n, g(1:n), 1, g(1:n), 1))
    !xnorm = sqrt (ddot (n, x(1:n), 1, x(1:n), 1))
    !xnorm = max (d1, xnorm)
    
    !if ( gnorm / xnorm <= eps )  finish = .true.
    if ( gnorm <= eps )  finish = .true.
    
    if ( iprint(1) >= 0 ) then
      call lb_prnt (gnorm, x(1:n), f, g(1:n), stp, iter, nfun)
    end if
    
    if ( finish ) then
      iflag = 0
      return
    end if


  go to 100

  ! End of iteration.


  return
end subroutine lbfgs_drv



subroutine lb_prnt (gnorm, x, f, g, stp, iter, nfun)

  ! input variables

  integer, intent(in) :: iter, nfun

  real(kind=dp), intent(in) :: gnorm, f, stp
  real(kind=dp), dimension(n), intent(in) :: x, g

  ! other variables

  integer :: i


  if ( iter == 0 ) then

    if ( iprint(3) > 0 ) then
      open (unit = mp, file = outfile, status = 'old', &
          & form = 'formatted', position = 'append')
    end if

    write (mp, 10)
    write (mp, 20) f, gnorm

    if ( iprint(2) >= 1 ) then
      write (mp, 30)
      write (mp, 50) (x(i), i = 1, n)
      write (mp, 40)
      write (mp, 50) (g(i), i = 1, n)
    end if

  else

    if ( ( iprint(1) == 0 ) .and. &
       & ( iter /= 1 .and. .not. finish ) ) return

    if ( iprint(3) > 0 ) then
      open (unit = mp, file = outfile, status = 'old', &
          & form = 'formatted', position = 'append')
    end if

    if ( iprint(1) /= 0 ) then
      if ( mod(iter-1,iprint(1)) == 0 .or. finish ) then
        write (mp, 60) iter, nfun, f, gnorm, stp
      end if
    else
      write (mp, 60) iter, nfun, f, gnorm, stp
    end if

    if ( iprint(2) >= 2 ) then
      if ( finish ) then
        write (mp, 31)
      else
        write (mp, 30)
      end if

      write (mp, 50) (x(i), i = 1, n)

      if ( iprint(2) >= 3 ) then
        write (mp, 40)
        write (mp, 50) (g(i), i = 1, n)
      end if
    end if

    if ( finish ) write (mp, 80)
  end if

  if ( iprint(3) > 0 ) then
    close (unit = mp)
  end if


  ! format statements

  10 format (/' initial values:')
  20 format (/' f     = ', F20.12, /' gnorm = ', F20.12)
  30 format (//' vector x:'/)
  31 format (//' final vector x:'/)
  40 format (//' gradient vector g:'/)
  50 format (1P, 6(2X, E14.6))
  60 format (/' iter  = ', I20, /' nfunc = ', I20, &
           & /' f     = ', F20.12, /' gnorm = ', F20.12, &
           & /' step  = ', 1P, E20.6)
  80 format (/' The minimization terminated without errors.')


  return
end subroutine lb_prnt


end module lbfgs


