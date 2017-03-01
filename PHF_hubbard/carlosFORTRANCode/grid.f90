

module grid

  use constants

  implicit none
  public

! +----------------------------------------------------------------+
! |                                                                |
! | grid                                                           |
! |                                                                |
! |                                                                |
! | A collection of subroutines used to construct quadrature       |
! | integration grids.                                             |
! |                                                                |
! +----------------------------------------------------------------+

contains


subroutine gauleg (x1, x2, n, x, w)

! +----------------------------------------------------------------+
! |                                                                |
! | gauleg  --  CAJH, 11.2012                                      |
! |                                                                |
! |                                                                |
! | Subroutine to determine the abscissas and weights of the       |
! | Gauss-Legendre n-point quadrature formula.                     |
! |                                                                |
! | Given the lower and upper limits of integration x1 and x2,     |
! | and given n, this routine returns arrays x(1:n) and w(1:n)     |
! | of length n, containing the abscissas and weights of the       |
! | Gauss-Legendre n-point quadrature formula.                     |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The original code was obtained from                            |
! |                                                                |
! |   W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.        |
! |   Flannery, Numerical Recipes in Fortran 77: The Art of        |
! |   Scientific Computing. 2nd ed. Cambridge University Press,    |
! |   NY, 2005.                                                    |
! |                                                                |
! | CAJH - Modified it to use Fortran 90 constructs.               |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   x1 - lower limit of integration
  !   x2 - upper limit of integration
  !   n  - number of grid points
  !   x  - abscissas of Gauss-Legendre quadrature
  !   w  - weights of Gauss-Legendre quadrature

  integer, intent(in) :: n

  real(kind=dp), intent(in) :: x1, x2
  real(kind=dp), dimension(n), intent(inout) :: x, w


  ! other variables

  integer :: m, i, j

  real(kind=dp) :: xm, xl, p1, p2, p3, pp, z, z1


  ! some parameters

  real(kind=dp) :: pi, eps
  parameter ( pi  = 4.0e0_dp * atan(d1) )
  parameter ( eps = 1.0e-15_dp )


  ! Only half of the roots need to be found given that the polynomial
  ! is symmetric about the center of the interval.

  m = (n + 1)/2

  xm = 0.5e0_dp * (x2 + x1)
  xl = 0.5e0_dp * (x2 - x1)

  do i = 1, m

    ! Approximation to the i-th root.

    z = cos (pi * (real(i,dp)-0.25e0_dp) / (real(n,dp)+0.50e0_dp))

    ! Refine the polynomial root by Newton's method.

    z1 = z + d1
    pp = d1

    do while ( abs(z-z1) > eps )

      p1 = d1
      p2 = d0

      ! Loop over the recurrence relation to evaluate the Legendre
      ! polynomial at z.

      do j = 1, n
        p3 = p2
        p2 = p1
        p1 = ((d2*real(j,dp) - d1) * z * p2 - (real(j,dp)-d1) * p3) &
         & / real(j,dp)
      end do

      ! p1 is now the desired Legendre polynomial. We next compute
      ! pp, its derivative, by a standard relation involving also p2,
      ! the polynomial of one lower order.

      pp = real(n,dp) * (z*p1 - p2) / (z*z - d1)
      z1 = z

      z = z1 - p1 / pp  ! Newton's method...

    end do

    ! Scale the root to the desired interval, and put in its
    ! symmetric counterpart. Compute the weight and its symmetric
    ! counterpart.

    x(i)     = xm - xl*z
    x(n+1-i) = xm + xl*z

    w(i)     = d2 * xl / ((d1 - z*z)*pp*pp)
    w(n+1-i) = w(i)

  end do


  return
end subroutine gauleg



subroutine trpgrd (x1, x2, n, x, w)

! +----------------------------------------------------------------+
! |                                                                |
! | trpgrd  --  CAJH, 11.2012                                      |
! |                                                                |
! |                                                                |
! | Subroutine to determine the grid positions and weights for     |
! | a trapezoidal integration in the interval (x1,x2) using n      |
! | points.                                                        |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   x1 - lower limit of integration
  !   x2 - upper limit of integration
  !   n  - number of grid points
  !   x  - grid points to use in trapezoidal integration
  !   w  - weights to use in trapezoidal integration

  integer, intent(in) :: n

  real(kind=dp), intent(in) :: x1, x2
  real(kind=dp), dimension(n), intent(inout) :: x, w


  ! other variables

  integer :: i

  real(kind=dp) :: xl, spc, len


  ! determine initial position and spacing

  len = x2 - x1
  spc = len / real(n,dp)
  xl  = x1 + spc / d2

  ! fill array with positions and weights

  do i = 1, n
    x(i) = xl + real(i-1,dp)*spc
  end do

  w(1:n) = len / real(n,dp)


  return
end subroutine trpgrd


end module grid


