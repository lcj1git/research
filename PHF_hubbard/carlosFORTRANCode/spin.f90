

module spin

  use constants

  implicit none
  private

  public :: wignerd, spinrot

! +----------------------------------------------------------------+
! |                                                                |
! | spin                                                           |
! |                                                                |
! |                                                                |
! | A collection of subroutines performing tasks related to spin   |
! | projection.                                                    |
! |                                                                |
! +----------------------------------------------------------------+

contains


subroutine wignerd (zval, lint, j, m, n, beta, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | wignerd  --  CAJH, 11.2012                                     |
! |                                                                |
! |                                                                |
! | Evaluate matrix elements of Wigner's small d-matrix            |
! |                                                                |
! |   d^j_{m,n}  =  <j;m | exp(-i*beta*Jy) | j;n>,                 |
! |                                                                |
! | for physical values (integers or half-integers) of the indices |
! | j, m, n. Here, Jy is an angular momentum operator and |j;m>    |
! | is an eigenstate of J^2 and Jz with eigenvalues j*(j+1) and m, |
! | respectively.                                                  |
! |                                                                |
! |   WARNING:  Be aware of the sign convention used!              |
! |             See discussion in Brink and Satchler.              |
! |                                                                |
! | NOTE: The input parameters j,m,n are all integers. The logical |
! |   variable lint should be used to control whether these        |
! |   coefficients should be interpreted as actual integers or     |
! |   as half-integers:                                            |
! |                                                                |
! |   lint = .true.,   integer coefficients                        |
! |        = .false.,  half-integer coefficients                   |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  lint = .false. but even j, m, or n              |
! |         =  2,  unphysical value: j < 0                         |
! |         =  3,  unphysical value: m > j or n > j                |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The evaluation of matrix elements is done using explicit       |
! | expressions for j <= 2 (obtained from Wikipedia). For j > 2,   |
! | the following general expression is used                       |
! |                                                                |
! |   d^j_{m,n} (beta) = sqrt [ (j+n)!(j-n)! / {(j+m)!(j-m)!} ] *  |
! |                      sum_s (-1)^s { binom (j+m,s) *            |
! |                                     binom (j-m,j-n-s) *        |
! |                                     cos(beta/2)^(2j+m-n-2s) *  |
! |                                     sin(beta/2)^(n-m+2s) },    |
! |                                                                |
! | where binom (j,k) is the corresponding binomial coefficient.   |
! | The sum over s is restricted to the argument of every          |
! | factorial being non-negative.                                  |
! |                                                                |
! | The above expression can be found in, e.g.,                    |
! |                                                                |
! |   D.M. Brink and G.R. Satchler, Angular Momentum [ p. 22 ].    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   lint - whether j, m, n should be interpreted as integers
  !            (.true.) or half-integers (.false.)
  !   j    - index j of Wigner matrix
  !   m    - index m of Wigner matrix
  !   n    - index n of Wigner matrix
  !   beta - angle to use in evaluation

  integer, intent(in) :: j, m, n
  logical, intent(in) :: lint

  real(kind=dp), intent(in) :: beta


  ! output variables

  !   zval - evaluated function, i.e., d^j_{m,n} (beta)
  !   iflg - return code

  integer, intent(out) :: iflg

  real(kind=dp), intent(out) :: zval


  ! other variables

  integer :: jwrk, mwrk, nwrk, sgf
  integer :: fac1, fac2, exp1, exp2, j1, j2, mins, maxs, s
  logical :: tstm, tstn, ltest

  real(kind=dp) :: cosbt, sinbt, sums, z, fz1, fz2


  ! constants

  real(kind=dp) :: f12
  parameter ( f12 = d1 / d2 )


  iflg = 0
  zval = d0

  ! check that j, m, n, are all odd if lint = .false.

  if ( .not. lint ) then
    if ( mod (abs(j),2) == 0 .or. &
       & mod (abs(m),2) == 0 .or. &
       & mod (abs(n),2) == 0 ) then
      iflg = 1
      return
    end if
  end if

  ! verify that values of j, m, n are physical

  if ( j < 0 ) then
    iflg = 2
    return
  end if

  tstm = abs(m) <= j
  tstn = abs(n) <= j

  if ( .not. tstm .or. .not. tstn ) then
    iflg = 3
    return
  end if


  ! Not all matrix elements are coded. The following relations are
  ! used to relate certain elements:

  !   d^j_{m,n} (beta)  =  (-1)^(m-n) * d^j_{n,m} (beta),
  !                     =  (-1)^(m-n) * d^j_{-m,-n} (beta),
  !                     =  d^j_{-n,-m} (beta).

  ! We choose to work with m >= 0, m >= abs(n).

  !   sgf - used to keep track of sign

  sgf  = 1

  ! set defaults

  jwrk = j
  mwrk = m
  nwrk = n

  ltest = ( lint .and. mod (n-m,2) == 1 ) .or. &
        & ( .not. lint .and. mod ((n-m)/2,2) == 1 ) 


  if ( m >= 0 ) then

    if ( m >= abs(n) ) then

      ! nothing to do

    else if ( m < n ) then      ! n > 0, n > m

      mwrk = n    ! use  d^j_{m,n}  =  (-1)^(m-n) * d^j_{n,m}
      nwrk = m

      if ( ltest )  sgf = -sgf

    else if ( m > n ) then      ! n < 0, -n > m

      mwrk = -n   ! use  d^j_{m,n}  =  d^j_{-n,-m}
      nwrk = -m

    end if

  else

    if ( -m >= abs(n) ) then

      mwrk = -m   ! use  d^j_{m,n}  =  (-1)^(m-n) * d^j_{-m,-n}
      nwrk = -n

      if ( ltest )  sgf = -sgf

    else if ( m < n ) then      ! n > 0, n > -m

      mwrk = n    ! use  d^j_{m,n}  =  (-1)^(m-n) * d^j_{n,m}
      nwrk = m

      if ( ltest )  sgf = -sgf

    else if ( m > n ) then      ! n < 0, -n > -m

      mwrk = -n   ! use  d^j_{m,n}  =  d^j_{-n,-m}
      nwrk = -m

    end if

  end if


  ! +---------------------------+
  ! |  evaluate matrix element  |
  ! +---------------------------+

  ! Explicit expressions are used for j <= 2. The expressions used can
  ! be found in Wikipedia or in Brink and Satchler.

  ! From now on, we adopt the notation
  !   d (j,m,n,beta)  ==  d^j_{m,n} (beta)

  z = d1


  ! ==  j = 0  ==

  !   d (0, 0, 0)  =  1


  ! ==  j = 1/2  ==

  !   d (1/2, 1/2, +1/2)  = +cos(beta/2)
  !   d (1/2, 1/2, -1/2)  = -sin(beta/2)

  if ( .not. lint .and. jwrk == 1 ) then
    if ( nwrk == +1 ) then
      z = + cos (beta/d2)
    else if ( nwrk == -1 ) then
      z = - sin (beta/d2)
    end if
  end if


  ! ==  j = 1  ==

  !   d (1, 1, +1)  = +1/2 * (1 + cos(beta))
  !   d (1, 1,  0)  = -sin(beta) / sqrt(2)
  !   d (1, 1, -1)  = +1/2 * (1 - cos(beta))
  !   d (1, 0,  0)  = +cos(beta)

  if ( lint .and. jwrk == 1 ) then
    if ( mwrk == 1 ) then
      if ( nwrk == +1 ) then
        z = + f12 * (d1 + cos (beta))
      else if ( nwrk == 0 ) then
        z = - sin (beta) / sqrt (d2)
      else if ( nwrk == -1 ) then
        z = + f12 * (d1 - cos (beta))
      end if

    else if ( mwrk == 0 ) then
      z = + cos (beta)
    end if
  end if


  ! ==  j = 3/2  ==

  !   d (3/2, 3/2, +3/2)  = +1/2 * (1 + cos(beta)) * cos(beta/2)
  !   d (3/2, 3/2, +1/2)  = -sqrt(3)/2 * (1 + cos(beta)) * sin(beta/2)
  !   d (3/2, 3/2, -1/2)  = +sqrt(3)/2 * (1 - cos(beta)) * cos(beta/2)
  !   d (3/2, 3/2, -3/2)  = -1/2 * (1 - cos(beta)) * sin(beta/2)
  !   d (3/2, 1/2, +1/2)  = +1/2 * (3*cos(beta) - 1) * cos(beta/2)
  !   d (3/2, 1/2, -1/2)  = -1/2 * (3*cos(beta) + 1) * sin(beta/2)

  if ( .not. lint .and. jwrk == 3 ) then
    if ( mwrk == 3 ) then
      if ( nwrk == +3 ) then
        z = + f12 * cos (beta/d2) * (d1 + cos (beta))
      else if ( nwrk == +1 ) then
        z = - f12 * sqrt(d3) * sin (beta/d2) * (d1 + cos (beta))
      else if ( nwrk == -1 ) then
        z = + f12 * sqrt(d3) * cos (beta/d2) * (d1 - cos (beta))
      else if ( nwrk == -3 ) then
        z = - f12 * sin (beta/d2) * (d1 - cos (beta))
      end if

    else if ( mwrk == 1 ) then
      if ( nwrk == +1 ) then
        z = + f12 * cos (beta/d2) * (d3 * cos (beta) - d1)
      else if ( nwrk == -1 ) then
        z = - f12 * sin (beta/d2) * (d3 * cos (beta) + d1)
      end if
    end if
  end if


  ! ==  j = 2  ==

  !   d (2, 2, +2)  = +1/4 * (1 + cos(beta))^2
  !   d (2, 2, +1)  = -1/2 * (1 + cos(beta)) * sin(beta)
  !   d (2, 2,  0)  = +sqrt(6)/4 * sin^2(beta)
  !   d (2, 2, -1)  = -1/2 * (1 - cos(beta)) * sin(beta)
  !   d (2, 2, -2)  = +1/4 * (1 - cos(beta))^2
  !   d (2, 1, +1)  = +1/2 * (1 + cos(beta)) * (2*cos(beta) - 1)
  !   d (2, 1,  0)  = -sqrt(3/2) * sin(beta) * cos(beta)
  !   d (2, 1, -1)  = +1/2 * (1 - cos(beta)) * (2*cos(beta) + 1)
  !   d (2, 0,  0)  = +1/2 * (3*cos^2(beta) - 1)

  if ( lint .and. jwrk == 2 ) then
    if ( mwrk == 2 ) then
      if ( nwrk == +2 ) then
        z = + f12*f12 * (d1 + cos (beta))**2
      else if ( nwrk == +1 ) then
        z = - f12 * (d1 + cos (beta)) * sin (beta)
      else if ( nwrk == 0 ) then
        z = + f12*f12 * sqrt (d3*d2) * (sin (beta))**2
      else if ( nwrk == -1 ) then
        z = - f12 * (d1 - cos (beta)) * sin (beta)
      else if ( nwrk == -2 ) then
        z = + f12*f12 * (d1 - cos (beta))**2
      end if

    else if ( mwrk == 1 ) then
      if ( nwrk == +1 ) then
        z = + f12 * (d1 + cos (beta)) * (d2 * cos (beta) - d1)
      else if ( nwrk == 0 ) then
        z = - sqrt(d3/d2) * sin (beta) * cos (beta)
      else if ( nwrk == -1 ) then
        z = + f12 * (d1 - cos (beta)) * (d2 * cos (beta) + d1)
      end if

    else if ( mwrk == 0 ) then
      z = + f12 * (d3 * (cos (beta))**2 - d1)
    end if
  end if


  ! ==  all other j:  j > 2  ==

  ! The matrix element is computed using the following expression:

  !   d (j,m,n,beta)  =  sqrt [ (j+n)! (j-n)! / { (j+m)! (j-m)! } ] *
  !                      sum_s (-1)^s { binom (j+m,s) *
  !                                     binom (j-m,j-n-s) *
  !                                     cos(beta/2)^(2j+m-n-2s) *
  !                                     sin(beta/2)^(n-m+2s) },

  ! where binom(j,k) is the corresponding binomial coefficient. The
  ! sum over s is restricted to the argument of every factorial being
  ! non-negative.

  if ( ( lint .and. jwrk >= 3 ) .or. &
     & ( .not. lint .and. jwrk >= 5 ) ) then

     ! compute pre-factor

     !   sqrt [ (j+n)! (j-n)! / { (j+m)! (j-m)! } ]

     ! We have chosen to work with m >= 0, which implies j+m >= j-m.
     ! Let us determine whether j+n >= j-n  or  j-n >= j+n.

     if ( nwrk >= 0 ) then
       j1 = jwrk + nwrk       ! j1 >= j2 in all cases
       j2 = jwrk - nwrk
     else
       j1 = jwrk - nwrk
       j2 = jwrk + nwrk
     end if

     ! Compute the ratios of factorials:

     !   fz1 = j1! / (j+m)!
     !   fz2 = j2! / (j-m)!

     if ( lint ) then
       fz1 = factorial_ratio (j1, jwrk+mwrk)
       fz2 = factorial_ratio (j2, jwrk-mwrk)
     else
       fz1 = factorial_ratio (j1/2, (jwrk+mwrk)/2)
       fz2 = factorial_ratio (j2/2, (jwrk-mwrk)/2)
     end if

     ! Store pre-factor in z.

     z = sqrt (fz1 * fz2)


     ! We are ready to compute the sum over s.
     ! There are four conditions on s from the factorials:

     !   s <= j + m,    s >= 0,
     !   s <= j - n,    s >= m - n.

     ! We have chosen to work with m >= 0, m >= abs(n).
     ! Since m-n >= 0, s >= m-n already implies s >= 0.

     ! Since m >= abs(n), j+m >= j-n in all cases.
     ! Therefore, s <= j-n already implies s <= j+m.

     ! We conclude that the interval for s is defined by

     !   m-n <= s <= j-n.

     if ( lint ) then
       mins = mwrk-nwrk
       maxs = jwrk-nwrk
     else
       mins = (mwrk-nwrk)/2
       maxs = (jwrk-nwrk)/2
     end if


     ! prepare loop to perform sum over s

     cosbt = cos (beta/d2)
     sinbt = sin (beta/d2)

       if ( lint ) then
         exp1 = 2*jwrk + mwrk - nwrk
         exp2 = nwrk - mwrk
       else
         exp1 = jwrk + (mwrk - nwrk)/2
         exp2 = (nwrk - mwrk)/2
       end if


     ! loop over s

     sums = d0

     do s = mins, maxs

       ! binomial coefficients

       if ( lint ) then
         fac1 = binomial (jwrk+mwrk, s)
         fac2 = binomial (jwrk-mwrk, jwrk-nwrk-s)
       else
         fac1 = binomial ((jwrk+mwrk)/2, s)
         fac2 = binomial ((jwrk-mwrk)/2, (jwrk-nwrk)/2-s)
       end if

       ! add term to sums

       if ( mod(s,2) == 0 ) then
         sums = sums + real (fac1,dp) * real (fac2,dp) &
                   & * cosbt**(exp1-2*s) * sinbt**(exp2+2*s)
       else
         sums = sums - real (fac1,dp) * real (fac2,dp) &
                   & * cosbt**(exp1-2*s) * sinbt**(exp2+2*s)
       end if
     end do

     z = z * sums
  end if


  ! account for factor of (-1)^(m-n) if needed

  z = real (sgf, dp) * z

  ! return computed value

  zval = z


  return
end subroutine wignerd



function binomial (j, k)

! +----------------------------------------------------------------+
! |                                                                |
! | binomial  --  CAJH, 11.2012                                    |
! |                                                                |
! |                                                                |
! | Compute the binomial coefficient                               |
! |                                                                |
! |   binomial (j,k)  =  j! / [ k! (j-k)! ]                        |
! |                                                                |
! | for integer, non-negative j and k, such that j >= k.           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   j, k - integers defining the binomial to evaluate

  integer, intent(in) :: j, k


  ! output variables

  !   binomial - binomial coefficient

  integer :: binomial


  ! other variables

  integer :: mink, maxk, i, z


  ! error checking:

  !   - j, k should be non-negative
  !   - k <= j

  if ( j < 0 .or. k < 0 ) then
    write (6, *) 'error: Negative j or k in binomial.'
    stop
  end if

  if ( k > j ) then
    write (6, *) 'error: k > j in binomial.'
    stop
  end if

  ! a few trivial cases, where binomial = 1

  if ( ( j == k ) .or. ( j <= 1 ) .or. ( k == 0 ) ) then
    binomial = 1
    return
  end if

  ! determine which factor is larger: k or j-k

  maxk = max (k, j-k)
  mink = min (k, j-k)

  ! compute z = j! / maxk!

  z = 1

  do i = j, maxk+1, -1
    z = z * i
  end do

  ! now compute z = z / mink!

  do i = 1, mink
    z = z / i
  end do

  binomial = z


  return
end function binomial



function factorial_ratio (j, k)

! +----------------------------------------------------------------+
! |                                                                |
! | factorial_ratio  --  CAJH, 11.2012                             |
! |                                                                |
! |                                                                |
! | Compute the ratio of factorials                                |
! |                                                                |
! |   factorial_ratio (j,k)  =  j! / k!                            |
! |                                                                |
! | for integer, non-negative j and k.                             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   j, k - integers defining the factorials to evaluate

  integer, intent(in) :: j, k


  ! output variables

  !   factorial_ratio - ratio of factorials j! / k!

  real(kind=dp) :: factorial_ratio


  ! other variables

  integer :: ij, ik, i

  real(kind=dp) :: z


  ! error checking:

  !   - j, k should be non-negative

  if ( j < 0 .or. k < 0 ) then
    write (6, *) 'error: Negative j or k in factorial_ratio.'
    stop
  end if

  ij = j
  ik = k

  if ( ij == 0 )  ij = 1
  if ( ik == 0 )  ik = 1


  ! if j == k, the ratio is equal to 1

  if ( ( j == k ) ) then
    factorial_ratio = d1
    return
  end if

  ! handle j > k and k > j independently

  z = d1

  if ( ij > ik ) then
    do i = ij, ik+1, -1
      z = z * real (i,dp)
    end do

  else
    do i = ik, ij+1, -1
      z = z / real (i,dp)
    end do
  end if

  factorial_ratio = z


  return
end function factorial_ratio



subroutine spinrot (nbs, alpha, beta, gamma, xmat, ymat, itype)

! +----------------------------------------------------------------+
! |                                                                |
! | spinrot  --  CAJH, 11.2012                                     |
! |                                                                |
! |                                                                |
! | Perform the matrix multiplication of an arbitrary matrix xmat  |
! | with the matrix representation of the spin rotation operator.  |
! |                                                                |
! | The following operations are supported:                        |
! |                                                                |
! |   itype  =  1,  ymat  =  R(Omega) . xmat                       |
! |          =  2,  ymat  =  R(Omega)* . xmat                      |
! |                                                                |
! | Both xmat and ymat are assumed to be full spin-orbital         |
! | matrices (dimension 2*nbs x 2*nbs).                            |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The spin rotation operator R(Omega) is given by                |
! |                                                                |
! |   R(Omega)  =  exp[-i*alpha*Sz] . exp[-i*beta*Sy] .            |
! |                exp[-i*gamma*Sz]                                |
! |                                                                |
! | The matrix representation of such an operator in a basis of    |
! | orthonormal single-fermion (s=1/2) states is then given by     |
! |                                                                |
! |   R(Omega)  =  R_alpha . R_beta . R_gamma,                     |
! |                                                                |
! |   R_alpha  =  [  exp(-i*alpha/2)          0       ]            |
! |               [         0         exp(+i*alpha/2) ]            |
! |                                                                |
! |   R_beta   =  [    cos(beta/2)     -sin(beta/2)   ]            |
! |               [    sin(beta/2)      cos(beta/2)   ]            |
! |                                                                |
! |   R_gamma  =  [  exp(-i*gamma/2)          0       ]            |
! |               [         0         exp(+i*gamma/2) ]            |
! |                                                                |
! | This subroutine exploits the very simple structure of the      |
! | rotation matrix to compute the products efficiently.           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nbs   - dimension of spin-blocks
  !   alpha - Euler angle alpha
  !   beta  - Euler angle beta
  !   gamma - Euler angle gamma
  !   xmat  - double complex matrix xmat
  !   itype - type of operation to perform (see above)
  !   ymat  - resulting double complex matrix ymat [ out ]

  integer, intent(in) :: nbs, itype

  real(kind=dp), intent(in) :: alpha, beta, gamma

  complex(kind=dp), dimension(2*nbs,2*nbs), intent(in) :: xmat
  complex(kind=dp), dimension(2*nbs,2*nbs), intent(inout) :: ymat


  ! other variables

  integer :: idum

  real(kind=dp), dimension(2,2) :: wigmat

  complex(kind=dp) :: f1a, f2a, f1g, f2g


  ! build Wigner's small d-matrix (1/2,1/2)

  call wignerd (wigmat(1,1), .false., 1, +1, +1, beta, idum)
  call wignerd (wigmat(1,2), .false., 1, +1, -1, beta, idum)
  call wignerd (wigmat(2,1), .false., 1, -1, +1, beta, idum)
  call wignerd (wigmat(2,2), .false., 1, -1, -1, beta, idum)

  ! calculate alpha and gamma exponentials

  f1a = exp (-zi*alpha/d2)
  f2a = exp (+zi*alpha/d2)

  f1g = exp (-zi*gamma/d2)
  f2g = exp (+zi*gamma/d2)


  select case (itype)

    case (1)

      ! ymat  =  R(Omega) . xmat

      ! uu block

      ymat(1:nbs,1:nbs) = &
           &   f1a * f1g * wigmat(1,1) * xmat(1:nbs,1:nbs) &
           & + f1a * f2g * wigmat(1,2) * xmat(nbs+1:2*nbs,1:nbs)

      ! ud block

      ymat(1:nbs,nbs+1:2*nbs) = &
           &   f1a * f1g * wigmat(1,1) * xmat(1:nbs,nbs+1:2*nbs) &
           & + f1a * f2g * wigmat(1,2) * xmat(nbs+1:2*nbs,nbs+1:2*nbs)

      ! du block

      ymat(nbs+1:2*nbs,1:nbs) = &
           &   f2a * f1g * wigmat(2,1) * xmat(1:nbs,1:nbs) &
           & + f2a * f2g * wigmat(2,2) * xmat(nbs+1:2*nbs,1:nbs)

      ! ud block

      ymat(nbs+1:2*nbs,nbs+1:2*nbs) = &
           &   f2a * f1g * wigmat(2,1) * xmat(1:nbs,nbs+1:2*nbs) &
           & + f2a * f2g * wigmat(2,2) * xmat(nbs+1:2*nbs,nbs+1:2*nbs)

    case (2)

      ! ymat  =  R(Omega)* . xmat

      ! uu block

      ymat(1:nbs,1:nbs) = &
           &   f2a * f2g * wigmat(1,1) * xmat(1:nbs,1:nbs) &
           & + f2a * f1g * wigmat(1,2) * xmat(nbs+1:2*nbs,1:nbs)

      ! ud block

      ymat(1:nbs,nbs+1:2*nbs) = &
           &   f2a * f2g * wigmat(1,1) * xmat(1:nbs,nbs+1:2*nbs) &
           & + f2a * f1g * wigmat(1,2) * xmat(nbs+1:2*nbs,nbs+1:2*nbs)

      ! du block

      ymat(nbs+1:2*nbs,1:nbs) = &
           &   f1a * f2g * wigmat(2,1) * xmat(1:nbs,1:nbs) &
           & + f1a * f1g * wigmat(2,2) * xmat(nbs+1:2*nbs,1:nbs)

      ! ud block

      ymat(nbs+1:2*nbs,nbs+1:2*nbs) = &
           &   f1a * f2g * wigmat(2,1) * xmat(1:nbs,nbs+1:2*nbs) &
           & + f1a * f1g * wigmat(2,2) * xmat(nbs+1:2*nbs,nbs+1:2*nbs)

  end select


  return
end subroutine spinrot


end module spin


