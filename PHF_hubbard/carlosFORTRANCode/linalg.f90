

module linalg

  use constants

  implicit none
  public

! +----------------------------------------------------------------+
! |                                                                |
! | linalg                                                         |
! |                                                                |
! |                                                                |
! | A collection of subroutines that perform linear algebra        |
! | operations not included in BLAS or LaPack.                     |
! |                                                                |
! +----------------------------------------------------------------+

contains


subroutine zmat_square (ioper, uplo, idim, amat)

! +----------------------------------------------------------------+
! |                                                                |
! | zmat_square  --  CAJH, 12.2012                                 |
! |                                                                |
! |                                                                |
! | Expand either the upper or lower triangle of a double complex  |
! | square matrix [dimension idim x idim] amat into the full       |
! | matrix. ioper determines how amat should be expanded:          |
! |                                                                |
! |   ioper  =  1,  amat Hermitian                                 |
! |          =  2,  amat symmetric                                 |
! |          =  3,  amat anti-Hermitian                            |
! |          =  4,  amat anti-symmetric                            |
! |                                                                |
! | The character input variable uplo controls the behavior of     |
! | this subroutine:                                               |
! |                                                                |
! |   uplo  =  'u' or 'U'                                          |
! |                                                                |
! |     The upper triangle of amat is referenced.                  |
! |     The lower triangle is updated.                             |
! |                                                                |
! |   uplo  =  'l' or 'L'                                          |
! |                                                                |
! |     The lower triangle of amat is referenced.                  |
! |     The upper triangle is updated.                             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   ioper - type of operation to perform
  !   uplo  - determines which triangle should be updated
  !   idim  - dimension of the matrix amat
  !   amat  - matrix to be squared [ updated ]

  integer, intent(in) :: ioper, idim

  character(len=*), intent(in) :: uplo

  complex(kind=dp), dimension(idim,idim), intent(inout) :: amat


  ! other variables

  integer :: j, k


  if ( uplo(1:1) == 'u' .or. uplo(1:1) == 'U' ) then

    select case (ioper)

      case (1)
        do j = 1, idim
          do k = 1, j-1
            amat(j,k) = conjg(amat(k,j))
          end do
        end do

      case (2)
        do j = 1, idim
          do k = 1, j-1
            amat(j,k) = amat(k,j)
          end do
        end do

      case (3)
        do j = 1, idim
          do k = 1, j-1
            amat(j,k) = -conjg(amat(k,j))
          end do
        end do

      case (4)
        do j = 1, idim
          do k = 1, j-1
            amat(j,k) = -amat(k,j)
          end do
        end do

      case default
        write (6, *) 'error: Unrecognized ioper in zmat_square.'
        stop
    end select

  else if ( uplo(1:1) == 'l' .or. uplo(1:1) == 'L' ) then

    select case (ioper)

      case (1)
        do j = 1, idim
          do k = 1, j-1
            amat(k,j) = conjg(amat(j,k))
          end do
        end do

      case (2)
        do j = 1, idim
          do k = 1, j-1
            amat(k,j) = amat(j,k)
          end do
        end do

      case (3)
        do j = 1, idim
          do k = 1, j-1
            amat(k,j) = -conjg(amat(j,k))
          end do
        end do

      case (4)
        do j = 1, idim
          do k = 1, j-1
            amat(k,j) = -amat(j,k)
          end do
        end do

      case default
        write (6, *) 'error: Unrecognized ioper in zmat_square.'
        stop
    end select

  else
    write (6, *) 'error: Unrecognized uplo in zmat_square.'
    stop
  end if


  return
end subroutine zmat_square



subroutine zmat_symm (ioper, uplo, idim, amat)

! +----------------------------------------------------------------+
! |                                                                |
! | zmat_symm  --  CAJH, 12.2012                                   |
! |                                                                |
! |                                                                |
! | Perform either of the following operations to the double       |
! | complex square matrix [dimension idim x idim] amat:            |
! |                                                                |
! |   ioper  =  1,  amat = 1/2 * (amat + amat!)   [Hermitian]      |
! |          =  2,  amat = 1/2 * (amat + amat^T)  [symmetric]      |
! |          =  3,  amat = 1/2 * (amat - amat!)   [anti-Hermitian] |
! |          =  4,  amat = 1/2 * (amat - amat^T)  [anti-symmetric] |
! |                                                                |
! | The character input variable uplo controls the behavior of     |
! | this subroutine:                                               |
! |                                                                |
! |   uplo  =  'u' or 'U'                                          |
! |                                                                |
! |     The upper triangle of amat is returned with the            |
! |     symmetrized version. The lower triangle is untouched.      |
! |                                                                |
! |   uplo  =  'l' or 'L'                                          |
! |                                                                |
! |     The lower triangle of amat is returned with the            |
! |     symmetrized version. The upper triangle is untouched.      |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   ioper - type of operation to perform
  !   uplo  - determines which triangle should be symmetrized
  !   idim  - dimension of the matrix amat
  !   amat  - matrix to symmetrize [ updated ]

  integer, intent(in) :: ioper, idim

  character(len=*), intent(in) :: uplo

  complex(kind=dp), dimension(idim,idim), intent(inout) :: amat


  ! other variables

  integer :: j, k


  ! constants

  real(kind=dp) :: f12
  parameter ( f12 = 0.50_dp )


  if ( uplo(1:1) == 'u' .or. uplo(1:1) == 'U' ) then

    select case (ioper)

      case (1)
        do j = 1, idim
          do k = 1, j
            amat(k,j) = f12 * (amat(k,j) + conjg(amat(j,k)))
          end do
        end do

      case (2)
        do j = 1, idim
          do k = 1, j
            amat(k,j) = f12 * (amat(k,j) + amat(j,k))
          end do
        end do

      case (3)
        do j = 1, idim
          do k = 1, j
            amat(k,j) = f12 * (amat(k,j) - conjg(amat(j,k)))
          end do
        end do

      case (4)
        do j = 1, idim
          do k = 1, j
            amat(k,j) = f12 * (amat(k,j) - amat(j,k))
          end do
        end do

      case default
        write (6, *) 'error: Unrecognized ioper in zmat_symm.'
        stop
    end select

  else if ( uplo(1:1) == 'l' .or. uplo(1:1) == 'L' ) then

    select case (ioper)

      case (1)
        do j = 1, idim
          do k = 1, j
            amat(j,k) = f12 * (amat(j,k) + conjg(amat(k,j)))
          end do
        end do

      case (2)
        do j = 1, idim
          do k = 1, j
            amat(j,k) = f12 * (amat(j,k) + amat(k,j))
          end do
        end do

      case (3)
        do j = 1, idim
          do k = 1, j
            amat(j,k) = f12 * (amat(j,k) - conjg(amat(k,j)))
          end do
        end do

      case (4)
        do j = 1, idim
          do k = 1, j
            amat(j,k) = f12 * (amat(j,k) - amat(k,j))
          end do
        end do

      case default
        write (6, *) 'error: Unrecognized ioper in zmat_symm.'
        stop
    end select

  else
    write (6, *) 'error: Unrecognized uplo in zmat_symm.'
    stop
  end if


  return
end subroutine zmat_symm



function trab_zmat (n, amat, lda, bmat, ldb, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | trab_zmat  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Given two double complex square (dimension n x n) matrices     |
! | A and B, compute the trace of the product:                     |
! |                                                                |
! |   trab_zmat = trace ( A . B )                                  |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  lda < n or ldb < n                              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   n    - dimension of matrices A, B
  !   amat - double complex square matrix
  !   bmat - double complex square matrix
  !   lda  - leading dimension of A
  !   ldb  - leading dimension of B
  !   iflg - return code
  !   trab_zmat - trace of A . B

  integer, intent(in) :: n, lda, ldb
  integer, intent(out) :: iflg

  complex(kind=dp) :: trab_zmat
  complex(kind=dp), dimension(lda,n), intent(in) :: amat
  complex(kind=dp), dimension(ldb,n), intent(in) :: bmat


  ! other variables

  integer :: j, k
  complex(kind=dp) :: ztr


  iflg = 0
  ztr = z0

  if ( lda < n .or. ldb < n ) then
    iflg = 1
    trab_zmat = ztr
    return
  end if

  do j = 1, n
    do k = 1, n
      ztr = ztr + amat(j,k) * bmat(k,j)
    end do
  end do

  trab_zmat = ztr


  return
end function trab_zmat



function det_zmat (n, xlu, ipiv)

! +----------------------------------------------------------------+
! |                                                                |
! | det_zmat  --  CAJH, 11.2012                                    |
! |                                                                |
! |                                                                |
! | Compute the determinant of a (square) double complex matrix X  |
! | [ dimension n x n ] given its LU decomposition.                |
! |                                                                |
! | det_zmat requires a previous call to LAPACK's zgetrf to obtain |
! | the LU decomposition of X, along with the vector of pivoting   |
! | indices ipiv. The LU decomposition is returned by zgetrf in    |
! | the same place that X was originally stored in. The lower      |
! | triangle corresponds to L, while the upper triangle (including |
! | the diagonal elements) corresponds to U. Note that the         |
! | diagonal elements of L (all equal to unity) are not stored.    |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The LU decomposition of X generated by LAPACK is of the form   |
! |                                                                |
! |   X = P . L . U,                                               |
! |                                                                |
! | where P is a permutation matrix (stored as a vector of         |
! | pivoting indices), and L and U are, respectively, lower and    |
! | upper triangular matrices.                                     |
! |                                                                |
! | The determinant of X can then be computed as                   |
! |                                                                |
! |   det(X) = det(P) * det(L) * det(U).                           |
! |                                                                |
! | The determinants of L and U are easily computed due to their   |
! | triangular nature. The determinant of a triangular matrix T    |
! | is given by the product of the diagonal elements, i.e.,        |
! |                                                                |
! |   det(T) = prod_i T(i,i).                                      |
! |                                                                |
! | Because the diagonal elements of L are all unity, we have that |
! | det(L) = 1.                                                    |
! |                                                                |
! | Given that P is a permutation matrix, we have det(P) = +/- 1,  |
! | depending on the signature of the permutation.                 |
! |                                                                |
! | det_zmat then computes the determinant of X as follows         |
! |                                                                |
! |   det(X)  =  +/- det(U)                                        |
! |           =  +/- prod_i  U(i,i)                                |
! |                                                                |
! | In order to avoid overflows or underflows as much as possible, |
! | we accumulate the logs of the absolute values.                 |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   n    - dimension of matrix X
  !   xlu  - LU decomposition of matrix X (use LAPACK)
  !   ipiv - vector of pivoting indices (use LAPACK)
  !   det_zmat - determinant of matrix X

  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: ipiv

  complex(kind=dp) :: det_zmat
  complex(kind=dp), dimension(n,n), intent(in) :: xlu


  ! other variables

  integer :: j, sgn

  real(kind=dp) :: xr, xi
  real(kind=dp) :: xphs, xabs


  ! determine the overall sign of the determinant by computing the
  ! signature of the permutation (from permutation matrix)

  sgn = 1

  do j = 1, n
    if ( ipiv(j) /= j )  sgn = -sgn
  end do


  xphs = d0  ! phase of determinant
  xabs = d0  ! absolute value of determinant

  ! loop over diagonal elements of the LU decomposition of X (i.e.,
  ! the diagonal elements in U)

  do j = 1, n

    xr = real  (xlu(j,j))  ! split into real and imaginary
    xi = aimag (xlu(j,j))

    ! accumulate phase and absolute value

    xphs = xphs + atan2 (xi,xr)
    xabs = xabs + log (sqrt (xr*xr + xi*xi))

  end do


  ! compute determinant

  det_zmat = exp (xabs) * exp (zi*xphs)

  if ( sgn == -1 )  det_zmat = -det_zmat


  return
end function det_zmat



subroutine zgenevs (uplo, idim, amat, bmat, ldim, xvec, aval, bval, &
     & ascr, bscr, sscr, work, lwork, rwork, thrsh, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | zgenevs  --  CAJH, 11.2012                                     |
! |                                                                |
! |                                                                |
! | Given a Hermitian matrix A and a Hermitian positive semi-      |
! | definite matrix B (square matrices, dimension idim x idim),    |
! | solve the generalized eigenvalue problem                       |
! |                                                                |
! |   A . X  =  B . X . a,                                         |
! |                                                                |
! | where a is the diagonal matrix of generalized eigenvalues.     |
! |                                                                |
! | This routine is non-destructive:                               |
! |   A (amat), B (bmat) are left intact.                          |
! |                                                                |
! | The parameter thrsh (input) determines the spread of           |
! | eigenvalues of B allowed for linear independency.              |
! |                                                                |
! | The following quantities are returned on output:               |
! |   ldim - dimension of linearly independent subspace            |
! |   xvec - generalized eigenvectors of A                         |
! |   aval - generalized eigenvalues of A                          |
! |   bval - eigenvalues of B (in descending order)                |
! |                                                                |
! | The character variable uplo determines whether the upper or    |
! | lower triangles of A, B should be referenced:                  |
! |                                                                |
! |   uplo = 'u' or 'U'                                            |
! |     The upper triangles of A, B are referenced.                |
! |     The lower triangles are not used.                          |
! |                                                                |
! |   uplo = 'l' or 'L'                                            |
! |     The lower triangles of A, B are referenced.                |
! |     The upper triangles are not used.                          |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  B is not positive semi-definite                 |
! |         =  2,  allocated dimension lwork is too small          |
! |         =  3,  diagonalization of B failed                     |
! |         =  4,  diagonalization of A failed                     |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | Several scratch arrays are required:                           |
! |   ascr - dimension idim x idim                                 |
! |   bscr - dimension idim x idim                                 |
! |   sscr - dimension idim x idim                                 |
! |                                                                |
! | Additionally, we require two scratch arrays for LAPACK's       |
! | diagonalization:                                               |
! |   work  - complex array of dimension lwork                     |
! |   rwork - double precision array of dimension 7*idim           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   uplo  - determines which triangle should be referenced
  !   idim  - dimension of matrices A, B
  !   amat  - matrix A
  !   bmat  - matrix B
  !   ascr  - idim*idim scratch matrix
  !   bscr  - idim*idim scratch matrix
  !   sscr  - idim*idim scratch matrix
  !   work  - lwork scratch space
  !   lwork - allocated space for scratch space work
  !   rwork - double precision scratch space (dim 7*idim)
  !   thrsh - maximum allowed spread of eigenvalues of B to
  !           denote linear independency

  character(len=*), intent(in) :: uplo

  integer, intent(in) :: idim, lwork

  real(kind=dp), intent(in) :: thrsh
  real(kind=dp), dimension(7*idim), intent(inout) :: rwork

  complex(kind=dp), dimension(lwork), intent(inout) :: work
  complex(kind=dp), dimension(idim,idim), intent(in) :: amat, bmat
  complex(kind=dp), dimension(idim*idim), intent(inout) :: ascr, bscr, sscr


  ! output variables

  !   ldim  - dimension of linearly independent subspace
  !   xvec  - ldim eigenvectors X of A (allocated dim idim x idim)
  !   aval  - sorted generalized eigenvalues of A
  !   bval  - sorted eigenvalues of B (descending order)
  !   iflg  - return code

  integer, intent(out) :: ldim, iflg

  real(kind=dp), dimension(idim), intent(out) :: aval, bval

  complex(kind=dp), dimension(idim*idim), intent(out) :: xvec


  ! other variables

  integer :: idum, inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  integer, dimension(5*idim) :: iwork
  integer, dimension(idim) :: ifail

  integer :: k, ind1, ind2, idimsq
  integer :: nb, lwrk, info

  real(kind=dp) :: s


  ! constants

  !   eps  - estimate of machine precision

  real(kind=dp) :: f12, eps

  parameter ( f12 = d1/d2, eps = 1.0e-16_dp )


  ! external functions

  integer, external :: ilaenv
  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  ldim = idim
  iflg = 0

  ! solve the trivial case of idim = 1

  if ( idim == 1 ) then
    s = d1 / sqrt (real(bmat(1,1),dp))

    xvec(1) = cmplx (s, d0, dp)
    aval(1) = s * real(amat(1,1),dp) * s
    bval(1) = real(bmat(1,1),dp)

    return
  end if


  ! check that allocated dimension of work array is sufficient

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', idim, -1, -1, -1)
  lwrk = max (1, (nb+1)*idim)

  if ( lwork < lwrk ) then
    iflg = 2
    return
  end if


  ! diagonalize B matrix
  ! set to negative to recover largest eigenvalues on top

  idimsq = idim*idim

  bscr(1:idimsq) = reshape (bmat, (/ idim*idim /))
  bscr(1:idimsq) = -bscr(1:idimsq)

  call zheevx ('v', 'a', uplo, idim, bscr, idim, ddum, ddum, &
       & idum, idum, abstol, inum, bval, sscr, idim, work, &
       & lwork, rwork, iwork, ifail, info)

  bscr(1:idimsq) = sscr(1:idimsq)

  ! call zheev ('v', uplo, idim, bscr, idim, bval, work, lwork, &
  !      & rwork, info)

  if ( info /= 0 ) then
    iflg = 3
    return
  end if

  bval(1:idim) = -bval(1:idim)

  ! is B positive semi-definite?

  do k = 1, idim
    if ( bval(k) < -eps ) then
      iflg = 1
      return
    end if
  end do


  ! determine dimension of linearly-independent subspace
  ! prepare b^(-1/2) . vec(B)

  do k = 2, idim
    if ( abs(bval(k)) / abs(bval(1)) < thrsh ) then
      ldim = k-1
      exit
    end if
  end do

  do k = 1, ldim
    ind1 = (k-1)*idim + 1
    ind2 = ind1 + idim-1

    bscr(ind1:ind2) = bscr(ind1:ind2) &
                  & * d1/sqrt(bval(k))
  end do


  ! project A into linearly independent subspace

  call zhemm ('l', uplo, idim, ldim, z1, amat, idim, &
       & bscr, idim, z0, sscr, idim)
  call zgemm ('c', 'n', ldim, ldim, idim, z1, bscr, idim, &
       & sscr, idim, z0, ascr, ldim)


  ! diagonalize A

  call zheevx ('v', 'a', 'u', ldim, ascr, ldim, ddum, ddum, &
       & idum, idum, abstol, inum, aval, sscr, ldim, work, &
       & lwork, rwork, iwork, ifail, info)

  ascr(1:ldim*ldim) = sscr(1:ldim*ldim)

  ! call zheev ('v', 'u', ldim, ascr, ldim, aval, work, lwork, &
  !      & rwork, info)

  aval(ldim+1:idim) = d0

  if ( info /= 0 ) then
    iflg = 4
    return
  end if


  ! transform the eigenvectors of A into original basis
  !   Avec  =  b^(-1/2) . Bvec . Avec

  call zgemm ('n', 'n', idim, ldim, ldim, z1, bscr, idim, &
       & ascr, ldim, z0, sscr, idim)

  xvec(1:idim*ldim) = sscr(1:idim*ldim)
  xvec(idim*ldim+1:idimsq) = z0


  return
end subroutine zgenevs


end module linalg


