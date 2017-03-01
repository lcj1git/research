

module trdnmt

  use constants
  use util, only : error_alloc
  use linalg

  implicit none
  private
  save

  public :: setup_trdnmt, shutdown_trdnmt
  public :: trdnmt_ghf, trdnmt_uhf, trdnmt_rhf

! +----------------------------------------------------------------+
! |                                                                |
! | trdnmt                                                         |
! |                                                                |
! |                                                                |
! | This module constructs the transition density matrix between   |
! | two HF states |phi_1> and |phi_2>, given by                    |
! |                                                                |
! |                    <phi_1 | c!_i c_k |phi_2>                   |
! |   rho^[12]_ki  =  ---------------------------                  |
! |                         <phi_1 | phi_2>                        |
! |                                                                |
! |                =  sum_h   D1_{ih} . D1*_{kh}                   |
! |                +  sum_ph  D1_{ih} . Zc_{ph} . D1*_{kp},        |
! |                                                                |
! |                =  sum_hh' D1_{ih} . D2*_{kh'} . (Lc^-1)_{h'h}  |
! |                                                                |
! | where                                                          |
! |                                                                |
! |   Zc_{ph}  =  sum_h'  (D1^T . D2*)_{ph'} . (Lc^-1)_{h'h},      |
! |   Lc_{h'h} =  (D1^T . D2*)_{h'h}.                              |
! |                                                                |
! | This module assumes that an orthonormal basis is used.         |
! |                                                                |
! |                                                                |
! | NOTE:                                                          |
! |   In order to use this module,                                 |
! |     * setup_trdnmt *                                           |
! |   should be invoked before any calls to trdnmt_***.            |
! |   When the module is no longer needed, one may call            |
! |     * shutdown_trdnmt *                                        |
! |   to deallocate all scratch arrays.                            |
! |                                                                |
! +----------------------------------------------------------------+


  ! scratch arrays

  integer :: lwrkl
  integer, dimension(:), allocatable :: ipiv

  complex(kind=dp), dimension(:), allocatable :: wrkl
  complex(kind=dp), dimension(:,:), allocatable :: lcmt, scrl
  complex(kind=dp), dimension(:,:), allocatable :: lcmt_up, lcmt_dn

  !$omp  threadprivate(lwrkl, ipiv, wrkl, &
  !$omp&               lcmt, scrl, lcmt_up, lcmt_dn)


contains


subroutine setup_trdnmt (iwfnty, nbs, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_trdnmt  --  CAJH, 11.2012                                |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in trdnmt.         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   nbs    - dimension of spin-blocks
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons

  integer, intent(in) :: iwfnty, nbs, nup, ndn


  ! other variables

  integer :: istatus
  integer :: ist1, ist2, nb, ntot


  ! external functions

  integer, external :: ilaenv


  ist1 = 0
  ist2 = 0
  istatus = 0

  ntot = nup + ndn

  ! determine lwrkl

  if ( iwfnty == 1 .or. iwfnty == 2 ) then
    ! nb = 100
    nb = ilaenv (1, 'zgetri', ' ', nbs, -1, -1, -1)
    lwrkl = nbs*nb
  else if ( iwfnty == 3 ) then
    ! nb = 100
    nb = ilaenv (1, 'zgetri', ' ', 2*nbs, -1, -1, -1)
    lwrkl = 2*nbs*nb
  end if

  ! allocate space for lcmt

  if ( iwfnty == 1 ) then
    allocate (lcmt(nup,nup), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (lcmt(ntot,ntot), stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    allocate (lcmt_up(nup,nup), stat=ist1)
    end if

    if ( ndn > 0 ) then
    allocate (lcmt_dn(ndn,ndn), stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'lcmt', 'setup_trdnmt')

  ! allocate space for ipiv

  if ( iwfnty == 1 ) then
    allocate (ipiv(nup), stat=istatus)
  else if ( iwfnty == 2 .or. iwfnty == 3 ) then
    allocate (ipiv(ntot), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'ipiv', 'setup_trdnmt')

  ! allocate space for wrkl

  allocate (wrkl(lwrkl), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wrkl', 'setup_trdnmt')

  ! allocate space for scrl

  if ( iwfnty == 1 ) then
    allocate (scrl(nbs,nup), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (scrl(nbs,ntot), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scrl(2*nbs,ntot), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scrl', 'setup_trdnmt')

  return
end subroutine setup_trdnmt



subroutine shutdown_trdnmt (iwfnty, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_trdnmt  --  CAJH, 11.2012                             |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in trdnmt.       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons

  integer, intent(in) :: iwfnty, nup, ndn


  ! other variables

  integer :: istatus
  integer :: ist1, ist2

  ist1 = 0
  ist2 = 0
  istatus = 0

  ! deallocate space for lcmt

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (lcmt, stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    deallocate (lcmt_up, stat=ist1)
    end if

    if ( ndn > 0 ) then
    deallocate (lcmt_dn, stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (2, 'lcmt', 'shutdown_trdnmt')

  ! deallocate space for ipiv

  deallocate (ipiv, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ipiv', 'shutdown_trdnmt')

  ! deallocate space for wrkl

  deallocate (wrkl, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'wrkl', 'shutdown_trdnmt')

  ! deallocate space for scrl

  deallocate (scrl, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scrl', 'shutdown_trdnmt')

  return
end subroutine shutdown_trdnmt



subroutine trdnmt_ghf (nbs, nup, ndn, d1st, d2st, rho, ovpX, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | trdnmt_ghf  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | Build the transition density matrix rho between two GHF-type   |
! | determinants. That is, we build                                |
! |                                                                |
! |                     <phi_1 | c!_i c_k | phi_2>                 |
! |   rho_{ki}^[12]  =  -------------------------- ,               |
! |                          <phi_1 | phi_2>                       |
! |                                                                |
! | where |phi1>, |phi2> are HF-type determinants characterized by |
! | the orbital coefficients D1, D2. Additionally, the overlap     |
! | <phi1 | phi2> is evaluated.                                    |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  states |phi1>, |phi2> are orthogonal;           |
! |                rho was not constructed                         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   nbs  - dimension of spin-blocks
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons
  !   d1st - matrix D1* of orbital coefficients for |phi_1>
  !   d2st - matrix D2* of orbital coefficients for |phi_2>
  !   rho  - transition density matrix
  !   ovpX - overlap <D1 | D2>
  !   iflg - return code

  integer, intent(in) :: nbs, nup, ndn
  integer, intent(out) :: iflg

  complex(kind=dp), intent(out) :: ovpX
  complex(kind=dp), dimension(4*nbs*nbs), intent(in) :: d1st, d2st
  complex(kind=dp), dimension(2*nbs,2*nbs), intent(inout) :: rho


  ! other variables

  integer :: info, nbssq, ntot


  iflg = 0
  nbssq = nbs*nbs
  ovpX = z1

  ntot = nup + ndn

  ! compute Lc:  Lc = D1(:,1:N)^T . D2(:,1:N)*

  call zgemm ('c', 'n', ntot, ntot, 2*nbs, z1, d1st, 2*nbs, &
       & d2st, 2*nbs, z0, lcmt, ntot)

  ! LU factorization of Lc
  ! compute determinant of Lc
  ! invert Lc

  call zgetrf (ntot, ntot, lcmt, ntot, ipiv, info)

    ovpX = ovpX * det_zmat (ntot, lcmt, ipiv)

  call zgetri (ntot, lcmt, ntot, ipiv, wrkl, lwrkl, info)

    if ( info > 0 ) then
      iflg = 1
      return
    end if

  ! compute rho:  rho = D2(:,1:N)* . inv(Lc) . D1(:,1:N)^T

  call zgemm ('n', 'n', 2*nbs, ntot, ntot, z1, d2st, 2*nbs, &
       & lcmt, ntot, z0, scrl, 2*nbs)

  call zgemm ('n', 'c', 2*nbs, 2*nbs, ntot, z1, scrl, 2*nbs, &
       & d1st, 2*nbs, z0, rho, 2*nbs)


  return
end subroutine trdnmt_ghf



subroutine trdnmt_uhf (nbs, nup, ndn, d1st, d2st, rho, ovpX, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | trdnmt_uhf  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | ( UHF version of trdnmt_ghf. )                                 |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nbs  - dimension of spin-blocks
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons
  !   d1st - matrix D1* of orbital coefficients for |phi_1>
  !   d2st - matrix D2* of orbital coefficients for |phi_2>
  !   rho  - transition density matrix
  !   ovpX - overlap <D1 | D2>
  !   iflg - return code

  integer, intent(in) :: nbs, nup, ndn
  integer, intent(out) :: iflg

  complex(kind=dp), intent(out) :: ovpX
  complex(kind=dp), dimension(2*nbs*nbs), intent(in) :: d1st, d2st
  complex(kind=dp), dimension(nbs,2*nbs), intent(inout) :: rho


  ! other variables

  integer :: info_up, info_dn, nbssq, ntot


  iflg = 0
  nbssq = nbs*nbs
  ovpX = z1

  ntot = nup + ndn

  info_up = 0
  info_dn = 0

  ! compute Lc:  Lc = D1(:,1:N)^T . D2(:,1:N)*

  if ( nup > 0 ) then
  call zgemm ('c', 'n', nup, nup, nbs, z1, d1st, nbs, &
       & d2st, nbs, z0, lcmt_up, nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('c', 'n', ndn, ndn, nbs, z1, d1st(nbssq+1), nbs, &
       & d2st(nbssq+1), nbs, z0, lcmt_dn, ndn)
  end if

  ! LU factorization of Lc
  ! compute determinant of Lc
  ! invert Lc

  if ( nup > 0 ) then
  call zgetrf (nup, nup, lcmt_up, nup, ipiv, info_up)
  end if

  if ( ndn > 0 ) then
  call zgetrf (ndn, ndn, lcmt_dn, ndn, ipiv(nup+1), info_dn)
  end if

    if ( nup > 0 ) then
    ovpX = ovpX * det_zmat (nup, lcmt_up, ipiv)
    end if

    if ( ndn > 0 ) then
    ovpX = ovpX * det_zmat (ndn, lcmt_dn, ipiv(nup+1))
    end if

  if ( nup > 0 ) then
  call zgetri (nup, lcmt_up, nup, ipiv, wrkl, lwrkl, info_up)
  end if

  if ( ndn > 0 ) then
  call zgetri (ndn, lcmt_dn, ndn, ipiv(nup+1), wrkl, lwrkl, info_dn)
  end if

    if ( info_up > 0 .or. info_dn > 0 ) then
      iflg = 1
      return
    end if

  ! compute rho:  rho = D2(:,1:N)* . inv(Lc) . D1(:,1:N)^T

  if ( nup > 0 ) then
  call zgemm ('n', 'n', nbs, nup, nup, z1, d2st, nbs, &
       & lcmt_up, nup, z0, scrl(1,1), nbs)

  call zgemm ('n', 'c', nbs, nbs, nup, z1, scrl(1,1), nbs, &
       & d1st, nbs, z0, rho(1,1), nbs)
  else
    rho(1:nbs,1:nbs) = z0
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', nbs, ndn, ndn, z1, d2st(nbssq+1), nbs, &
       & lcmt_dn, ndn, z0, scrl(1,nup+1), nbs)

  call zgemm ('n', 'c', nbs, nbs, ndn, z1, scrl(1,nup+1), nbs, &
       & d1st(nbssq+1), nbs, z0, rho(1,nbs+1), nbs)
  else
    rho(1:nbs,nbs+1:2*nbs) = z0
  end if


  return
end subroutine trdnmt_uhf



subroutine trdnmt_rhf (nbs, nup, d1st, d2st, rho, ovpX, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | trdnmt_rhf  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | ( RHF version of trdnmt_ghf. )                                 |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nbs  - dimension of spin-blocks
  !   nup  - number of spin-up electrons
  !   d1st - matrix D1* of orbital coefficients for |phi_1>
  !   d2st - matrix D2* of orbital coefficients for |phi_2>
  !   rho  - transition density matrix
  !   ovpX - overlap <D1 | D2>
  !   iflg - return code

  integer, intent(in) :: nbs, nup
  integer, intent(out) :: iflg

  complex(kind=dp), intent(out) :: ovpX
  complex(kind=dp), dimension(nbs*nbs), intent(in) :: d1st, d2st
  complex(kind=dp), dimension(nbs,nbs), intent(inout) :: rho


  ! other variables

  integer :: info, nbssq

  complex(kind=dp) :: ztr


  iflg = 0
  nbssq = nbs*nbs
  ovpX = z1

  ! compute Lc:  Lc = D1(:,1:N)^T . D2(:,1:N)*

  call zgemm ('c', 'n', nup, nup, nbs, z1, d1st, nbs, &
       & d2st, nbs, z0, lcmt, nup)

  ! LU factorization of Lc
  ! compute determinant of Lc
  ! invert Lc

  call zgetrf (nup, nup, lcmt, nup, ipiv, info)

    ztr  = det_zmat (nup, lcmt, ipiv)
    ovpX = ovpX * ztr * ztr

  call zgetri (nup, lcmt, nup, ipiv, wrkl, lwrkl, info)

    if ( info > 0 ) then
      iflg = 1
      return
    end if

  ! compute rho:  rho = D2(:,1:N)* . inv(Lc) . D1(:,1:N)^T

  call zgemm ('n', 'n', nbs, nup, nup, z1, d2st, nbs, &
       & lcmt, nup, z0, scrl, nbs)

  call zgemm ('n', 'c', nbs, nbs, nup, z1, scrl, nbs, &
       & d1st, nbs, z0, rho, nbs)


  return
end subroutine trdnmt_rhf


end module trdnmt


