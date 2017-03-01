

module diagh11

  use constants
  use i2sint
  use util
  use linalg
  use trdnmt
  use erictr

  implicit none
  private
  save

  public :: setup_diagh11, shutdown_diagh11
  public :: diagh11_rhf, diagh11_uhf, diagh11_ghf

! +----------------------------------------------------------------+
! |                                                                |
! | diagh11                                                        |
! |                                                                |
! |                                                                |
! | A collection of subroutines that diagonalize the H^11          |
! | Hamiltonian in order to uniquely define the occupied (and      |
! | virtual) orbitals in a determinant, as well as to associate an |
! | energy with each of them.                                      |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | Given a Slater determinant |Phi>, The hole-hole sector of H^11 |
! | is defined as                                                  |
! |                                                                |
! |                  <Phi | b_h!  H  b_h' |Phi>                    |
! |   H^11_{hh'}  =  --------------------------                    |
! |                         <Phi | Phi>                            |
! |                                                                |
! |               -  d_{hh'} * <Phi | H | Phi> / <Phi | Phi>       |
! |                                                                |
! | Similarly, the particle-particle sector is defined as          |
! |                                                                |
! |                  <Phi | b_p  H  b_p'! |Phi>                    |
! |   H^11_{pp'}  =  --------------------------                    |
! |                         <Phi | Phi>                            |
! |                                                                |
! |               -  d_{pp'} * <Phi | H | Phi> / <Phi | Phi>       |
! |                                                                |
! | Here, b_h! and b_p are, respectively, the creation operator    |
! | of a hole state and the annihilation operator of a particle    |
! | state defined in terms of |Phi>.                               |
! |                                                                |
! | The hole-hole sector of H^11 can be evaluated as               |
! |                                                                |
! |   H^11_{hh'}  =  [ D! . rho^T . D ]_{hh'} * E                  |
! |         -  [ D! . rho^T . (T + Gam)^T . rho^T . D ]_{hh'}      |
! |         -  d_{hh'}*E.                                          |
! |                                                                |
! | The particle-particle sector of H^11 can be evaluated as       |
! |                                                                |
! |   H^11_{pp'}  =  [ D^T . (I-rho) . D* ]_{pp'} * E              |
! |         +  [ D^T . (I-rho) . (T + Gam) . (I-rho) . D* ]_{pp'}  |
! |         -  d_{pp'}*E.                                          |
! |                                                                |
! | Here, E is the energy of the state |Phi>.                      |
! |                                                                |
! +----------------------------------------------------------------+


  ! global variables
  ! scratch arrays

  integer :: lwork

  integer, dimension(:), allocatable :: iwork, ifail

  real(kind=dp), dimension(:), allocatable :: rwork
  complex(kind=dp), dimension(:), allocatable :: work

  complex(kind=dp), dimension(:), allocatable :: dst, dnost, scrt
  complex(kind=dp), dimension(:,:), allocatable :: rho, gam, tpgam
  complex(kind=dp), dimension(:,:), allocatable :: scr_spb

  complex(kind=dp), dimension(:,:), allocatable :: hmath, hmatp
  complex(kind=dp), dimension(:,:), allocatable :: hmath_up, hmatp_up
  complex(kind=dp), dimension(:,:), allocatable :: hmath_dn, hmatp_dn

  complex(kind=dp), dimension(:,:), allocatable :: hmt1h, hmt2h
  complex(kind=dp), dimension(:,:), allocatable :: hmt1h_up, hmt2h_up
  complex(kind=dp), dimension(:,:), allocatable :: hmt1h_dn, hmt2h_dn

  complex(kind=dp), dimension(:,:), allocatable :: hmt1p, hmt2p
  complex(kind=dp), dimension(:,:), allocatable :: hmt1p_up, hmt2p_up
  complex(kind=dp), dimension(:,:), allocatable :: hmt1p_dn, hmt2p_dn

  real(kind=dp), dimension(:), allocatable :: valh, valp
  real(kind=dp), dimension(:), allocatable :: valh_up, valp_up
  real(kind=dp), dimension(:), allocatable :: valh_dn, valp_dn

  complex(kind=dp), dimension(:,:), allocatable :: vech, vecp
  complex(kind=dp), dimension(:,:), allocatable :: vech_up, vecp_up
  complex(kind=dp), dimension(:,:), allocatable :: vech_dn, vecp_dn

  ! scratch arrays for hhmat, ppmat

  complex(kind=dp), dimension(:,:), allocatable :: scr1, scr2


contains


subroutine setup_diagh11 (iwfnty, norb, nbas, nbct, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_diagh11  --  CAJH, 01.2013                               |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in diagh11_*.      |
! |                                                                |
! | This subroutine should be called before diagh11. After the     |
! | calling program is done with all calls to diagh11_*, the       |
! | shutdown routine can be used to clear all scratch arrays.      |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty  - type of wavefunction to use
  !   norb    - number of orbitals
  !   nbas    - number of basis functions
  !   nbct    - number of Cartesian basis functions
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons

  integer, intent(in) :: iwfnty, norb, nbas, nbct, nup, ndn


  ! other variables

  integer :: nosq, ntot, nb
  integer :: nscr, nscr_up, nscr_dn
  integer :: istatus


  ! functions

  integer, external :: ilaenv


  ! set some important variables

  nosq = norb*norb
  ntot = nup + ndn

  nscr = 0
  nscr_up = 0
  nscr_dn = 0
  istatus = 0

  if ( iwfnty == 1 ) then
    nscr = max (nup, norb-nup)
  else if ( iwfnty == 2 ) then
    nscr_up = max (nup, norb-nup)
    nscr_dn = max (ndn, norb-ndn)
    nscr    = nscr_up + nscr_dn
  else if ( iwfnty == 3 ) then
    nscr = max (ntot, 2*norb-ntot)
  end if


  ! allocate space for dst

  if ( iwfnty == 1 ) then
    allocate (dst(nosq), dnost(nosq), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (dst(2*nosq), dnost(2*nosq), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (dst(4*nosq), dnost(4*nosq), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'dst, dnost', 'setup_diagh11')

  ! allocate space for scrt

  if ( iwfnty == 1 ) then
    allocate (scrt(norb*nscr), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (scrt(norb*max(nscr_up,nscr_dn)), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scrt(2*norb*nscr), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scrt', 'setup_diagh11')

  ! allocate space for rho

  if ( iwfnty == 1 ) then
    allocate (rho(norb,norb), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (rho(norb,2*norb), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (rho(2*norb,2*norb), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'rho', 'setup_diagh11')

  ! allocate space for gam, tpgam

  if ( iwfnty == 1 ) then
    allocate (gam(norb,norb), tpgam(norb,norb), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (gam(norb,2*norb), tpgam(norb,2*norb), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (gam(2*norb,2*norb), tpgam(2*norb,2*norb), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'gam', 'setup_diagh11')

  ! scratch spin-block for GHF

  if ( iwfnty == 3 ) then
    allocate (scr_spb(norb,norb), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'scr_spb', 'setup_diagh11')
  end if


  ! setup trdnmt scratch

  call setup_trdnmt (iwfnty, norb, nup, ndn)

  ! setup erictr scratch

  call setup_erictr (iwfnty, norb, nbas, nbct)


  ! allocate space for hmath, hmatp

  if ( iwfnty == 1 ) then
    allocate (hmath(nup,nup), hmatp(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmat?', 'setup_diagh11')

  else if ( iwfnty == 2 ) then
    allocate (hmath_up(nup,nup), hmatp_up(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmat?_up', 'setup_diagh11')

    allocate (hmath_dn(ndn,ndn), hmatp_dn(norb-ndn,norb-ndn), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmat?_dn', 'setup_diagh11')

  else if ( iwfnty == 3 ) then
    allocate (hmath(ntot,ntot), hmatp(2*norb-ntot,2*norb-ntot), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmat?', 'setup_diagh11')
  end if


  ! allocate space for hmt1, hmt2

  if ( iwfnty == 1 ) then
    allocate (hmt1h(nup,nup), hmt2h(nup,nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?h', 'setup_diagh11')

    allocate (hmt1p(norb-nup,norb-nup), &
            & hmt2p(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?p', 'setup_diagh11')

  else if ( iwfnty == 2 ) then
    allocate (hmt1h_up(nup,nup), hmt2h_up(nup,nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?h_up', 'setup_diagh11')

    allocate (hmt1p_up(norb-nup,norb-nup), &
            & hmt2p_up(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?p_up', 'setup_diagh11')

    allocate (hmt1h_dn(ndn,ndn), hmt2h_dn(ndn,ndn), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?h_dn', 'setup_diagh11')

    allocate (hmt1p_dn(norb-ndn,norb-ndn), &
            & hmt2p_dn(norb-ndn,norb-ndn), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?p_dn', 'setup_diagh11')

  else if ( iwfnty == 3 ) then
    allocate (hmt1h(ntot,ntot), hmt2h(ntot,ntot), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?h', 'setup_diagh11')

    allocate (hmt1p(2*norb-ntot,2*norb-ntot), &
            & hmt2p(2*norb-ntot,2*norb-ntot), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmt?p', 'setup_diagh11')
  end if


  ! allocate space for vec, val

  if ( iwfnty == 1 ) then
    allocate (vech(nup,nup), vecp(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'vec?', 'setup_diagh11')

    allocate (valh(nup), valp(norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'val?', 'setup_diagh11')

  else if ( iwfnty == 2 ) then
    allocate (vech_up(nup,nup), vecp_up(norb-nup,norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'vec?_up', 'setup_diagh11')

    allocate (valh_up(nup), valp_up(norb-nup), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'val?_up', 'setup_diagh11')

    allocate (vech_dn(ndn,ndn), vecp_dn(norb-ndn,norb-ndn), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'vec?_dn', 'setup_diagh11')

    allocate (valh_dn(ndn), valp_dn(norb-ndn), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'val?_dn', 'setup_diagh11')

  else if ( iwfnty == 3 ) then
    allocate (vech(ntot,ntot), vecp(2*norb-ntot,2*norb-ntot), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'vec?', 'setup_diagh11')

    allocate (valh(ntot), valp(2*norb-ntot), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'val?', 'setup_diagh11')
  end if


  ! scratch arrays used in hhmat, ppmat
  ! -----------------------------------

  ! allocate space for scr1, scr2

  if ( iwfnty == 1 ) then
    allocate (scr2(norb,nscr), scr1(norb,nscr), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (scr2(norb,nscr), scr1(norb,nscr), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scr2(2*norb,nscr), scr1(2*norb,nscr), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'setup_diagh11')


  ! scratch arrays used for diagonalization
  ! ---------------------------------------

  if ( iwfnty == 2 ) then
    nscr = max (nscr_up, nscr_dn)
  end if

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', nscr, -1, -1, -1)
  lwork = max (1, (nb+1)*nscr)

  allocate (iwork(5*nscr), ifail(nscr), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'setup_diagh11')

  allocate (work(lwork), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'work', 'setup_diagh11')

  ! allocate (rwork(3*nscr-2), stat=istatus)
  allocate (rwork(7*nscr), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rwork', 'setup_diagh11')


  return
end subroutine setup_diagh11



subroutine shutdown_diagh11 (iwfnty, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_diagh11  --  CAJH, 01.2013                            |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in diagh11_*.    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty  - type of wavefunction to use
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons

  integer, intent(in) :: iwfnty, nup, ndn


  ! other variables

  integer :: istatus


  deallocate (dst, dnost, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dst, dnost', 'shutdown_diagh11')

  deallocate (scrt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scrt', 'shutdown_diagh11')

  deallocate (rho, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rho', 'shutdown_diagh11')

  deallocate (gam, tpgam, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'gam', 'shutdown_diagh11')

  if ( iwfnty == 3 ) then
    deallocate (scr_spb, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'scr_spb', 'shutdown_diagh11')
  end if


  ! shutdown trdnmt scratch

  call shutdown_trdnmt (iwfnty, nup, ndn)

  ! shutdown erictr scratch

  call shutdown_erictr


  ! deallocate space for hmath, hmatp

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (hmath, hmatp, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmat?', 'shutdown_diagh11')
  end if

  if ( iwfnty == 2 ) then
    deallocate (hmath_up, hmatp_up, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmat?_up', 'shutdown_diagh11')

    deallocate (hmath_dn, hmatp_dn, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmat?_dn', 'shutdown_diagh11')
  end if


  ! deallocate space for hmt1, hmt2

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (hmt1h, hmt2h, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?h', 'shutdown_diagh11')

    deallocate (hmt1p, hmt2p, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?p', 'shutdown_diagh11')
  end if

  if ( iwfnty == 2 ) then
    deallocate (hmt1h_up, hmt2h_up, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?h_up', 'shutdown_diagh11')

    deallocate (hmt1p_up, hmt2p_up, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?p_up', 'shutdown_diagh11')

    deallocate (hmt1h_dn, hmt2h_dn, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?h_dn', 'shutdown_diagh11')

    deallocate (hmt1p_dn, hmt2p_dn, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'hmt?p_dn', 'shutdown_diagh11')
  end if


  ! deallocate space for vec, val

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (vech, vecp, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'vec?', 'shutdown_diagh11')

    deallocate (valh, valp, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'val?', 'shutdown_diagh11')
  end if

  if ( iwfnty == 2 ) then
    deallocate (vech_up, vecp_up, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'vec?_up', 'shutdown_diagh11')

    deallocate (valh_up, valp_up, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'val?_up', 'shutdown_diagh11')

    deallocate (vech_dn, vecp_dn, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'vec?_dn', 'shutdown_diagh11')

    deallocate (valh_dn, valp_dn, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'val?_dn', 'shutdown_diagh11')
  end if


  ! scratch arrays used in hhmat, ppmat
  ! -----------------------------------

  deallocate (scr1, scr2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scr?', 'shutdown_diagh11')


  ! scratch arrays used for diagonalization
  ! ---------------------------------------

  deallocate (iwork, ifail, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'iwork', 'shutdown_diagh11')

  deallocate (work, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'work', 'shutdown_diagh11')

  deallocate (rwork, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rwork', 'shutdown_diagh11')


  return
end subroutine shutdown_diagh11



subroutine diagh11_ghf (norb, nbas, nbct, nup, ndn, darr, enarr, xdim, &
     & edim, ndet, indk, xmat, hmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | diagh11_ghf  --  CAJH, 01.2013                                 |
! |                                                                |
! |                                                                |
! | Diagonalize the H^11 Hamiltonian given a Slater determinant    |
! | of GHF type. Reconstruct the determinant using the             |
! | eigenvectors of H^11.                                          |
! |                                                                |
! | This routine accepts an array of D matrices [darr]. The index  |
! | indk indicates which determinant to use.                       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb  - number of orbitals
  !   nbas  - number of basis functions
  !   nbct  - number of Cartesian basis functions
  !   nup   - number of spin-up electrons
  !   ndn   - number of spin-dn electrons

  integer, intent(in) :: norb, nbas, nbct, nup, ndn

  !   darr  - array of determinants [updated]
  !   enarr - array with orbital energies [updated]
  !   xdim  - leading dimension of array darr
  !   edim  - leading dimension of array enarr
  !   ndet  - total number of determinants stored
  !   indk  - index of determinant to use

  integer, intent(in) :: xdim, edim, ndet, indk

  real(kind=dp), dimension(edim,ndet), intent(inout) :: enarr
  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr

  !   xmat  - transformation matrix [ = S^(-1/2) ]
  !   hmat  - core Hamiltonian matrix ( ort AO basis )

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,norb), intent(in) :: hmat

  !   i2sv  - vector of 2-electron integrals
  !   ni2s  - number of 2-electron integrals stored in i2sv
  !   lsf   - whether stored integrals include symmetry factor

  integer, intent(in) :: ni2s
  logical, intent(in) :: lsf

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: nosq, ntot, nntot
  integer :: k, iflg, info, idum

  integer :: inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  complex(kind=dp) :: ovp, hmovp
  complex(kind=dp) :: ztr, ztr_up, ztr_dn


  ! external functions

  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  ntot = nup + ndn
  nosq = norb*norb

  nntot = 2*norb - ntot


  ! load D, D*

  dst(1:4*nosq) = conjg (darr(1:4*nosq,indk))
  dnost(1:4*nosq) = darr(1:4*nosq,indk)
  
  ! compute overlap [should be 1]
  ! build density matrix

  call trdnmt_ghf (norb, nup, ndn, dst, dst, rho, ovp, iflg)

  if ( iflg /= 0 ) then
    write (6, *) 'error: trdnmt_ghf failed in diagh11_ghf.'
    stop
  end if

  ! contract density matrix with two-particle integrals

  call erictr_ghf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

  ! compute energy

  hmovp = z0

  ztr_up = trab_zmat (norb, hmat, norb, rho(1,1), 2*norb, idum)
  ztr_dn = trab_zmat (norb, hmat, norb, rho(norb+1,norb+1), 2*norb, idum)
  hmovp = hmovp + ztr_up + ztr_dn

  ztr = trab_zmat (2*norb, rho, 2*norb, gam, 2*norb, idum)
  hmovp = hmovp + 0.50e0_dp * ztr


  ! build tpgam

  tpgam(1:2*norb,1:2*norb) = gam(1:2*norb,1:2*norb)

  tpgam(1:norb,1:norb) &
       & = tpgam(1:norb,1:norb) &
       & + hmat(1:norb,1:norb)
  tpgam(norb+1:2*norb,norb+1:2*norb) &
       & = tpgam(norb+1:2*norb,norb+1:2*norb) &
       & + hmat(1:norb,1:norb)


  ! build hmt1h, hmt2h
  ! build hmt1p, hmt2p

  call hhmat_ghf (norb, nup, ndn)
  call ppmat_ghf (norb, nup, ndn)


  ! construct hmath, hmatp

  hmath(1:ntot,1:ntot) = &
       & hmovp * hmt1h(1:ntot,1:ntot) - hmt2h(1:ntot,1:ntot)

  hmatp(1:nntot,1:nntot) = &
       & hmovp * hmt1p(1:nntot,1:nntot) + hmt2p(1:nntot,1:nntot)

  do k = 1, ntot
    hmath(k,k) = hmath(k,k) - hmovp
  end do

  do k = 1, nntot
    hmatp(k,k) = hmatp(k,k) - hmovp
  end do


  ! diagonalize hmath, hmatp
  ! set hmat to negative to recover negative MO energies

  call zmat_symm (1, 'u', ntot,  hmath)
  call zmat_symm (1, 'u', nntot, hmatp)

  ! vech(1:ntot,1:ntot)   = -hmath(1:ntot,1:ntot)
  ! vecp(1:nntot,1:nntot) = hmatp(1:nntot,1:nntot)

  hmath(1:ntot,1:ntot) = -hmath(1:ntot,1:ntot)

  call zheevx ('v', 'a', 'u', ntot, hmath, ntot, ddum, ddum, &
       & idum, idum, abstol, inum, valh, vech, ntot, work, &
       & lwork, rwork, iwork, ifail, info)

  ! call zheev ('v', 'u', ntot, vech, ntot, valh, work, lwork, &
  !      & rwork, info)

  if ( info /= 0 ) then
    write (6, *) 'error: Diagonalization of hmath failed in diagh11_ghf.'
    stop
  end if

  call zheevx ('v', 'a', 'u', nntot, hmatp, nntot, ddum, ddum, &
       & idum, idum, abstol, inum, valp, vecp, nntot, work, &
       & lwork, rwork, iwork, ifail, info)

  ! call zheev ('v', 'u', nntot, vecp, nntot, valp, work, lwork, &
  !      & rwork, info)

  if ( info /= 0 ) then
    write (6, *) 'error: Diagonalization of hmatp failed in diagh11_ghf.'
    stop
  end if


  ! redefine vecp

  vecp(1:nntot,1:nntot) = conjg (vecp(1:nntot,1:nntot))


  ! form new determinant D

  call zgemm ('n', 'n', 2*norb, ntot, ntot, z1, darr(1,indk), &
       & 2*norb, vech, ntot, z0, scrt, 2*norb)

  darr(1:2*norb*ntot,indk) = scrt(1:2*norb*ntot)

  call zgemm ('n', 'n', 2*norb, nntot, nntot, z1, darr(2*norb*ntot+1,indk), &
       & 2*norb, vecp, nntot, z0, scrt, 2*norb)

  darr(2*norb*ntot+1:4*nosq,indk) = scrt(1:2*norb*nntot)


  ! store orbital energies

  enarr(1:ntot,indk)        = valh(1:ntot)
  enarr(ntot+1:2*norb,indk) = valp(1:nntot)


  return
end subroutine diagh11_ghf



subroutine hhmat_ghf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | hhmat_ghf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Build contributions to hole-hole block of H^11 for GHF-type    |
! | determinants.                                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: ntot


  ntot = nup + ndn

  ! compute:  scr1 = rho^T . D(:,1:N)

  call zgemm ('t', 'n', 2*norb, ntot, 2*norb, z1, rho, 2*norb, &
       & dnost, 2*norb, z0, scr1, 2*norb)

  ! compute:  hmt1h = D(:,1:N)! . scr1

  call zgemm ('c', 'n', ntot, ntot, 2*norb, z1, dnost, &
       & 2*norb, scr1, 2*norb, z0, hmt1h, ntot)

  ! compute:  scr2 = (T + Gam)^T . scr1

  call zgemm ('t', 'n', 2*norb, ntot, 2*norb, z1, tpgam, 2*norb, &
       & scr1, 2*norb, z0, scr2, 2*norb)

  ! compute:  scr1 = rho^T . scr2

  call zgemm ('t', 'n', 2*norb, ntot, 2*norb, z1, rho, 2*norb, &
       & scr2, 2*norb, z0, scr1, 2*norb)

  ! compute:  hmt2h = D(:,1:N)! . scr1

  call zgemm ('c', 'n', ntot, ntot, 2*norb, z1, dnost, &
       & 2*norb, scr1, 2*norb, z0, hmt2h, ntot)


  return
end subroutine hhmat_ghf



subroutine ppmat_ghf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | ppmat_ghf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Build contributions to particle-particle block of H^11 for     |
! | GHF-type determinants.                                         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: k, ntot, nntot


  ntot  = nup + ndn
  nntot = 2*norb - ntot

  ! build rho-I; use rho as scratch space

  do k = 1, 2*norb
    rho(k,k) = rho(k,k) - z1
  end do

  ! compute:  scr1 = (I - rho) . D(:,N+1:end)*

  call zgemm ('n', 'n', 2*norb, nntot, 2*norb, -z1, rho, 2*norb, &
       & dst(2*norb*ntot+1), 2*norb, z0, scr1, 2*norb)

  ! compute:  hmt1p = D(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', nntot, nntot, 2*norb, z1, dst(2*norb*ntot+1), &
       & 2*norb, scr1, 2*norb, z0, hmt1p, nntot)

  ! compute:  scr2 = (T + Gam) . scr1

  call zgemm ('n', 'n', 2*norb, nntot, 2*norb, z1, tpgam, 2*norb, &
       & scr1, 2*norb, z0, scr2, 2*norb)

  ! compute:  scr1 = (I - rho) . scr2

  call zgemm ('n', 'n', 2*norb, nntot, 2*norb, -z1, rho, 2*norb, &
       & scr2, 2*norb, z0, scr1, 2*norb)

  ! compute:  hmt2p = D(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', nntot, nntot, 2*norb, z1, dst(2*norb*ntot+1), &
       & 2*norb, scr1, 2*norb, z0, hmt2p, nntot)


  return
end subroutine ppmat_ghf



subroutine diagh11_uhf (norb, nbas, nbct, nup, ndn, darr, enarr, xdim, &
     & edim, ndet, indk, xmat, hmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | diagh11_uhf  --  CAJH, 01.2013                                 |
! |                                                                |
! |                                                                |
! | ( UHF version of diagh11_ghf. )                                |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb  - number of orbitals
  !   nbas  - number of basis functions
  !   nbct  - number of Cartesian basis functions
  !   nup   - number of spin-up electrons
  !   ndn   - number of spin-dn electrons

  integer, intent(in) :: norb, nbas, nbct, nup, ndn

  !   darr  - array of determinants [updated]
  !   enarr - array with orbital energies [updated]
  !   xdim  - leading dimension of array darr
  !   edim  - leading dimension of array enarr
  !   ndet  - total number of determinants stored
  !   indk  - index of determinant to use

  integer, intent(in) :: xdim, edim, ndet, indk

  real(kind=dp), dimension(edim,ndet), intent(inout) :: enarr
  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr

  !   xmat  - transformation matrix [ = S^(-1/2) ]
  !   hmat  - core Hamiltonian matrix ( ort AO basis )

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,norb), intent(in) :: hmat

  !   i2sv  - vector of 2-electron integrals
  !   ni2s  - number of 2-electron integrals stored in i2sv
  !   lsf   - whether stored integrals include symmetry factor

  integer, intent(in) :: ni2s
  logical, intent(in) :: lsf

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: nosq, nnup, nndn
  integer :: k, iflg, info_up, info_dn, idum

  integer :: inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  complex(kind=dp) :: ovp, hmovp
  complex(kind=dp) :: ztr_up, ztr_dn


  ! external functions

  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  nosq = norb*norb

  nnup = norb-nup
  nndn = norb-ndn


  ! load D, D*

  dst(1:2*nosq) = conjg (darr(1:2*nosq,indk))
  dnost(1:2*nosq) = darr(1:2*nosq,indk)
  
  ! compute overlap [should be 1]
  ! build density matrix

  call trdnmt_uhf (norb, nup, ndn, dst, dst, rho, ovp, iflg)

  if ( iflg /= 0 ) then
    write (6, *) 'error: trdnmt_uhf failed in diagh11_uhf.'
    stop
  end if

  ! contract density matrix with two-particle integrals

  call erictr_uhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

  ! compute energy

  hmovp = z0

  ztr_up = trab_zmat (norb, hmat, norb, rho(1,1), norb, idum)
  ztr_dn = trab_zmat (norb, hmat, norb, rho(1,norb+1), norb, idum)
  hmovp = hmovp + ztr_up + ztr_dn

  ztr_up = trab_zmat (norb, rho(1,1), norb, gam(1,1), norb, idum)
  ztr_dn = trab_zmat (norb, rho(1,norb+1), norb, gam(1,norb+1), norb, idum)
  hmovp = hmovp + 0.50e0_dp * (ztr_up + ztr_dn)


  ! build tpgam

  tpgam(1:norb,1:norb) = hmat(1:norb,1:norb) &
                     & + gam(1:norb,1:norb)

  tpgam(1:norb,norb+1:2*norb) = hmat(1:norb,1:norb) &
                            & + gam(1:norb,norb+1:2*norb)


  ! build hmt1h, hmt2h
  ! build hmt1p, hmt2p

  call hhmat_uhf (norb, nup, ndn)
  call ppmat_uhf (norb, nup, ndn)


  ! construct hmath, hmatp

  if ( nup > 0 ) then
  hmath_up(1:nup,1:nup) = &
       & hmovp * hmt1h_up(1:nup,1:nup) - hmt2h_up(1:nup,1:nup)
  end if

  if ( ndn > 0 ) then
  hmath_dn(1:ndn,1:ndn) = &
       & hmovp * hmt1h_dn(1:ndn,1:ndn) - hmt2h_dn(1:ndn,1:ndn)
  end if

  if ( nnup > 0 ) then
  hmatp_up(1:nnup,1:nnup) = &
       & hmovp * hmt1p_up(1:nnup,1:nnup) + hmt2p_up(1:nnup,1:nnup)
  end if

  if ( nndn > 0 ) then
  hmatp_dn(1:nndn,1:nndn) = &
       & hmovp * hmt1p_dn(1:nndn,1:nndn) + hmt2p_dn(1:nndn,1:nndn)
  end if

  do k = 1, nup
    hmath_up(k,k) = hmath_up(k,k) - hmovp
  end do
  do k = 1, ndn
    hmath_dn(k,k) = hmath_dn(k,k) - hmovp
  end do

  do k = 1, nnup
    hmatp_up(k,k) = hmatp_up(k,k) - hmovp
  end do
  do k = 1, nndn
    hmatp_dn(k,k) = hmatp_dn(k,k) - hmovp
  end do


  ! diagonalize hmath, hmatp
  ! set hmat to negative to recover negative MO energies

  if ( nup > 0 )   call zmat_symm (1, 'u', nup,  hmath_up)
  if ( ndn > 0 )   call zmat_symm (1, 'u', ndn,  hmath_dn)
  if ( nnup > 0 )  call zmat_symm (1, 'u', nnup, hmatp_up)
  if ( nndn > 0 )  call zmat_symm (1, 'u', nndn, hmatp_dn)

  ! if ( nup > 0 )   vech_up(1:nup,1:nup)   = -hmath_up(1:nup,1:nup)
  ! if ( ndn > 0 )   vech_dn(1:ndn,1:ndn)   = -hmath_dn(1:ndn,1:ndn)
  ! if ( nnup > 0 )  vecp_up(1:nnup,1:nnup) = hmatp_up(1:nnup,1:nnup)
  ! if ( nndn > 0 )  vecp_dn(1:nndn,1:nndn) = hmatp_dn(1:nndn,1:nndn)

  if ( nup > 0 )  hmath_up(1:nup,1:nup) = -hmath_up(1:nup,1:nup)
  if ( ndn > 0 )  hmath_dn(1:ndn,1:ndn) = -hmath_dn(1:ndn,1:ndn)

  info_up = 0;  info_dn = 0

  if ( nup > 0 ) then
  call zheevx ('v', 'a', 'u', nup, hmath_up, nup, ddum, ddum, &
       & idum, idum, abstol, inum, valh_up, vech_up, nup, work, &
       & lwork, rwork, iwork, ifail, info_up)

  ! call zheev ('v', 'u', nup, vech_up, nup, valh_up, work, lwork, &
  !      & rwork, info_up)
  end if

  if ( ndn > 0 ) then
  call zheevx ('v', 'a', 'u', ndn, hmath_dn, ndn, ddum, ddum, &
       & idum, idum, abstol, inum, valh_dn, vech_dn, ndn, work, &
       & lwork, rwork, iwork, ifail, info_dn)

  ! call zheev ('v', 'u', ndn, vech_dn, ndn, valh_dn, work, lwork, &
  !      & rwork, info_dn)
  end if

  if ( info_up /= 0 .or. info_dn /= 0 ) then
    write (6, *) 'error: Diagonalization of hmath failed in diagh11_uhf.'
    stop
  end if

  if ( nnup > 0 ) then
  call zheevx ('v', 'a', 'u', nnup, hmatp_up, nnup, ddum, ddum, &
       & idum, idum, abstol, inum, valp_up, vecp_up, nnup, work, &
       & lwork, rwork, iwork, ifail, info_up)

  ! call zheev ('v', 'u', nnup, vecp_up, nnup, valp_up, work, lwork, &
  !      & rwork, info_up)
  end if

  if ( nndn > 0 ) then
  call zheevx ('v', 'a', 'u', nndn, hmatp_dn, nndn, ddum, ddum, &
       & idum, idum, abstol, inum, valp_dn, vecp_dn, nndn, work, &
       & lwork, rwork, iwork, ifail, info_dn)

  ! call zheev ('v', 'u', nndn, vecp_dn, nndn, valp_dn, work, lwork, &
  !      & rwork, info_dn)
  end if

  if ( info_up /= 0 .or. info_dn /= 0 ) then
    write (6, *) 'error: Diagonalization of hmatp failed in diagh11_uhf.'
    stop
  end if


  ! redefine vecp

  if ( nnup > 0 )  vecp_up(1:nnup,1:nnup) = conjg (vecp_up(1:nnup,1:nnup))
  if ( nndn > 0 )  vecp_dn(1:nndn,1:nndn) = conjg (vecp_dn(1:nndn,1:nndn))


  ! form new determinant D

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nup, nup, z1, darr(1,indk), &
       & norb, vech_up, nup, z0, scrt, norb)

  darr(1:norb*nup,indk) = scrt(1:norb*nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, ndn, ndn, z1, darr(nosq+1,indk), &
       & norb, vech_dn, ndn, z0, scrt, norb)

  darr(nosq+1:nosq+norb*ndn,indk) = scrt(1:norb*ndn)
  end if

  if ( nnup > 0 ) then
  call zgemm ('n', 'n', norb, nnup, nnup, z1, darr(norb*nup+1,indk), &
       & norb, vecp_up, nnup, z0, scrt, norb)

  darr(norb*nup+1:nosq,indk) = scrt(1:norb*nnup)
  end if

  if ( nndn > 0 ) then
  call zgemm ('n', 'n', norb, nndn, nndn, z1, darr(nosq+norb*ndn+1,indk), &
       & norb, vecp_dn, nndn, z0, scrt, norb)

  darr(nosq+norb*ndn+1:2*nosq,indk) = scrt(1:norb*nndn)
  end if


  ! store orbital energies

  if ( nup > 0 )   enarr(1:nup,indk)             = valh_up(1:nup)
  if ( ndn > 0 )   enarr(norb+1:norb+ndn,indk)   = valh_dn(1:ndn)
  if ( nnup > 0 )  enarr(nup+1:norb,indk)        = valp_up(1:nnup)
  if ( nndn > 0 )  enarr(norb+ndn+1:2*norb,indk) = valp_dn(1:nndn)


  return
end subroutine diagh11_uhf



subroutine hhmat_uhf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | hhmat_uhf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | ( UHF version of hhmat_ghf. )                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: nosq


  nosq = norb*norb

  ! compute:  scr1 = rho^T . D(:,1:N)

  if ( nup > 0 ) then
  call zgemm ('t', 'n', norb, nup, norb, z1, rho(1,1), norb, &
       & dnost, norb, z0, scr1(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('t', 'n', norb, ndn, norb, z1, rho(1,norb+1), norb, &
       & dnost(nosq+1), norb, z0, scr1(1,nup+1), norb)
  end if

  ! compute:  hmt1h = D(:,1:N)! . scr1

  if ( nup > 0 ) then
  call zgemm ('c', 'n', nup, nup, norb, z1, dnost, &
       & norb, scr1(1,1), norb, z0, hmt1h_up, nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('c', 'n', ndn, ndn, norb, z1, dnost(nosq+1), &
       & norb, scr1(1,nup+1), norb, z0, hmt1h_dn, ndn)
  end if

  ! compute:  scr2 = (T + Gam)^T . scr1

  if ( nup > 0 ) then
  call zgemm ('t', 'n', norb, nup, norb, z1, tpgam(1,1), norb, &
       & scr1(1,1), norb, z0, scr2(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('t', 'n', norb, ndn, norb, z1, tpgam(1,norb+1), norb, &
       & scr1(1,nup+1), norb, z0, scr2(1,nup+1), norb)
  end if

  ! compute:  scr1 = rho^T . scr2

  if ( nup > 0 ) then
  call zgemm ('t', 'n', norb, nup, norb, z1, rho(1,1), norb, &
       & scr2(1,1), norb, z0, scr1(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('t', 'n', norb, ndn, norb, z1, rho(1,norb+1), norb, &
       & scr2(1,nup+1), norb, z0, scr1(1,nup+1), norb)
  end if

  ! compute:  hmt2h = D(:,1:N)! . scr1

  if ( nup > 0 ) then
  call zgemm ('c', 'n', nup, nup, norb, z1, dnost, &
       & norb, scr1(1,1), norb, z0, hmt2h_up, nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('c', 'n', ndn, ndn, norb, z1, dnost(nosq+1), &
       & norb, scr1(1,nup+1), norb, z0, hmt2h_dn, ndn)
  end if


  return
end subroutine hhmat_uhf



subroutine ppmat_uhf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | ppmat_uhf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | ( UHF version of ppmat_ghf. )                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: k, nosq, nnup, nndn


  nosq = norb*norb
  nnup = norb-nup
  nndn = norb-ndn

  ! build rho-I; use rho as scratch space

  do k = 1, norb
    rho(k,k) = rho(k,k) - z1
  end do

  do k = 1, norb
    rho(k,k+norb) = rho(k,k+norb) - z1
  end do

  ! compute:  scr1 = (I - rho) . D(:,N+1:end)*

  if ( nnup > 0 ) then
  call zgemm ('n', 'n', norb, nnup, norb, -z1, rho(1,1), norb, &
       & dst(norb*nup+1), norb, z0, scr1(1,1), norb)
  end if

  if ( nndn > 0 ) then
  call zgemm ('n', 'n', norb, nndn, norb, -z1, rho(1,norb+1), norb, &
       & dst(nosq+norb*ndn+1), norb, z0, scr1(1,nnup+1), norb)
  end if

  ! compute:  hmt1p = D(:,N+1:end)^T . scr1

  if ( nnup > 0 ) then
  call zgemm ('c', 'n', nnup, nnup, norb, z1, dst(norb*nup+1), &
       & norb, scr1(1,1), norb, z0, hmt1p_up, nnup)
  end if

  if ( nndn > 0 ) then
  call zgemm ('c', 'n', nndn, nndn, norb, z1, dst(nosq+norb*ndn+1), &
       & norb, scr1(1,nnup+1), norb, z0, hmt1p_dn, nndn)
  end if

  ! compute:  scr2 = (T + Gam) . scr1

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nnup, norb, z1, tpgam(1,1), norb, &
       & scr1(1,1), norb, z0, scr2(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, nndn, norb, z1, tpgam(1,norb+1), norb, &
       & scr1(1,nnup+1), norb, z0, scr2(1,nnup+1), norb)
  end if

  ! compute:  scr1 = (I - rho) . scr2

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nnup, norb, -z1, rho(1,1), norb, &
       & scr2(1,1), norb, z0, scr1(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, nndn, norb, -z1, rho(1,norb+1), norb, &
       & scr2(1,nnup+1), norb, z0, scr1(1,nnup+1), norb)
  end if

  ! compute:  hmt2p = D(:,N+1:end)^T . scr1

  if ( nnup > 0 ) then
  call zgemm ('c', 'n', nnup, nnup, norb, z1, dst(norb*nup+1), &
       & norb, scr1(1,1), norb, z0, hmt2p_up, nnup)
  end if

  if ( nndn > 0 ) then
  call zgemm ('c', 'n', nndn, nndn, norb, z1, dst(nosq+norb*ndn+1), &
       & norb, scr1(1,nnup+1), norb, z0, hmt2p_dn, nndn)
  end if


  return
end subroutine ppmat_uhf



subroutine diagh11_rhf (norb, nbas, nbct, nup, darr, enarr, xdim, &
     & edim, ndet, indk, xmat, hmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | diagh11_rhf  --  CAJH, 01.2013                                 |
! |                                                                |
! |                                                                |
! | ( RHF version of diagh11_ghf. )                                |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb  - number of orbitals
  !   nbas  - number of basis functions
  !   nbct  - number of Cartesian basis functions
  !   nup   - number of spin-up electrons

  integer, intent(in) :: norb, nbas, nbct, nup

  !   darr  - array of determinants [updated]
  !   enarr - array with orbital energies [updated]
  !   xdim  - leading dimension of array darr
  !   edim  - leading dimension of array enarr
  !   ndet  - total number of determinants stored
  !   indk  - index of determinant to use

  integer, intent(in) :: xdim, edim, ndet, indk

  real(kind=dp), dimension(edim,ndet), intent(inout) :: enarr
  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr

  !   xmat  - transformation matrix [ = S^(-1/2) ]
  !   hmat  - core Hamiltonian matrix ( ort AO basis )

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,norb), intent(in) :: hmat

  !   i2sv  - vector of 2-electron integrals
  !   ni2s  - number of 2-electron integrals stored in i2sv
  !   lsf   - whether stored integrals include symmetry factor

  integer, intent(in) :: ni2s
  logical, intent(in) :: lsf

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: nosq, nnup
  integer :: k, iflg, info, idum

  integer :: inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  complex(kind=dp) :: ovp, hmovp
  complex(kind=dp) :: ztr


  ! external functions

  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  nosq = norb*norb

  nnup = norb-nup


  ! load D, D*

  dst(1:nosq) = conjg (darr(1:nosq,indk))
  dnost(1:nosq) = darr(1:nosq,indk)
  
  ! compute overlap [should be 1]
  ! build density matrix

  call trdnmt_rhf (norb, nup, dst, dst, rho, ovp, iflg)

  if ( iflg /= 0 ) then
    write (6, *) 'error: trdnmt_rhf failed in diagh11_rhf.'
    stop
  end if

  ! contract density matrix with two-particle integrals

  call erictr_rhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

  ! compute energy

  hmovp = z0

  ztr = trab_zmat (norb, hmat, norb, rho, norb, idum)
  hmovp = hmovp + d2*ztr

  ztr = trab_zmat (norb, rho, norb, gam, norb, idum)
  hmovp = hmovp + ztr


  ! build tpgam

  tpgam(1:norb,1:norb) = hmat(1:norb,1:norb) &
                     & + gam(1:norb,1:norb)


  ! build hmt1h, hmt2h
  ! build hmt1p, hmt2p

  call hhmat_rhf (norb, nup)
  call ppmat_rhf (norb, nup)


  ! construct hmath, hmatp

  hmath(1:nup,1:nup) = &
       & hmovp * hmt1h(1:nup,1:nup) - hmt2h(1:nup,1:nup)

  hmatp(1:nnup,1:nnup) = &
       & hmovp * hmt1p(1:nnup,1:nnup) + hmt2p(1:nnup,1:nnup)

  do k = 1, nup
    hmath(k,k) = hmath(k,k) - hmovp
  end do

  do k = 1, nnup
    hmatp(k,k) = hmatp(k,k) - hmovp
  end do


  ! diagonalize hmath, hmatp
  ! set hmat to negative to recover negative MO energies

  call zmat_symm (1, 'u', nup,  hmath)
  call zmat_symm (1, 'u', nnup, hmatp)

  ! vech(1:nup,1:nup)   = -hmath(1:nup,1:nup)
  ! vecp(1:nnup,1:nnup) = hmatp(1:nnup,1:nnup)

  hmath(1:nup,1:nup) = -hmath(1:nup,1:nup)

  call zheevx ('v', 'a', 'u', nup, hmath, nup, ddum, ddum, &
       & idum, idum, abstol, inum, valh, vech, nup, work, &
       & lwork, rwork, iwork, ifail, info)

  ! call zheev ('v', 'u', nup, vech, nup, valh, work, lwork, &
  !      & rwork, info)

  if ( info /= 0 ) then
    write (6, *) 'error: Diagonalization of hmath failed in diagh11_rhf.'
    stop
  end if

  call zheevx ('v', 'a', 'u', nnup, hmatp, nnup, ddum, ddum, &
       & idum, idum, abstol, inum, valp, vecp, nnup, work, &
       & lwork, rwork, iwork, ifail, info)

  ! call zheev ('v', 'u', nnup, vecp, nnup, valp, work, lwork, &
  !      & rwork, info)

  if ( info /= 0 ) then
    write (6, *) 'error: Diagonalization of hmatp failed in diagh11_rhf.'
    stop
  end if


  ! redefine vecp

  vecp(1:nnup,1:nnup) = conjg (vecp(1:nnup,1:nnup))


  ! form new determinant D

  call zgemm ('n', 'n', norb, nup, nup, z1, darr(1,indk), &
       & norb, vech, nup, z0, scrt, norb)

  darr(1:norb*nup,indk) = scrt(1:norb*nup)

  call zgemm ('n', 'n', norb, nnup, nnup, z1, darr(norb*nup+1,indk), &
       & norb, vecp, nnup, z0, scrt, norb)

  darr(norb*nup+1:nosq,indk) = scrt(1:norb*nnup)


  ! store orbital energies

  enarr(1:nup,indk)      = valh(1:nup)
  enarr(nup+1:norb,indk) = valp(1:nnup)


  return
end subroutine diagh11_rhf



subroutine hhmat_rhf (norb, nup)

! +----------------------------------------------------------------+
! |                                                                |
! | hhmat_rhf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | ( RHF version of hhmat_ghf. )                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons

  integer, intent(in) :: norb, nup


  ! compute:  scr1 = rho^T . D(:,1:N)

  call zgemm ('t', 'n', norb, nup, norb, z1, rho, norb, &
       & dnost, norb, z0, scr1, norb)

  ! compute:  hmt1h = D(:,1:N)! . scr1

  call zgemm ('c', 'n', nup, nup, norb, z1, dnost, &
       & norb, scr1, norb, z0, hmt1h, nup)

  ! compute:  scr2 = (T + Gam)^T . scr1

  call zgemm ('t', 'n', norb, nup, norb, z1, tpgam, norb, &
       & scr1, norb, z0, scr2, norb)

  ! compute:  scr1 = rho^T . scr2

  call zgemm ('t', 'n', norb, nup, norb, z1, rho, norb, &
       & scr2, norb, z0, scr1, norb)

  ! compute:  hmt2h = D(:,1:N)! . scr1

  call zgemm ('c', 'n', nup, nup, norb, z1, dnost, &
       & norb, scr1, norb, z0, hmt2h, nup)


  return
end subroutine hhmat_rhf



subroutine ppmat_rhf (norb, nup)

! +----------------------------------------------------------------+
! |                                                                |
! | ppmat_rhf  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | ( RHF version of ppmat_ghf. )                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons

  integer, intent(in) :: norb, nup


  ! other variables

  integer :: k, nnup


  nnup = norb-nup

  ! build rho-I; use rho as scratch space

  do k = 1, norb
    rho(k,k) = rho(k,k) - z1
  end do

  ! compute:  scr1 = (I - rho) . D(:,N+1:end)*

  call zgemm ('n', 'n', norb, nnup, norb, -z1, rho, norb, &
       & dst(norb*nup+1), norb, z0, scr1, norb)

  ! compute:  hmt1p = D(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', nnup, nnup, norb, z1, dst(norb*nup+1), &
       & norb, scr1, norb, z0, hmt1p, nnup)

  ! compute:  scr2 = (T + Gam) . scr1

  call zgemm ('n', 'n', norb, nnup, norb, z1, tpgam, norb, &
       & scr1, norb, z0, scr2, norb)

  ! compute:  scr1 = (I - rho) . scr2

  call zgemm ('n', 'n', norb, nnup, norb, -z1, rho, norb, &
       & scr2, norb, z0, scr1, norb)

  ! compute:  hmt2p = D(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', nnup, nnup, norb, z1, dst(norb*nup+1), &
       & norb, scr1, norb, z0, hmt2p, nnup)


  return
end subroutine ppmat_rhf


end module diagh11


