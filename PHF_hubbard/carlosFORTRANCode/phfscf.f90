

module phfscf

  use constants
  use util

  implicit none
  private
  save

  public :: scfhf_lbfgs

! +----------------------------------------------------------------+
! |                                                                |
! | phfscf                                                         |
! |                                                                |
! |                                                                |
! | Driver module to handle SCF optimizations of Hartree-Fock      |
! | or projected Hartree-Fock wavefunctions.                       |
! |                                                                |
! +----------------------------------------------------------------+


  ! scratch arrays

  complex(kind=dp), dimension(:), allocatable :: dscr, zscr
  complex(kind=dp), dimension(:), allocatable :: zscr_up, zscr_dn

  complex(kind=dp), dimension(:), allocatable :: scrh, scrp
  complex(kind=dp), dimension(:), allocatable :: scrh_up, scrp_up
  complex(kind=dp), dimension(:), allocatable :: scrh_dn, scrp_dn

  complex(kind=dp), dimension(:,:), allocatable :: lmat, mmat
  complex(kind=dp), dimension(:,:), allocatable :: lmat_up, mmat_up
  complex(kind=dp), dimension(:,:), allocatable :: lmat_dn, mmat_dn


contains


subroutine scfhf_lbfgs (root, procid, nproc, iout, outfile, iprint, imethd, &
     & iwfnty, norb, nbas, nbct, nup, ndn, xmat, hmat, smat, i2sv, ni2s, lsf, &
     & enr, imult, ngrda, ngrdb, ngrdg, igrda, igrdb, igrdg, sptgr, sptirr, &
     & frot, ndspin, ndsptl, ndcplx, darr0, enarr, xdim, edim, ndet, ndtt, &
     & ovdet, hmdet, nstat, idetst, mxdtst, ovstat, hmstat, stvec, stval, &
     & mdvec, iopt, ndtopt, maxit, mxfunc, nbuf, gtol, nflag)

  use i2sint
  use grid
  use spatpr
  use linalg
  use lbfgs
  use diagh11
  use hamilt

! +----------------------------------------------------------------+
! |                                                                |
! | scfhf_lbfgs  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Driver routine for SCF (self-consistent field) optimizations   |
! | of Hartree-Fock or projected Hartree-Fock wavefunctions using  |
! | a limited-memory Quasi-Newton (BFGS) method.                   |
! |                                                                |
! | Given a reference state expressed as a linear combination of   |
! | Slater determinants {|phi_j>}, we parametrize the energy in    |
! | terms of the Thouless coefficients Z_{ph,j} as                 |
! |                                                                |
! |           sum_{j,k} c_j* c_k <phi_j | e^Z!_j H e^Z_k | phi_k>  |
! |   E[Z] =  ---------------------------------------------------  |
! |            sum_{j,k} c_j* c_k <phi_j | e^Z!_j e^Z_k | phi_k>   |
! |                                                                |
! | where Z_j = sum_{ph} Z_{ph,j} b!_p b_h is a Thouless rotation  |
! | operator written in terms of hole-annihilation (b_h) and       |
! | particle-creation (b!_p) operators acting on the reference     |
! | determinant |phi_j>.                                           |
! |                                                                |
! | We start the optimization by setting Z=0. The optimization is  |
! | handled with the LBFGS Quasi-Newton method.                    |
! |                                                                |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   nflag  =  0,  optimization completed                         |
! |          =  1,  optimization failed - LBGFS error              |
! |          =  2,  maximum number of iterations exceeded          |
! |          =  3,  maximum number of fcn evals exceeded           |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The Thouless theorem states that                               |
! |                                                                |
! | Given a Slater determinant |phi> which is a vacuum to the      |
! | operators {b!_h, b_p}, any Slater determinant |phi'> which is  |
! | not orthogonal to |phi> can be written as                      |
! |                                                                |
! |   |phi'>  =  N . exp [sum_ph Z_{ph} b!_p b_h] |phi>,           |
! |                                                                |
! | where N is a normalization constant and the coefficients Z_ph  |
! | are uniquely determined. Conversely, any wavefunction of the   |
! | form of the r.h.s. is also a Slater determinant.               |
! |                                                                |
! | Hence, at any point during the optimization, we construct an   |
! | updated set of Slater determinants by using Z and the set of   |
! | reference states {|phi_j>}.                                    |
! |                                                                |
! +----------------------------------------------------------------+

  include 'mpif.h'

  ! input / output variables

  !   root    - id of root process
  !   procid  - id of each process
  !   nproc   - total number of processes running
  !   iout    - output file unit number
  !   outfile - output file name
  !   iprint  - printing level

  integer, intent(in) :: root, procid, nproc, iout, iprint

  character(len=*), intent(in) :: outfile

  !   imethd  - PHF method to use
  !   iwfnty  - type of wavefunction to use
  !   norb    - number of orbitals
  !   nbas    - number of basis functions
  !   nbct    - number of Cartesian basis functions
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons
  !   xmat    - transformation matrix [ = S^(-1/2) ]
  !   hmat    - core Hamiltonian matrix ( std AO basis )
  !   smat    - overlap matrix ( std AO basis )

  integer, intent(in) :: imethd, iwfnty, norb, nbas, nbct, nup, ndn

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: hmat, smat

  !   i2sv    - vector of 2-electron integrals
  !   ni2s    - number of 2-electron integrals stored in i2sv
  !   lsf     - whether stored integrals include symmetry factor
  !   enr     - nuclear repulsion energy

  integer, intent(in) :: ni2s
  logical, intent(in) :: lsf

  real(kind=dp) :: enr

  type (i2s), dimension(ni2s), intent(in) :: i2sv

  !   imult   - multiplicity to use in spin projection
  !   ngrda   - number of grid points (alpha)
  !   ngrdb   - number of grid points (beta )
  !   ngrdg   - number of grid points (gamma)
  !   igrda   - quadrature selection (alpha integration)
  !   igrdb   - quadrature selection (beta  integration)
  !   igrdg   - quadrature selection (gamma integration)
  !   sptgr   - point group label for spatial projection
  !   sptirr  - irrep to use in spatial projection
  !   frot    - name of file with rotation matrices

  integer, intent(in) :: imult, sptirr
  integer, intent(in) :: ngrda, ngrdb, ngrdg
  integer, intent(in) :: igrda, igrdb, igrdg

  character(len=4), intent(in) :: sptgr
  character(len=*), intent(in) :: frot

  !   ndspin  - number of configs per bs determinant (spin)
  !   ndsptl  - number of configs per bs determinant (spatial)
  !   ndcplx  - number of configs per bs determinant (complex)

  integer, intent(in) :: ndspin, ndsptl, ndcplx

  !   darr0   - array of determinants [updated]
  !   enarr   - array with orbital energies [updated]
  !   xdim    - leading dimension of array darr0
  !   edim    - leading dimension of array enarr
  !   ndet    - total number of determinants stored
  !   ndtt    - dimension of ovdet, hmdet
  !             (ndtt = ndet*ndspin*ndsptl*ndcplx)
  !   ovdet   - overlap matrix among dets [updated]
  !   hmdet   - hamiltonian matrix among dets [updated]

  integer, intent(in) :: xdim, edim, ndet, ndtt

  real(kind=dp), dimension(edim,ndet), intent(inout) :: enarr
  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr0
  complex(kind=dp), dimension(ndtt,ndtt), intent(inout) :: ovdet, hmdet

  !   nstat   - number of states in wavefunction
  !   idetst  - array with number of determinants per state
  !   mxdtst  - maximum number of determinants in a state
  !   ovstat  - overlap matrix among states [updated]
  !   hmstat  - hamiltonian matrix among states [updated]
  !   stvec   - generalized eigenvectors of hmstat [updated]
  !   stval   - generalized eigenvalues of hmstat [updated]
  !   mdvec   - CI expansion of states among its dets [updated]

  integer, intent(in) :: nstat, mxdtst
  integer, dimension(nstat), intent(in) :: idetst

  real(kind=dp), dimension(nstat), intent(inout) :: stval

  complex(kind=dp), dimension(nstat,nstat), intent(inout) :: stvec
  complex(kind=dp), dimension(nstat,nstat), intent(inout) :: ovstat, hmstat
  complex(kind=dp), dimension(ndspin*ndsptl*ndcplx*mxdtst,nstat), &
       & intent(inout) :: mdvec

  !   iopt    - type of optimization (1 - FED, 2 - RES)
  !   ndtopt  - number of determinants to update

  integer, intent(in) :: iopt, ndtopt

  !   maxit   - maximum number of iterations allowed
  !   mxfunc  - maxumum number of function evaluations allowed
  !   nbuf    - number of BFGS updates to use in L-BFGS
  !   gtol    - convergence criterion (norm of gradient)
  !   nflag   - return code (see above)

  integer, intent(in) :: maxit, mxfunc, nbuf
  integer, intent(out) :: nflag

  real(kind=dp), intent(in) :: gtol


  ! spatial symmetry projection

  complex(kind=dp), dimension(:,:,:), allocatable :: rotmat, wgtmat

  ! spin integration grid

  real(kind=dp), dimension(:), allocatable :: grda, grdb, grdg
  real(kind=dp), dimension(:), allocatable :: wgta, wgtb, wgtg

  ! lbfgs variables

  integer :: iter, nfunc, iflag

  integer, dimension(3) :: iprtv
  integer, dimension(100) :: ibuf

  real(kind=dp), dimension(:), allocatable :: zvec, gvec, diag


  ! other variables

  integer :: zdim, ldim, mdim, nop, ngrdx
  integer :: j, ik, ipass, ipos, ierr, iflg, istatus
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: mint, mdc
  parameter ( mint = mpi_integer, mdc = mpi_double_complex )

  integer :: mworld
  parameter ( mworld = mpi_comm_world )

  real(kind=dp) :: energy

  complex(kind=dp), dimension(:), allocatable :: zarr, larr, marr, garr
  complex(kind=dp), dimension(:,:), allocatable :: darr1

  complex(kind=dp), dimension(:,:), allocatable :: inovst, inhmst

  complex(kind=dp), dimension(:,:), allocatable :: hmatt, scrt


  ! constants

  real(kind=dp) :: pi
  parameter ( pi = 4.0e0_dp * atan(d1) )

  !   thrsh - threshold to use in determining linear independency

  real(kind=dp) :: thrsh
  parameter ( thrsh = 1.0e-8_dp )


  nflag = 0

  lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
  lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
  lspin = mod(imethd-1,2) == 1

  lsuhf = .false.
  lsghf = .false.

  if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
  if ( lspin .and. iwfnty == 3 )  lsghf = .true.

    ik = 0
  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do

  if ( iopt == 1 ) then          ! for FED optimizations, only the last
    ik = ik + idetst(nstat)-1    ! determinant is updated
  end if


  ! transform core Hamiltonian to ort AO basis

  allocate (hmatt(norb,norb), scrt(nbas,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'hmatt', 'scfhf_lbfgs')

  call zgemm ('n', 'n', nbas, norb, nbas, z1, hmat, nbas, &
       & xmat, nbas, z0, scrt, nbas)
  call zgemm ('c', 'n', norb, norb, nbas, z1, xmat, nbas, &
       & scrt, nbas, z0, hmatt, norb)


  ! determine dimension of zdim, ldim, mdim

  !   zdim - dimension of Thouless matrices
  !   ldim - dimension of L matrix (nh x nh)
  !   mdim - dimension of M matrix (np x np)

  if ( iwfnty == 1 ) then
    zdim = nup*(norb-nup)
    ldim = nup*nup
    mdim = (norb-nup)*(norb-nup)

  else if ( iwfnty == 2 ) then
    zdim = nup*(norb-nup) + ndn*(norb-ndn)
    ldim = nup*nup + ndn*ndn
    mdim = (norb-nup)*(norb-nup) + (norb-ndn)*(norb-ndn)

  else if ( iwfnty == 3 ) then
    zdim = (nup+ndn)*(2*norb-nup-ndn)
    ldim = (nup+ndn)*(nup+ndn)
    mdim = (2*norb-nup-ndn)*(2*norb-nup-ndn)
  end if


  ! allocate space for zarr, larr, marr, garr

  if ( procid == root ) then
    allocate (zarr(zdim*ndtopt), larr(ldim*ndtopt), &
            & marr(mdim*ndtopt), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, '?arr', 'scfhf_lbfgs')
  end if

  allocate (garr(zdim*ndtopt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'garr', 'scfhf_lbfgs')


  ! prepare spatial symmetry projection

  if ( lsptl ) then

    if ( procid == root ) then
      select case (sptgr)
        case ('cs  ')
          call spatpr_cs  (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c2  ')
          call spatpr_c2  (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c2h ')
          call spatpr_c2h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c2v ')
          call spatpr_c2v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c4v ')
          call spatpr_c4v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c6v ')
          call spatpr_c6v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c8v ')
          call spatpr_c8v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c12v')
          call spatpr_c12v(sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('c16v')
          call spatpr_c16v(sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d2  ')
          call spatpr_d2  (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d2h ')
          call spatpr_d2h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d4h ')
          call spatpr_d4h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d6h ')
          call spatpr_d6h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d8h ')
          call spatpr_d8h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d12h')
          call spatpr_d12h(sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case ('d16h')
          call spatpr_d16h(sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
               & nbas, nbct, xmat, smat, frot)

        case default
          write (6, *) 'error: Point group not supported in phfscf.'
          stop
      end select
    end if

    call mpi_bcast (nop, 1, mint, root, mworld, ierr)

    if ( procid /= root ) then
      allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'scfhf_lbfgs')

      allocate (rotmat(norb,norb,nop), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'scfhf_lbfgs')
    end if

    call mpi_bcast (rotmat, norb*norb*nop, mdc, root, mworld, ierr)
    call mpi_bcast (wgtmat, ndsptl*ndsptl*nop, mdc, root, mworld, ierr)

  else
    nop = 1

    allocate (wgtmat(1,1,1), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'scfhf_lbfgs')

    allocate (rotmat(norb,norb,1), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'scfhf_lbfgs')

    wgtmat(1,1,1) = z0
    rotmat(1:norb,1:norb,1) = z0
  end if


  ! prepare spin projection grid

  allocate (grda(ngrda), grdb(ngrdb), grdg(ngrdg), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'grd?', 'scfhf_lbfgs')

  allocate (wgta(ngrda), wgtb(ngrdb), wgtg(ngrdg), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgt?', 'scfhf_lbfgs')

  if ( ngrda > 1 .and. igrda == 1 ) then
    call trpgrd (d0, d2*pi, ngrda, grda(1:ngrda), wgta(1:ngrda))

  else if ( ngrda > 1 .and. igrda == 2 ) then
    call gauleg (d0, d2*pi, ngrda, grda(1:ngrda), wgta(1:ngrda))

  else if ( ngrda == 1 ) then
    grda(:) = d0
    wgta(:) = d1
  end if

    if ( ngrda > 1 ) then
      wgta(1:ngrda) = wgta(1:ngrda) / (d2 * pi)
    end if

  if ( ngrdb > 1 .and. igrdb == 1 ) then
    call trpgrd (d0, pi, ngrdb, grdb(1:ngrdb), wgtb(1:ngrdb))

  else if ( ngrdb > 1 .and. igrdb == 2 ) then
    call gauleg (d0, pi, ngrdb, grdb(1:ngrdb), wgtb(1:ngrdb))

  else if ( ngrdb == 1 ) then
    grdb(:) = d0
    wgtb(:) = d1
  end if

    if ( ngrdb > 1 ) then    ! account for sin(beta) in integral
      do j = 1, ngrdb
        wgtb(j) = wgtb(j) * sin(grdb(j)) / d2
      end do
    end if

  if ( ngrdg > 1 .and. igrdg == 1 ) then
    call trpgrd (d0, d2*pi, ngrdg, grdg(1:ngrdg), wgtg(1:ngrdg))

  else if ( ngrdg > 1 .and. igrdg == 2 ) then
    call gauleg (d0, d2*pi, ngrdg, grdg(1:ngrdg), wgtg(1:ngrdg))

  else if ( ngrdg == 1 ) then
    grdg(:) = d0
    wgtg(:) = d1
  end if

    if ( ngrdg > 1 ) then
      wgtg(1:ngrdg) = wgtg(1:ngrdg) / (d2 * pi)
    end if

  ngrdx = ngrda*ngrdb*ngrdg*nop


  ! invert ovstat among nstat-1 previous states

  if ( procid == root ) then
    call invert_ovst (nstat, ovstat, hmstat, inovst, inhmst, iflg)

    if ( iflg /= 0 ) then
      write (6, *) 'error: Problems inverting ovstat.'
      stop
    end if
  end if


  ! choose orbitals that diagonalize H^11
  !   ( only determinants set to be optimized )

  if ( procid == root ) then

    call setup_diagh11 (iwfnty, norb, nbas, nbct, nup, ndn)

    do j = ik+1, ndet
      if ( iwfnty == 1 ) then
        call diagh11_rhf (norb, nbas, nbct, nup, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)

      else if ( iwfnty == 2 ) then
        call diagh11_uhf (norb, nbas, nbct, nup, ndn, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)

      else if ( iwfnty == 3 ) then
        call diagh11_ghf (norb, nbas, nbct, nup, ndn, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)
      end if
    end do

    call shutdown_diagh11 (iwfnty, nup, ndn)

  end if


  ! print determinants to output file

  if ( iprint >= 1 .and. procid == root ) then
    open (unit = iout, file = outfile, status = 'old', &
        & form = 'formatted', position = 'append')

    call print_darr_ao (iout, iwfnty, norb, nbas, xmat, darr0, enarr, &
         & xdim, edim, ndet, ik+1, ndet)

    close (unit = iout)
  end if


  ! setup scratch arrays to update orbitals

  if ( procid == root ) then
    call setup_update (norb, nup, ndn, iwfnty)
  end if


  ! setup lbfgs

  if ( procid == root ) then
    allocate (zvec(2*zdim*ndtopt), gvec(2*zdim*ndtopt), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'zvec, gvec', 'scfhf_lbfgs')

    allocate (diag(2*zdim*ndtopt), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'diag', 'scfhf_lbfgs')

    zvec(1:2*zdim*ndtopt) = d0   ! initialize position vector to 0
    iprtv = (/ 1, 0, 1 /)        ! see lbfgs for more details

    call lbfgs_setup (iout, outfile, 2*zdim*ndtopt, nbuf, &
         & .false., iprtv, gtol)
  end if


  ! setup energy evaluation

  call setup_hf (root, procid, nproc, imethd, iwfnty, norb, nbas, nbct, &
       & nup, ndn, ndspin, ndsptl, ndcplx, ndet, ngrdx, iopt, ndtopt, &
       & nstat, idetst, zdim)


  ! +------------------+
  ! |  iteration loop  |
  ! +------------------+

  iflag = 0
  ipass = 0

  allocate (darr1(xdim,ndet), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'darr1', 'scfhf_lbfgs')

  darr1 = darr0   ! array with updated determinants

  do

    ! prepare zarr from zvec

    if ( procid == root .and. ipass /= 0 ) then
      do j = 1, zdim*ndtopt
        zarr(j) = cmplx (zvec(2*j-1), zvec(2*j), dp)
      end do
    end if

    ! update orbitals

    if ( procid == root .and. ipass /= 0 ) then
      if ( iwfnty == 1 ) then
        call update_orb_rhf (norb, nup, darr0, darr1, xdim, ndet, &
             & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
             & mdim, iflg)

      else if ( iwfnty == 2 ) then
        call update_orb_uhf (norb, nup, ndn, darr0, darr1, xdim, ndet, &
             & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
             & mdim, iflg)

      else if ( iwfnty == 3 ) then
        call update_orb_ghf (norb, nup, ndn, darr0, darr1, xdim, ndet, &
             & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
             & mdim, iflg)
      end if

      if ( iflg /= 0 ) then
        write (6, *) 'error: Update of orbitals failed.'
        stop
      end if
    end if


    ! broadcast updated determinants

    call mpi_bcast (darr1(1,ik+1), xdim*ndtopt, mdc, root, mworld, ierr)


    ! compute energy and local gradient

    call hamilt_hf (root, procid, iout, outfile, iprint, imethd, iwfnty, &
         & norb, nbas, nbct, nup, ndn, xmat, hmatt, i2sv, ni2s, lsf, enr, &
         & ndspin, ndsptl, ndcplx, darr1, xdim, ndet, ndtt, ovdet, hmdet, &
         & nstat, idetst, mxdtst, inovst, inhmst, mdvec, iopt, ndtopt, &
         & imult, ngrda, ngrdb, ngrdg, grda, grdb, grdg, wgta, wgtb, wgtg, &
         & nop, rotmat, wgtmat, energy, garr, zdim, thrsh, iflg)

    if ( procid == root ) then
      if ( iflg /= 0 ) then
        write (6, *) 'error: Energy evaluation failed.'
        stop
      end if
    end if


    ! transform local gradient into global gradient

    if ( procid == root ) then
      if ( ipass /= 0 ) then
        call build_global_grad (norb, nup, ndn, iwfnty, garr, zdim, &
             & ndtopt, larr, marr, ldim, mdim)
      end if

      do j = 1, zdim*ndtopt
        gvec(2*j-1) = real  (garr(j))
        gvec(2*j)   = aimag (garr(j))
      end do
    end if


    ! handle control to L-BFGS

    if ( procid == root ) then
      call lbfgs_drv (zvec, energy+enr, gvec, diag, iter, nfunc, iflag)

      ipos = 0

      call mpi_pack (iter,  1, mint, ibuf, 100, ipos, mworld, ierr)
      call mpi_pack (nfunc, 1, mint, ibuf, 100, ipos, mworld, ierr)
      call mpi_pack (iflag, 1, mint, ibuf, 100, ipos, mworld, ierr)
    end if

    ! broadcast iter, nfunc, iflag

    call mpi_bcast (ibuf, 100, mpi_packed, root, mworld, ierr)

    if ( procid /= root ) then
      ipos = 0

      call mpi_unpack (ibuf, 100, ipos, iter,  1, mint, mworld, ierr)
      call mpi_unpack (ibuf, 100, ipos, nfunc, 1, mint, mworld, ierr)
      call mpi_unpack (ibuf, 100, ipos, iflag, 1, mint, mworld, ierr)
    end if


    ! decide whether to exit

    if ( iflag <= 0 .or. iter >= maxit+1 .or. nfunc >= mxfunc+1 )  exit

    ipass = ipass + 1
  end do

  ! end of iteration

  if ( iflag < 0 ) then
    nflag = 1
  else if ( iter >= maxit ) then
    nflag = 2
  else if ( nfunc >= mxfunc ) then
    nflag = 3
  end if


  ! update darr0

  darr0 = darr1


  ! deallocate darr1

  deallocate (darr1, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'darr1', 'scfhf_lbfgs')


  ! shutdown energy evaluation

  call shutdown_hf (root, procid, imethd, iwfnty, nup, ndn)


  ! shutdown lbfgs

  if ( procid == root ) then
    deallocate (zvec, gvec, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'zvec, gvec', 'scfhf_lbfgs')

    deallocate (diag, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'diag', 'scfhf_lbfgs')

    call shutdown_lbfgs
  end if


  ! shutdown scratch arrays to update orbitals

  if ( procid == root ) then
    call shutdown_update (nup, ndn, iwfnty)
  end if


  ! choose orbitals that diagonalize H^11
  !   ( only determinants set to be optimized )

  if ( procid == root ) then

    call setup_diagh11 (iwfnty, norb, nbas, nbct, nup, ndn)

    do j = ik+1, ndet
      if ( iwfnty == 1 ) then
        call diagh11_rhf (norb, nbas, nbct, nup, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)

      else if ( iwfnty == 2 ) then
        call diagh11_uhf (norb, nbas, nbct, nup, ndn, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)

      else if ( iwfnty == 3 ) then
        call diagh11_ghf (norb, nbas, nbct, nup, ndn, darr0, enarr, xdim, &
             & edim, ndet, j, xmat, hmatt, i2sv, ni2s, lsf)
      end if
    end do

    call shutdown_diagh11 (iwfnty, nup, ndn)

  end if


  ! print determinants to output file

  if ( iprint >= 1 .and. procid == root ) then
    open (unit = iout, file = outfile, status = 'old', &
        & form = 'formatted', position = 'append')

    call print_darr_ao (iout, iwfnty, norb, nbas, xmat, darr0, enarr, &
         & xdim, edim, ndet, ik+1, ndet)

    close (unit = iout)
  end if


  ! call hamilt_hf with updated orbitals

    call setup_hf (root, procid, nproc, imethd, iwfnty, norb, nbas, nbct, &
         & nup, ndn, ndspin, ndsptl, ndcplx, ndet, ngrdx, iopt, ndtopt, &
         & nstat, idetst, zdim)

    call mpi_bcast (darr0(1,ik+1), xdim*ndtopt, mdc, root, mworld, ierr)

  call hamilt_hf (root, procid, iout, outfile, 0, imethd, iwfnty, norb, &
       & nbas, nbct, nup, ndn, xmat, hmatt, i2sv, ni2s, lsf, enr, ndspin, &
       & ndsptl, ndcplx, darr0, xdim, ndet, ndtt, ovdet, hmdet, nstat, &
       & idetst, mxdtst, inovst, inhmst, mdvec, iopt, ndtopt, imult, ngrda, &
       & ngrdb, ngrdg, grda, grdb, grdg, wgta, wgtb, wgtg, nop, rotmat, &
       & wgtmat, energy, garr, zdim, thrsh, iflg)

    call shutdown_hf (root, procid, imethd, iwfnty, nup, ndn)


  ! re-diagonalize Hamiltonian among states

  if ( procid == root ) then
    call final_diag (ndspin, ndsptl, ndcplx, ndtt, nstat, idetst, mxdtst, &
         & ovdet, hmdet, ovstat, hmstat, stvec, stval, mdvec, energy, &
         & thrsh, iflg)

    ! add nuclear repulsion energy

    stval(1:nstat) = stval(1:nstat) + enr

    if ( iflg /= 0 ) then
      write (6, *) 'error: Final diagonalization of states failed.'
      stop
    end if

    open (unit = iout, file = outfile, status = 'old', &
        & form = 'formatted', position = 'append')

    write (iout, *)
    write (iout, *) ' Final diagonalization of Hamiltonian matrix ', &
                  & '(hmstat) among states: '
    write (iout, *)

    ! print Hamiltonian and overlap among states

    call print_zmat (iout, 1, nstat, nstat, 'hmstat matrix', hmstat)
    call print_zmat (iout, 1, nstat, nstat, 'ovstat matrix', ovstat)

    ! print eigenvectors and eigenvalues

    call print_dmat (iout, 1, nstat, 1, 'hmstat eigenvalues', stval)
    call print_zmat (iout, 1, nstat, nstat, 'hmstat eigenvectors', stvec)

    close (unit = iout)
  end if


  ! shutdown spin projection grid

  deallocate (grda, grdb, grdg, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'grd?', 'scfhf_lbfgs')

  deallocate (wgta, wgtb, wgtg, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'wgt?', 'scfhf_lbfgs')


  ! shutdown spatial projection grid

  deallocate (wgtmat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'wgtmat', 'scfhf_lbfgs')

  deallocate (rotmat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rotmat', 'scfhf_lbfgs')


  ! deallocate space for zarr, larr, marr, garr

  if ( procid == root ) then
    deallocate (zarr, larr, marr, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, '?arr', 'scfhf_lbfgs')
  end if

  deallocate (garr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'garr', 'scfhf_lbfgs')


  ! deallocate space for inovst, inhmst

  if ( procid == root ) then
    if ( allocated (inovst) ) then
      deallocate (inovst, stat=istatus)
      if ( istatus /= 0 ) call error_alloc (2, 'inovst', 'scfhf_lbfgs')
    end if

    if ( allocated (inhmst) ) then
      deallocate (inhmst, stat=istatus)
      if ( istatus /= 0 ) call error_alloc (2, 'inhmst', 'scfhf_lbfgs')
    end if
  end if


  ! deallocate core Hamiltonian

  deallocate (hmatt, scrt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'hmatt', 'scfhf_lbfgs')


  return
end subroutine scfhf_lbfgs



subroutine invert_ovst (nstat, ovstat, hmstat, inovst, inhmst, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | invert_ovst  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Given a wavefunction with nstat states, compute the inverse    |
! | of the overlap matrix (ovstat) among the n-1 previous states.  |
! | This needs to be done in excited VAMP or excited FED/RES VAMP  |
! | calculations as the wavefunction for the last state needs to   |
! | be orthogonalized wrt all previous states.                     |
! |                                                                |
! | This subroutine produces two output arrays:                    |
! |                                                                |
! |   inovst - S^-1  among nstat-1 previous states                 |
! |   inhmst - S^-1 . H . S^-1 for nstat-1 previous states         |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  abort!, problems inverting ovstat found         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nstat  - number of states in wavefunction
  !   ovstat - overlap matrix among states
  !   hmstat - hamiltonian matrix among states
  !   inovst - inverse overlap among nstat-1 previous states
  !   inhmst - S-1 . H . S-1 for nstat-1 previous states
  !   iflg   - return code

  integer, intent(out) :: iflg
  integer, intent(in) :: nstat

  complex(kind=dp), dimension(nstat,nstat), intent(in) :: ovstat, hmstat
  complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: inovst, inhmst


  ! other variables

  integer :: nst1, istatus
  integer :: nb, lwrkl, info

  integer, dimension(:), allocatable :: ipiv

  complex(kind=dp), dimension(:), allocatable :: wrkl, scr


  ! functions

  integer, external :: ilaenv


  iflg = 0

  ! nothing to do if nstat = 1

  if ( nstat == 1 ) then
    return
  end if

  nst1 = nstat-1

  allocate (inovst(nst1,nst1), inhmst(nst1,nst1), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'inovst, inhmst', 'invert_ovst')


  ! trivial case if nstat = 2

  if ( nst1 == 1 ) then
    inovst(1,1) = d1 / ovstat(1,1)
    inhmst(1,1) = inovst(1,1) * inovst(1,1) * hmstat(1,1)

    return
  end if

  ! copy overlap and Hamiltonian matrices into inovst, inhmst

  inovst(1:nst1,1:nst1) = ovstat(1:nst1,1:nst1)
  inhmst(1:nst1,1:nst1) = hmstat(1:nst1,1:nst1)


  ! allocate memory for stratch arrays

  ! nb = 100
  nb = ilaenv (1, 'zgetri', ' ', nst1, -1, -1, -1)
  lwrkl = nst1*nb

  allocate (ipiv(nst1), wrkl(lwrkl), scr(nst1*nst1), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scratch arrays', 'invert_ovst')


  ! LU factorization of inovst
  ! invert inovst

  call zgetrf (nst1, nst1, inovst, nst1, ipiv, info)
  call zgetri (nst1, inovst, nst1, ipiv, wrkl, lwrkl, info)

  if ( info /= 0 ) then
    iflg = 1
    return
  end if

  ! compute
  !   inhmst = inovst . inhmst . inovst

  call zgemm ('n', 'n', nst1, nst1, nst1, z1, inovst, nst1, &
       & inhmst, nst1, z0, scr, nst1)

  call zgemm ('n', 'n', nst1, nst1, nst1, z1, scr, nst1, &
       & inovst, nst1, z0, inhmst, nst1)


  ! deallocate scratch arrays

  deallocate (ipiv, wrkl, scr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch arrays', 'invert_ovst')


  return
end subroutine invert_ovst



subroutine final_diag (ndspin, ndsptl, ndcplx, ndtt, nstat, idetst, &
     & mxdtst, ovdet, hmdet, ovstat, hmstat, stvec, stval, mdvec, &
     & energy, thrsh, iflg)

  use linalg

! +----------------------------------------------------------------+
! |                                                                |
! | final_diag  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | Perform a final diagonalization among states in an excited     |
! | VAMP or excited FED/RES VAMP calculation.                      |
! |                                                                |
! | The overlap and Hamiltonian matrix elements wrt the last       |
! | optimized state are evaluated and a diagonalization is         |
! | performed to obtain Hamiltonian eigenvectors and eigenvalues.  |
! |                                                                |
! | We assume that the rest of the overlap and Hamiltonian matrix  |
! | elements (among the nstat-1 previous states) are already       |
! | stored and shall remain unchanged.                             |
! |                                                                |
! | On output, the following arrays have been updated / modified:  |
! |                                                                |
! |   ovstat - overlap matrix among states                         |
! |   hmstat - Hamiltonian matrix among states                     |
! |   stvec  - eigenvectors of hmstat                              |
! |   stval  - eigenvalues of hmstat                               |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   ndspin - number of configs per bs determinant (spin)
  !   ndsptl - number of configs per bs determinant (spatial)
  !   ndcplx - number of configs per bs determinant (complex)
  !   ndtt   - dimension of ovdet, hmdet
  !            (ndtt = ndet*ndspin*ndsptl*ndcplx)
  !   nstat  - number of states in wavefunction
  !   idetst - array with number of determinants per state
  !   mxdtst - maximum number of determinants in a state
  !   ovdet  - overlap matrix among dets
  !   hmdet  - hamiltonian matrix among dets
  !   ovstat - overlap matrix among states [updated]
  !   hmstat - hamiltonian matrix among states [updated]
  !   stvec  - generalized eigenvectors of hmstat [updated]
  !   stval  - generalized eigenvalues of hmstat [updated]
  !   mdvec  - CI expansion of states among its dets [updated]
  !   energy - energy eigenvalue of optimized state
  !   thrsh  - threshold to use in determining linear independency
  !   iflg   - return code from zgenevs

  integer, intent(out) :: iflg
  integer, intent(in) :: ndspin, ndsptl, ndcplx, ndtt
  integer, intent(in) :: nstat, mxdtst
  integer, dimension(nstat), intent(in) :: idetst

  real(kind=dp), intent(in) :: energy, thrsh
  real(kind=dp), dimension(nstat), intent(inout) :: stval

  complex(kind=dp), dimension(nstat,nstat), intent(inout) :: stvec
  complex(kind=dp), dimension(ndtt,ndtt), intent(in) :: ovdet, hmdet
  complex(kind=dp), dimension(nstat,nstat), intent(inout) :: ovstat, hmstat
  complex(kind=dp), dimension(ndspin*ndsptl*ndcplx*mxdtst,nstat), &
       & intent(in) :: mdvec


  ! other variables

  integer :: j, ik, p1a, ndt2
  integer :: nstsq, nb, lwork, lldim
  integer :: iz1, izx1
  integer :: ist, idet, ij, ix
  integer :: istatus

  real(kind=dp), dimension(:), allocatable :: rwork, ovval

  complex(kind=dp), dimension(:), allocatable :: scr1, scr2, scr3, work


  ! functions

  integer, external :: ilaenv


  ! solve the trivial case of a single state

  iflg = 0

  if ( nstat == 1 ) then
    hmstat(1,1) = cmplx (energy, d0, dp)
    ovstat(1,1) = z1

    stval(1) = energy
    stvec(1,1) = z1
    return
  end if


  ! set some important variables

  ndt2 = ndspin*ndsptl*ndcplx*idetst(nstat)

    ik = 0
  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do
    p1a = ik+1;


  ovstat(nstat,1:nstat) = z0
  ovstat(1:nstat,nstat) = z0
  hmstat(nstat,1:nstat) = z0
  hmstat(1:nstat,nstat) = z0


  ! loop over the final row

    ij = 0
  do ist = 1, nstat
    idet = idetst(ist)
    do ix = 1, idet*ndspin*ndsptl*ndcplx
      ij = ij + 1

        izx1 = (p1a-1)*ndspin*ndsptl*ndcplx
      do iz1 = 1, ndt2
        izx1 = izx1 + 1

        ovstat(nstat,ist) = ovstat(nstat,ist) &
             & + conjg(mdvec(iz1,nstat)) * mdvec(ix,ist) * &
               & ovdet(izx1,ij)

        hmstat(nstat,ist) = hmstat(nstat,ist) &
             & + conjg(mdvec(iz1,nstat)) * mdvec(ix,ist) * &
               & hmdet(izx1,ij)
      end do
    end do
  end do

  ovstat(nstat,nstat) = real (ovstat(nstat,nstat))
  hmstat(nstat,nstat) = real (hmstat(nstat,nstat))

  do j = 1, nstat-1
    ovstat(j,nstat) = conjg (ovstat(nstat,j))
    hmstat(j,nstat) = conjg (hmstat(nstat,j))
  end do


  ! prepare for diagonalization

  nstsq = nstat*nstat

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', nstat, -1, -1, -1)
  lwork = max (1, (nb+1)*nstat)


  allocate (scr1(nstsq), scr2(nstsq), scr3(nstsq), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'final_diag')

  ! allocate (work(lwork), rwork(3*nstat-2), stat=istatus)
  allocate (work(lwork), rwork(7*nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'work', 'final_diag')

  allocate (ovval(nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'ovval', 'final_diag')


  ! perform diagonalization

  call zgenevs ('u', nstat, hmstat, ovstat, lldim, stvec, stval, &
       & ovval, scr1, scr2, scr3, work, lwork, rwork, thrsh, iflg)


  ! deallocate scratch arrays

  deallocate (ovval, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ovval', 'final_diag')

  deallocate (work, rwork, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'work', 'final_diag')

  deallocate (scr1, scr2, scr3, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scr?', 'final_diag')


  return
end subroutine final_diag



subroutine build_global_grad (norb, nup, ndn, iwfnty, garr, zdim, &
     & ndtopt, larr, marr, ldim, mdim)

! +----------------------------------------------------------------+
! |                                                                |
! | build_global_grad  --  CAJH, 11.2012                           |
! |                                                                |
! |                                                                |
! | Given the array garr with the local gradient of the energy     |
! | wrt the updated set of ndtopt determinants to optimize,        |
! | transform it to the global gradient of the energy wrt the      |
! | reference set of ndtopt determinants.                          |
! |                                                                |
! | The following transformation is performed for each det:        |
! |                                                                |
! |   G_gl  =  (M^T)^-1  . G_loc .  (L*)^-1,                       |
! |                                                                |
! | where M and L are lower triangular matrices determined from    |
! | Cholesky decompositions in update_orb*.                        |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   iwfnty - type of wavefunction to use
  !   zdim   - leading dimension of array garr
  !   ndtopt - number of determinants to update
  !   larr   - [L*^-1] array used to update determinants
  !   marr   - [M*^-1] array used to update determinants
  !   ldim   - leading dimension of array larr
  !   mdim   - leading dimension of array marr
  !   garr   - on input, local gradient for each determinant;
  !            on output, global gradient for each determinant

  integer, intent(in) :: norb, nup, ndn, iwfnty
  integer, intent(in) :: zdim, ldim, mdim, ndtopt

  complex(kind=dp), dimension(mdim,ndtopt), intent(in) :: marr
  complex(kind=dp), dimension(ldim,ndtopt), intent(in) :: larr
  complex(kind=dp), dimension(zdim,ndtopt), intent(inout) :: garr


  ! other variables

  integer :: k, nh, np
  integer :: nh_up, nh_dn, np_up, np_dn
  integer :: ldim_up, mdim_up, zdim_up


  ! some useful quantities

  nh = 0;     np = 0
  nh_up = 0;  np_up = 0
  nh_dn = 0;  np_dn = 0

  ldim_up = 0
  mdim_up = 0
  zdim_up = 0


  if ( iwfnty == 1 ) then
    nh = nup
    np = norb-nup

  else if ( iwfnty == 2 ) then
    nh_up = nup
    np_up = norb-nup
    nh_dn = ndn
    np_dn = norb-ndn

    ldim_up = nh_up*nh_up
    mdim_up = np_up*np_up
    zdim_up = np_up*nh_up

  else if ( iwfnty == 3 ) then
    nh = nup + ndn
    np = 2*norb - nup - ndn
  end if


  ! loop over determinants to update

  do k = 1, ndtopt

    ! compute  (M^T)^-1 . Z . L*^-1

    if ( iwfnty == 1 .or. iwfnty == 3 ) then

      call ztrmm ('r', 'l', 'n', 'n', np, nh, z1, larr(1,k), &
           & nh, garr(1,k), np)
      call ztrmm ('l', 'l', 'c', 'n', np, nh, z1, marr(1,k), &
           & np, garr(1,k), np)

    else if ( iwfnty == 2 ) then

      if ( nup > 0 ) then
      call ztrmm ('r', 'l', 'n', 'n', np_up, nh_up, z1, larr(1,k), &
           & nh_up, garr(1,k), np_up)
      call ztrmm ('l', 'l', 'c', 'n', np_up, nh_up, z1, marr(1,k), &
           & np_up, garr(1,k), np_up)
      end if

      if ( ndn > 0 ) then
      call ztrmm ('r', 'l', 'n', 'n', np_dn, nh_dn, z1, larr(ldim_up+1,k), &
           & nh_dn, garr(zdim_up+1,k), np_dn)
      call ztrmm ('l', 'l', 'c', 'n', np_dn, nh_dn, z1, marr(mdim_up+1,k), &
           & np_dn, garr(zdim_up+1,k), np_dn)
      end if

    end if

  end do


  return
end subroutine build_global_grad



subroutine setup_update (norb, nup, ndn, iwfnty)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_update  --  CAJH, 11.2012                                |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in update_orb.     |
! |                                                                |
! | This subroutine should be called before any calls to           |
! | update_orb_ghf, update_orb_uhf, or update_orb_rhf.             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons

  integer, intent(in) :: iwfnty, norb, nup, ndn


  ! other variables

  integer :: istatus, ntot, nosq
  integer :: ist1, ist2


  ist1 = 0
  ist2 = 0
  istatus = 0

  ntot = nup + ndn
  nosq = norb*norb


  ! allocate space for dscr

  if ( iwfnty == 1 ) then
    allocate (dscr(nosq), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (dscr(2*nosq), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (dscr(4*nosq), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'dscr', 'setup_update')

  ! allocate space for zscr

  if ( iwfnty == 1 ) then
    allocate (zscr((norb-nup)*nup), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (zscr((2*norb-ntot)*ntot), stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    allocate (zscr_up((norb-nup)*nup), stat=ist1)
    end if

    if ( ndn > 0 ) then
    allocate (zscr_dn((norb-ndn)*ndn), stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'zscr', 'setup_update')

  ! allocate space for lmat, mmat

  if ( iwfnty == 1 ) then
    allocate (lmat(nup,nup), mmat(norb-nup,norb-nup), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (lmat(ntot,ntot), mmat(2*norb-ntot,2*norb-ntot), stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    allocate (lmat_up(nup,nup), mmat_up(norb-nup,norb-nup), stat=ist1)
    end if

    if ( ndn > 0 ) then
    allocate (lmat_dn(ndn,ndn), mmat_dn(norb-ndn,norb-ndn), stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'lmat, mmat', 'setup_update')

  ! allocate space for scrh, scrp

  if ( iwfnty == 1 ) then
    allocate (scrh(norb*nup), scrp(norb*(norb-nup)), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scrh(2*norb*ntot), scrp(2*norb*(2*norb-ntot)), stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    allocate (scrh_up(norb*nup), scrp_up(norb*(norb-nup)), stat=ist1)
    end if

    if ( ndn > 0 ) then
    allocate (scrh_dn(norb*ndn), scrp_dn(norb*(norb-ndn)), stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scrh, scrp', 'setup_update')


  return
end subroutine setup_update



subroutine shutdown_update (nup, ndn, iwfnty)

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_update  --  CAJH, 11.2012                             |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in update_orb.   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons

  integer, intent(in) :: iwfnty, nup, ndn


  ! other variables

  integer :: istatus, ist1, ist2


  ist1 = 0
  ist2 = 0
  istatus = 0


  ! deallocate space for dscr

  deallocate (dscr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dscr', 'shutdown_update')

  ! deallocate space for zscr

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (zscr, stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    deallocate (zscr_up, stat=ist1)
    end if

    if ( ndn > 0 ) then
    deallocate (zscr_dn, stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (2, 'zscr', 'shutdown_update')

  ! deallocate space for lmat, mmat

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (lmat, mmat, stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    deallocate (lmat_up, mmat_up, stat=ist1)
    end if

    if ( ndn > 0 ) then
    deallocate (lmat_dn, mmat_dn, stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (2, 'lmat, mmat', 'shutdown_update')

  ! deallocate space for scrh, scrp

  if ( iwfnty == 1 .or. iwfnty == 3 ) then
    deallocate (scrh, scrp, stat=istatus)
  end if

  if ( iwfnty == 2 ) then
    if ( nup > 0 ) then
    deallocate (scrh_up, scrp_up, stat=ist1)
    end if

    if ( ndn > 0 ) then
    deallocate (scrh_dn, scrp_dn, stat=ist2)
    end if

    istatus = abs(ist1) + abs(ist2)
  end if

  if ( istatus /= 0 ) call error_alloc (2, 'scrh, scrp', 'shutdown_update')


  return
end subroutine shutdown_update



subroutine update_orb_ghf (norb, nup, ndn, darr0, darr1, xdim, ndet, &
     & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
     & mdim, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | update_orb_ghf  --  CAJH, 11.2012                              |
! |                                                                |
! |                                                                |
! | Given a Thouless matrix Z and a reference Slater determinant   |
! | |Phi>, construct an updated Slater determinant |Phi'> as       |
! |                                                                |
! |   |Phi'>  = N * exp [sum_{ph} Z_{ph} b!_p b_h] |Phi>,          |
! |                                                                |
! | where b!_p and b_h are particle-creation and hole-annihilation |
! | operators corresponding to state |Phi> (characterized by the   |
! | matrix of orbital coefficients D). N is a normalization        |
! | parameters.                                                    |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  abort!, problems with Cholesky decomposition    |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The operators b'!_p and b'_h (corresponding to |Phi'>) are     |
! | related to the operators b!_p and b_h by                       |
! |                                                                |
! |   b'!_h  =  sum_h'  [L^-1]_{hh'} *                             |
! |                     ( b!_h' + sum_p'  Z_{p'h'} b_p'! ),        |
! |                                                                |
! |   b'!_p  =  sum_p'  [M^-1]_{pp'} *                             |
! |                     ( b!_p' - sum_h' Z*_{p'h'} b_h'! ).        |
! |                                                                |
! | Here, L and M are triangular matrices obtained from Cholesky   |
! | decompositions:                                                |
! |                                                                |
! |   1  +  Z^T . Z*  =  L . L!,                                   |
! |   1  +  Z* . Z^T  =  M . M!.                                   |
! |                                                                |
! | We now use the fact that the matrix D relates the basis states |
! | to the HF hole and particle states as                          |
! |                                                                |
! |   b!_k  =  sum_i  D*_{ik}  c!_i.                               |
! |                                                                |
! | It is then straightforward to show that the matrix D'          |
! | characterizing the state |Phi'> is obtained from D by          |
! |                                                                |
! |   D'_{ih}  =  sum_h' [ D_{ih'} + sum_p' D_{ip'} Z*_{p'h'} ] *  |
! |                      ( inv(L) )*_{hh'},                        |
! |                                                                |
! |   D'_{ip}  =  sum_p' [ D_{ip'} - sum_h' D_{ih'}  Z_{h'p'} ] *  |
! |                      ( inv(M) )*_{pp'}.                        |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   darr0  - array of reference determinants
  !   xdim   - leading dimension of array darr0
  !   ndet   - total number of determinants stored

  integer, intent(in) :: norb, nup, ndn
  integer, intent(in) :: xdim, ndet

  complex(kind=dp), dimension(xdim,ndet), intent(in) :: darr0

  !   nstat  - number of states in wavefunction
  !   idetst - array with number of determinants per state

  integer, intent(in) :: nstat
  integer, dimension(nstat), intent(in) :: idetst

  !   iopt   - type of optimization (1 - FED, 2 - RES)
  !   zarr   - array of Thouless matrices to update determinants
  !   zdim   - leading dimension of array zarr
  !   ndtopt - number of determinants to update
  !   ldim   - leading dimension of array larr
  !   mdim   - leading dimension of array marr

  integer, intent(in) :: iopt, zdim, ldim, mdim, ndtopt

  complex(kind=dp), dimension(zdim,ndtopt), intent(in) :: zarr


  ! output variables

  !   darr1  - array of updated determinants
  !   larr   - [L*^-1] array used to update determinants
  !   marr   - [M*^-1] array used to update determinants
  !   iflg   - return code

  integer, intent(out) :: iflg

  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr1
  complex(kind=dp), dimension(ldim,ndtopt), intent(out) :: larr
  complex(kind=dp), dimension(mdim,ndtopt), intent(out) :: marr


  ! other variables

  integer :: nh, np, j, k, ik
  integer :: ierr, nosq
  integer :: infol, infom


  iflg = 0

  nh = nup + ndn
  np = 2*norb-nh

  nosq = norb*norb


  ! prepare initial index ik of determinants to update

  ik = 0

  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do

  if ( iopt == 1 ) then          ! for FED optimizations, only the last
    ik = ik + idetst(nstat)-1    ! determinant is updated
  end if


  ! loop over determinants to update

  ierr = 0

  do k = 1, ndtopt

    ik = ik + 1

    ! obtain zarr, darr

    zscr(1:zdim)   = conjg (zarr(1:zdim,k))
    dscr(1:4*nosq) = darr0(1:4*nosq,ik)


    ! build the matrices

    !   Lmat  =  1_{nh}  +  Z^T . Z*
    !   Mmat  =  1_{np}  +  Z* . Z^T

    call zgemm ('c', 'n', nh, nh, np, z1, zscr, np, zscr, np, &
         & z0, lmat, nh)
    call zgemm ('n', 'c', np, np, nh, z1, zscr, np, zscr, np, &
         & z0, mmat, np)

    do j = 1, nh
      lmat(j,j) = real(lmat(j,j)) + z1
    end do

    do j = 1, np
      mmat(j,j) = real(mmat(j,j)) + z1
    end do


    ! perform Cholesky decompositions to find the triangular matrices
    ! L and M:

    !   Lmat  =  L . L!
    !   Mmat  =  M . M!

    call zpotrf ('l', nh, lmat, nh, infol)
    call zpotrf ('l', np, mmat, np, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if


    ! invert L and M

    call ztrtri ('l', 'n', nh, lmat, nh, infol)
    call ztrtri ('l', 'n', np, mmat, np, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if


    ! build

    !   scrh  =  D(:,1:nh) + D(:,nh+1:nh+np) . Z*
    !   scrp  =  D(:,nh+1:nh+np) - D(:,1:nh) . Z^T

    scrh(1:2*norb*nh) = dscr(1:2*norb*nh)
    scrp(1:2*norb*np) = dscr(2*norb*nh+1:4*nosq)

    call zgemm ('n', 'n', 2*norb, nh, np, z1, dscr(2*norb*nh+1), &
         & 2*norb, zscr, np, z1, scrh, 2*norb)
    call zgemm ('n', 'c', 2*norb, np, nh, -z1, dscr, 2*norb, &
         & zscr, np, z1, scrp, 2*norb)


    ! compute

    !   scrh  =  scrh . (inv(L))!
    !   scrp  =  scrp . (inv(M))!

    call ztrmm ('r', 'l', 'c', 'n', 2*norb, nh, z1, lmat, nh, &
         & scrh, 2*norb)
    call ztrmm ('r', 'l', 'c', 'n', 2*norb, np, z1, mmat, np, &
         & scrp, 2*norb)


    ! update darr1 with new orbitals

    darr1(1:2*norb*nh,ik)        = scrh(1:2*norb*nh)
    darr1(2*norb*nh+1:4*nosq,ik) = scrp(1:2*norb*np)

    ! save L*^-1, M*^-1

    larr(1:nh*nh,k) = reshape (conjg(lmat(1:nh,1:nh)), (/ nh*nh /))
    marr(1:np*np,k) = reshape (conjg(mmat(1:np,1:np)), (/ np*np /))

  end do

  if ( ierr == 1 ) then
    iflg = 1
  end if


  return
end subroutine update_orb_ghf



subroutine update_orb_uhf (norb, nup, ndn, darr0, darr1, xdim, ndet, &
     & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
     & mdim, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | update_orb_uhf  --  CAJH, 11.2012                              |
! |                                                                |
! |                                                                |
! | ( UHF version of update_orb_ghf. )                             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   darr0  - array of reference determinants
  !   xdim   - leading dimension of array darr0
  !   ndet   - total number of determinants stored

  integer, intent(in) :: norb, nup, ndn
  integer, intent(in) :: xdim, ndet

  complex(kind=dp), dimension(xdim,ndet), intent(in) :: darr0

  !   nstat  - number of states in wavefunction
  !   idetst - array with number of determinants per state

  integer, intent(in) :: nstat
  integer, dimension(nstat), intent(in) :: idetst

  !   iopt   - type of optimization (1 - FED, 2 - RES)
  !   zarr   - array of Thouless matrices to update determinants
  !   zdim   - leading dimension of array zarr
  !   ndtopt - number of determinants to update
  !   ldim   - leading dimension of array larr
  !   mdim   - leading dimension of array marr

  integer, intent(in) :: iopt, zdim, ldim, mdim, ndtopt

  complex(kind=dp), dimension(zdim,ndtopt), intent(in) :: zarr


  ! output variables

  !   darr1  - array of updated determinants
  !   larr   - [L*^-1] array used to update determinants
  !   marr   - [M*^-1] array used to update determinants
  !   iflg   - return code

  integer, intent(out) :: iflg

  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr1
  complex(kind=dp), dimension(ldim,ndtopt), intent(out) :: larr
  complex(kind=dp), dimension(mdim,ndtopt), intent(out) :: marr


  ! other variables

  integer :: nh_up, np_up, nh_dn, np_dn, j, k, ik
  integer :: zdim_up, zdim_dn
  integer :: ierr, nosq
  integer :: infol, infom


  iflg = 0

  nh_up = nup
  np_up = norb-nh_up

  nh_dn = ndn
  np_dn = norb-nh_dn

  zdim_up = np_up * nh_up
  zdim_dn = np_dn * nh_dn

  nosq = norb*norb


  ! prepare initial index ik of determinants to update

  ik = 0

  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do

  if ( iopt == 1 ) then          ! for FED optimizations, only the last
    ik = ik + idetst(nstat)-1    ! determinant is updated
  end if


  ! loop over determinants to update

  ierr = 0

  do k = 1, ndtopt

    ik = ik + 1

    ! obtain zarr, darr

    if ( nup > 0 ) then
    zscr_up(1:zdim_up) = conjg (zarr(1:zdim_up,k))
    end if

    if ( ndn > 0 ) then
    zscr_dn(1:zdim_dn) = conjg (zarr(zdim_up+1:zdim_up+zdim_dn,k))
    end if

    dscr(1:2*nosq) = darr0(1:2*nosq,ik)


    ! build the matrices

    !   Lmat  =  1_{nh}  +  Z^T . Z*
    !   Mmat  =  1_{np}  +  Z* . Z^T

    if ( nup > 0 ) then
    call zgemm ('c', 'n', nh_up, nh_up, np_up, z1, zscr_up, np_up, &
         & zscr_up, np_up, z0, lmat_up, nh_up)
    call zgemm ('n', 'c', np_up, np_up, nh_up, z1, zscr_up, np_up, &
         & zscr_up, np_up, z0, mmat_up, np_up)

    do j = 1, nh_up
      lmat_up(j,j) = real(lmat_up(j,j)) + z1
    end do

    do j = 1, np_up
      mmat_up(j,j) = real(mmat_up(j,j)) + z1
    end do
    end if

    if ( ndn > 0 ) then
    call zgemm ('c', 'n', nh_dn, nh_dn, np_dn, z1, zscr_dn, np_dn, &
         & zscr_dn, np_dn, z0, lmat_dn, nh_dn)
    call zgemm ('n', 'c', np_dn, np_dn, nh_dn, z1, zscr_dn, np_dn, &
         & zscr_dn, np_dn, z0, mmat_dn, np_dn)

    do j = 1, nh_dn
      lmat_dn(j,j) = real(lmat_dn(j,j)) + z1
    end do

    do j = 1, np_dn
      mmat_dn(j,j) = real(mmat_dn(j,j)) + z1
    end do
    end if


    ! perform Cholesky decompositions to find the triangular matrices
    ! L and M:

    !   Lmat  =  L . L!
    !   Mmat  =  M . M!

    if ( nup > 0 ) then
    call zpotrf ('l', nh_up, lmat_up, nh_up, infol)
    call zpotrf ('l', np_up, mmat_up, np_up, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if
    end if

    if ( ndn > 0 ) then
    call zpotrf ('l', nh_dn, lmat_dn, nh_dn, infol)
    call zpotrf ('l', np_dn, mmat_dn, np_dn, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if
    end if


    ! invert L and M

    if ( nup > 0 ) then
    call ztrtri ('l', 'n', nh_up, lmat_up, nh_up, infol)
    call ztrtri ('l', 'n', np_up, mmat_up, np_up, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if
    end if

    if ( ndn > 0 ) then
    call ztrtri ('l', 'n', nh_dn, lmat_dn, nh_dn, infol)
    call ztrtri ('l', 'n', np_dn, mmat_dn, np_dn, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if
    end if


    ! build

    !   scrh  =  D(:,1:nh) + D(:,nh+1:nh+np) . Z*
    !   scrp  =  D(:,nh+1:nh+np) - D(:,1:nh) . Z^T

    if ( nup > 0 ) then
    scrh_up(1:norb*nh_up) = dscr(1:norb*nh_up)
    scrp_up(1:norb*np_up) = dscr(norb*nh_up+1:nosq)

    call zgemm ('n', 'n', norb, nh_up, np_up, z1, dscr(norb*nh_up+1), &
         & norb, zscr_up, np_up, z1, scrh_up, norb)
    call zgemm ('n', 'c', norb, np_up, nh_up, -z1, dscr, norb, &
         & zscr_up, np_up, z1, scrp_up, norb)
    end if

    if ( ndn > 0 ) then
    scrh_dn(1:norb*nh_dn) = dscr(nosq+1:nosq+norb*nh_dn)
    scrp_dn(1:norb*np_dn) = dscr(nosq+norb*nh_dn+1:2*nosq)

    call zgemm ('n', 'n', norb, nh_dn, np_dn, z1, dscr(nosq+norb*nh_dn+1), &
         & norb, zscr_dn, np_dn, z1, scrh_dn, norb)
    call zgemm ('n', 'c', norb, np_dn, nh_dn, -z1, dscr(nosq+1), norb, &
         & zscr_dn, np_dn, z1, scrp_dn, norb)
    end if


    ! compute

    !   scrh  =  scrh . (inv(L))!
    !   scrp  =  scrp . (inv(M))!

    if ( nup > 0 ) then
    call ztrmm ('r', 'l', 'c', 'n', norb, nh_up, z1, lmat_up, &
         & nh_up, scrh_up, norb)
    call ztrmm ('r', 'l', 'c', 'n', norb, np_up, z1, mmat_up, &
         & np_up, scrp_up, norb)
    end if

    if ( ndn > 0 ) then
    call ztrmm ('r', 'l', 'c', 'n', norb, nh_dn, z1, lmat_dn, &
         & nh_dn, scrh_dn, norb)
    call ztrmm ('r', 'l', 'c', 'n', norb, np_dn, z1, mmat_dn, &
         & np_dn, scrp_dn, norb)
    end if


    ! update darr1 with new orbitals

    if ( nup > 0 ) then
    darr1(1:norb*nh_up,ik)      = scrh_up(1:norb*nh_up)
    darr1(norb*nh_up+1:nosq,ik) = scrp_up(1:norb*np_up)
    end if

    if ( ndn > 0 ) then
    darr1(nosq+1:nosq+norb*nh_dn,ik)   = scrh_dn(1:norb*nh_dn)
    darr1(nosq+norb*nh_dn+1:2*nosq,ik) = scrp_dn(1:norb*np_dn)
    end if

    ! save L*^-1, M*^-1

    if ( nup > 0 ) then
    larr(1:nh_up*nh_up,k) = &
         & reshape (conjg(lmat_up(1:nh_up,1:nh_up)), (/ nh_up*nh_up /))
    marr(1:np_up*np_up,k) = &
         & reshape (conjg(mmat_up(1:np_up,1:np_up)), (/ np_up*np_up /))
    end if

    if ( ndn > 0 ) then
    larr(nh_up*nh_up+1:nh_up*nh_up+nh_dn*nh_dn,k) = &
         & reshape (conjg(lmat_dn(1:nh_dn,1:nh_dn)), (/ nh_dn*nh_dn /))
    marr(np_up*np_up+1:np_up*np_up+np_dn*np_dn,k) = &
         & reshape (conjg(mmat_dn(1:np_dn,1:np_dn)), (/ np_dn*np_dn /))
    end if

  end do

  if ( ierr == 1 ) then
    iflg = 1
  end if


  return
end subroutine update_orb_uhf



subroutine update_orb_rhf (norb, nup, darr0, darr1, xdim, ndet, &
     & nstat, idetst, iopt, zarr, zdim, ndtopt, larr, marr, ldim, &
     & mdim, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | update_orb_rhf  --  CAJH, 11.2012                              |
! |                                                                |
! |                                                                |
! | ( RHF version of update_orb_ghf. )                             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   darr0  - array of reference determinants
  !   xdim   - leading dimension of array darr0
  !   ndet   - total number of determinants stored

  integer, intent(in) :: norb, nup
  integer, intent(in) :: xdim, ndet

  complex(kind=dp), dimension(xdim,ndet), intent(in) :: darr0

  !   nstat  - number of states in wavefunction
  !   idetst - array with number of determinants per state

  integer, intent(in) :: nstat
  integer, dimension(nstat), intent(in) :: idetst

  !   iopt   - type of optimization (1 - FED, 2 - RES)
  !   zarr   - array of Thouless matrices to update determinants
  !   zdim   - leading dimension of array zarr
  !   ndtopt - number of determinants to update
  !   ldim   - leading dimension of array larr
  !   mdim   - leading dimension of array marr

  integer, intent(in) :: iopt, zdim, ldim, mdim, ndtopt

  complex(kind=dp), dimension(zdim,ndtopt), intent(in) :: zarr


  ! output variables

  !   darr1  - array of updated determinants
  !   larr   - [L*^-1] array used to update determinants
  !   marr   - [M*^-1] array used to update determinants
  !   iflg   - return code

  integer, intent(out) :: iflg

  complex(kind=dp), dimension(xdim,ndet), intent(inout) :: darr1
  complex(kind=dp), dimension(ldim,ndtopt), intent(out) :: larr
  complex(kind=dp), dimension(mdim,ndtopt), intent(out) :: marr


  ! other variables

  integer :: nh, np, j, k, ik
  integer :: ierr, nosq
  integer :: infol, infom


  iflg = 0

  nh = nup
  np = norb-nh

  nosq = norb*norb


  ! prepare initial index ik of determinants to update

  ik = 0

  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do

  if ( iopt == 1 ) then          ! for FED optimizations, only the last
    ik = ik + idetst(nstat)-1    ! determinant is updated
  end if


  ! loop over determinants to update

  ierr = 0

  do k = 1, ndtopt

    ik = ik + 1

    ! obtain zarr, darr

    zscr(1:zdim) = conjg (zarr(1:zdim,k))
    dscr(1:nosq) = darr0(1:nosq,ik)


    ! build the matrices

    !   Lmat  =  1_{nh}  +  Z^T . Z*
    !   Mmat  =  1_{np}  +  Z* . Z^T

    call zgemm ('c', 'n', nh, nh, np, z1, zscr, np, zscr, np, &
         & z0, lmat, nh)
    call zgemm ('n', 'c', np, np, nh, z1, zscr, np, zscr, np, &
         & z0, mmat, np)

    do j = 1, nh
      lmat(j,j) = real(lmat(j,j)) + z1
    end do

    do j = 1, np
      mmat(j,j) = real(mmat(j,j)) + z1
    end do


    ! perform Cholesky decompositions to find the triangular matrices
    ! L and M:

    !   Lmat  =  L . L!
    !   Mmat  =  M . M!

    call zpotrf ('l', nh, lmat, nh, infol)
    call zpotrf ('l', np, mmat, np, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if


    ! invert L and M

    call ztrtri ('l', 'n', nh, lmat, nh, infol)
    call ztrtri ('l', 'n', np, mmat, np, infom)

    if ( infol > 0 .or. infom > 0 ) then
      ierr = 1
      exit
    end if


    ! build

    !   scrh  =  D(:,1:nh) + D(:,nh+1:nh+np) . Z*
    !   scrp  =  D(:,nh+1:nh+np) - D(:,1:nh) . Z^T

    scrh(1:norb*nh) = dscr(1:norb*nh)
    scrp(1:norb*np) = dscr(norb*nh+1:nosq)

    call zgemm ('n', 'n', norb, nh, np, z1, dscr(norb*nh+1), &
         & norb, zscr, np, z1, scrh, norb)
    call zgemm ('n', 'c', norb, np, nh, -z1, dscr, norb, &
         & zscr, np, z1, scrp, norb)


    ! compute

    !   scrh  =  scrh . (inv(L))!
    !   scrp  =  scrp . (inv(M))!

    call ztrmm ('r', 'l', 'c', 'n', norb, nh, z1, lmat, nh, &
         & scrh, norb)
    call ztrmm ('r', 'l', 'c', 'n', norb, np, z1, mmat, np, &
         & scrp, norb)


    ! update darr1 with new orbitals

    darr1(1:norb*nh,ik)      = scrh(1:norb*nh)
    darr1(norb*nh+1:nosq,ik) = scrp(1:norb*np)

    ! save L*^-1, M*^-1

    larr(1:nh*nh,k) = reshape (conjg(lmat(1:nh,1:nh)), (/ nh*nh /))
    marr(1:np*np,k) = reshape (conjg(mmat(1:np,1:np)), (/ np*np /))

  end do

  if ( ierr == 1 ) then
    iflg = 1
  end if


  return
end subroutine update_orb_rhf



subroutine print_darr_ao (iout, iwfnty, norb, nbas, xmat, darr, enarr, &
     & xdim, edim, ndet, i1, i2)

! +----------------------------------------------------------------+
! |                                                                |
! | print_darr_ao  --  CAJH, 01.2013                               |
! |                                                                |
! |                                                                |
! | Print the orbital coefficients and energies of some dets in    |
! | the array darr to the output file in a nice format.            |
! |                                                                |
! | The orbital coefficients are printed in the std AO basis.      |
! |                                                                |
! | The indices i1 and i2 determine which determinants to print:   |
! |   i1 - index of 1st determinant to print                       |
! |   i2 - index of last determinant to print                      |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iout   - output file unit number
  !   iwfnty - type of wavefunction to use
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]

  integer, intent(in) :: iout, iwfnty, norb, nbas

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat

  !   darr   - array of determinants
  !   enarr  - array with orbital energies
  !   xdim   - leading dimension of array darr0
  !   edim   - leading dimension of array enarr
  !   ndet   - total number of determinants stored
  !   i1     - index of 1st determinant to print
  !   i2     - index of last determinant to print

  integer, intent(in) :: xdim, edim, ndet, i1, i2

  real(kind=dp), dimension(edim,ndet), intent(in) :: enarr
  complex(kind=dp), dimension(xdim,ndet), intent(in) :: darr


  ! other variables

  integer :: j, xdimx, nosq
  integer :: istatus

  complex(kind=dp), dimension(:), allocatable :: dscr

  complex(kind=dp), dimension(:,:), allocatable :: scr_s1, scr_s2
  complex(kind=dp), dimension(:,:), allocatable :: dmat_s1, dmat_s2

  character(len=3) :: idet
  character(len=23) :: ennam
  character(len=27) :: monam
  character(len=27) :: ennam_up, ennam_dn
  character(len=31) :: monam_up, monam_dn


  nosq = norb*norb

  xdimx = 0

  if ( iwfnty == 1 ) then
    xdimx = nbas*norb
  else if ( iwfnty == 2 ) then
    xdimx = 2*nbas*norb
  else if ( iwfnty == 3 ) then
    xdimx = 4*nbas*norb
  end if


  ! allocate space for D in std AO basis

  allocate (dscr(xdimx), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'dscr', 'print_darr_ao')

  ! scratch arrays for GHF wfns

  if ( iwfnty == 3 ) then
    allocate (dmat_s1(2*nbas,2*norb), dmat_s2(2*norb,2*norb), &
            & scr_s1(nbas,norb), scr_s2(norb,norb), stat=istatus)
  else
    allocate (dmat_s1(1,1), dmat_s2(1,1), &
            & scr_s1(1,1), scr_s2(1,1), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scratch', 'print_darr_ao')


  ! loop over determinants in darr

  write (iout, *)

  do j = i1, i2

    ! transform to std AO basis
    ! perform  D  <-  X . D

    if ( iwfnty == 1 ) then
        call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
             & darr(1,j), norb, z0, dscr, nbas)

    else if ( iwfnty == 2 ) then
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & darr(1,j), norb, z0, dscr(1), nbas)
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & darr(nosq+1,j), norb, z0, dscr(nbas*norb+1), nbas)

    else if ( iwfnty == 3 ) then
      dmat_s2(1:2*norb,1:2*norb) = &
           & reshape (darr(1:4*nosq,j), (/ 2*norb, 2*norb /))

      scr_s2(1:norb,1:norb) = dmat_s2(1:norb,1:norb)
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & scr_s2, norb, z0, scr_s1, nbas)
      dmat_s1(1:nbas,1:norb) = scr_s1(1:nbas,1:norb)

      scr_s2(1:norb,1:norb) = dmat_s2(1:norb,norb+1:2*norb)
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & scr_s2, norb, z0, scr_s1, nbas)
      dmat_s1(1:nbas,norb+1:2*norb) = scr_s1(1:nbas,1:norb)

      scr_s2(1:norb,1:norb) = dmat_s2(norb+1:2*norb,1:norb)
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & scr_s2, norb, z0, scr_s1, nbas)
      dmat_s1(nbas+1:2*nbas,1:norb) = scr_s1(1:nbas,1:norb)

      scr_s2(1:norb,1:norb) = dmat_s2(norb+1:2*norb,norb+1:2*norb)
      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & scr_s2, norb, z0, scr_s1, nbas)
      dmat_s1(nbas+1:2*nbas,norb+1:2*norb) = scr_s1(1:nbas,1:norb)

      dscr(1:4*nbas*norb) = &
           & reshape (dmat_s1(1:2*nbas,1:2*norb), (/ 4*nbas*norb /))
    end if


    ! print matrices

    write (idet, '(I3)') j

    if ( iwfnty == 1 .or. iwfnty == 3 ) then
      ennam = 'MO energies [det = ' // idet // ']'
      monam = 'MO coefficients [det = ' // idet // ']'

    else if ( iwfnty == 2 ) then
      ennam_up = 'MO energies [up, det = ' // idet // ']'
      monam_up = 'MO coefficients [up, det = ' // idet // ']'
      ennam_dn = 'MO energies [dn, det = ' // idet // ']'
      monam_dn = 'MO coefficients [dn, det = ' // idet // ']'
    end if

    if ( iwfnty == 1 ) then
      call print_dmat (iout, 1, norb, 1, ennam, enarr(1,j))
      call print_zmat (iout, 1, nbas, norb, monam, dscr)

    else if ( iwfnty == 2 ) then
      call print_dmat (iout, 1, norb, 1, ennam_up, enarr(1,j))
      call print_zmat (iout, 1, nbas, norb, monam_up, dscr(1))

      call print_dmat (iout, 1, norb, 1, ennam_dn, enarr(norb+1,j))
      call print_zmat (iout, 1, nbas, norb, monam_dn, dscr(nbas*norb+1))

    else if ( iwfnty == 3 ) then
      call print_dmat (iout, 1, 2*norb, 1, ennam, enarr(1,j))
      call print_zmat (iout, 1, 2*nbas, 2*norb, monam, dscr)
    end if
  end do

  write (iout, *)


  ! deallocate space for D in std AO basis

  deallocate (dscr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dscr', 'print_darr_ao')

  ! deallocate scratch arrays

  deallocate (dmat_s1, dmat_s2, scr_s1, scr_s2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch', 'print_darr_ao')


  return
end subroutine print_darr_ao


end module phfscf


