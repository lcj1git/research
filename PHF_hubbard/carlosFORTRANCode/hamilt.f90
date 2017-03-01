

module hamilt

  use constants
  use util
  use linalg
  use trdnmt
  use erictr

  implicit none
  private
  save

  public :: setup_hf, shutdown_hf
  public :: hamilt_hf

! +----------------------------------------------------------------+
! |                                                                |
! | hamilt                                                         |
! |                                                                |
! |                                                                |
! | A collection of subroutines related to the evaluation of the   |
! | energy of a wavefunction expressed as a linear combination of  |
! | projected Hartree-Fock states.                                 |
! |                                                                |
! | The gradient of the energy wrt the underlying deformed Slater  |
! | determinants is also computed.                                 |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The evaluation of the energy with projected HF states relies   |
! | ultimately on the computation of overlap and Hamiltonian       |
! | matrix elements among non-orthogonal Slater determinants.      |
! |                                                                |
! | Given two non-orthogonal Slater determinants |phi_1> and       |
! | |phi_2>, the overlap is given by                               |
! |                                                                |
! |   ov(1,2)  =  <phi_1 | phi_2>  =  det_N  (D1^T . D2*),         |
! |                                                                |
! | where det_N implies that the determinant is over the NxN       |
! | occupied space. Here, D1 and D2 are the matrices of orbital    |
! | coefficients characterizing the Slater determinants |phi_1>    |
! | and |phi_2>, respectively.                                     |
! |                                                                |
! | Matrix elements of arbitrary operators can be computed using   |
! | transition density matrices. The transition density matrix     |
! | rho^[12] is computed as                                        |
! |                                                                |
! |                    <phi_1 | c!_i c_k |phi_2>                   |
! |   rho^[12]_ki  =  ---------------------------                  |
! |                         <phi_1 | phi_2>                        |
! |                                                                |
! |                =  sum_h   D1_{ih} . D1*_{kh}                   |
! |                +  sum_ph  D1_{ih} . Zc_{ph} . D1*_{kp},        |
! |                                                                |
! | where                                                          |
! |                                                                |
! |   Zc_{ph}  =  sum_h'  (D1^T . D2*)_{ph'} . (Lc^-1)_{h'h},      |
! |   Lc_{h'h} =  (D1^T . D2*)_{h'h}.                              |
! |                                                                |
! | The transition density matrix rho^[12] can be contracted with  |
! | antisymmetrized two-particle integrals to yield the matrix     |
! | Gam^[12] as                                                    |
! |                                                                |
! |   Gam^[12]_{ik}  =  sum_jl  <ij | v | kl>  rho^[12]_{lj}.      |
! |                                                                |
! | The Hamiltonian matrix element can then be evaluated as        |
! |                                                                |
! |                                      <phi_1 | H | phi_2>       |
! |   hmov(1,2)  =  <phi_1 | phi_2>  .  ---------------------      |
! |                                        <phi_1 | phi_2>         |
! |                                                                |
! |              =  ov(1,2) * ( Tr [T.rho] + 1/2 Tr [Gam.rho] ),   |
! |                                                                |
! | where T is the matrix of one-particle integrals.               |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | LOCAL GRADIENT EVALUATION                                      |
! |                                                                |
! | Let us assume that our energy expression is written as         |
! |                                                                |
! |          sum_{j,k}  c_j* c_k  <phi_j | H | phi_k>              |
! |   E  =  ------------------------------------------ .           |
! |            sum_{j,k}  c_j* c_k  <phi_j | phi_k>                |
! |                                                                |
! | The local gradient wrt the i-th determinant is defined as      |
! |                                                                |
! |   G^i_{ph}  =  d E / d(Z*^i_{ph}) | [Z^i_{ph} = 0],            |
! |                                                                |
! | where Z^i is the Thouless parametrization of the i-th det.     |
! | One can show that                                              |
! |                                                                |
! |              sum_k c_i* c_k <phi_i | b!_h b_p (H - E) |phi_k>  |
! |   G^i_{ph} = ------------------------------------------------  |
! |                    sum_{jk}  c_j* c_k  <phi_j | phi_k>         |
! |                                                                |
! | where b!_h and b_p are, respectively, the annihilation of a    |
! | hole and the creation of a particle operators, defined wrt     |
! | |phi_i>. E is the energy (see equation above).                 |
! |                                                                |
! | The gradient can be evaluated as                               |
! |                                                                |
! |   G^i_{ph} =  sum_k c_i* c_k . yg(i,k) . [ hmov(i,k) - E ] .   |
! |                 [ (D^i)^T . rho^[ik] . (D^i)* ]_{ph}           |
! |            +  sum_k c_i* c_k . yg(i,k) .                       |
! |                 [ (D^i)^T . (I - rho^[ik]) .                   |
! |                     (T + Gam^[ik]) . rho^[ik] . (D^i)* ]_{ph}  | 
! |                                                                |
! | where D^i is the D-matrix of the i-th determinant and          |
! |                                                                |
! |   yg(i,k)  =  ov(i,k) / { sum_{jl} c_j* c_l  ov(j,l) }.        |
! |                                                                |
! +----------------------------------------------------------------+


  ! global variables
  ! scratch arrays

  integer :: lwrkx, jstrt, ntask
  integer, dimension(:,:), allocatable :: indx

  complex(kind=dp), dimension(:,:), allocatable :: ovlp, hmlt
  complex(kind=dp), dimension(:,:,:), allocatable :: garr1, garr2

  complex(kind=dp), dimension(:,:), allocatable :: ovlpr, hmltr
  complex(kind=dp), dimension(:,:), allocatable :: ovmt, hmmt
  complex(kind=dp), dimension(:,:,:), allocatable :: garr1r, garr2r

  complex(kind=dp), dimension(:), allocatable :: d1st, d2st
  complex(kind=dp), dimension(:,:), allocatable :: rho, gam, tpgam

  complex(kind=dp), dimension(:), allocatable :: d2stY, d2stX

  complex(kind=dp), dimension(:), allocatable :: gr1sc, gr2sc
  complex(kind=dp), dimension(:), allocatable :: gr1t, gr2t

  complex(kind=dp), dimension(:,:), allocatable :: xvec
  complex(kind=dp), dimension(:), allocatable :: scx1, scx2, scx3
  complex(kind=dp), dimension(:), allocatable :: wrkx

  real(kind=dp), dimension(:), allocatable :: xval, xoval
  real(kind=dp), dimension(:), allocatable :: rwrkx

  ! bldgrad scratch arrays

  complex(kind=dp), dimension(:,:), allocatable :: scr1, scr2

  !$omp  threadprivate(indx, ovlp, hmlt, garr1, garr2, &
  !$omp&               d1st, d2st, d2stY, d2stX, rho, gam, tpgam, &
  !$omp&               gr1sc, gr2sc, gr1t, gr2t, scr1, scr2)


contains



subroutine setup_hf (root, procid, nproc, imethd, iwfnty, norb, nbas, &
     & nbct, nup, ndn, ndspin, ndsptl, ndcplx, ndet, ngrdx, iopt, &
     & ndtopt, nstat, idetst, gdim)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_hf  --  CAJH, 11.2012                                    |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in hamilt_hf.      |
! | Some global variables are also set up.                         |
! |                                                                |
! | This subroutine should be called before hamilt_hf. Once the    |
! | program has done all calls to hamilt_hf, the shutdown routine  |
! | should be called to clear all scratch arrays.                  |
! |                                                                |
! +----------------------------------------------------------------+

  include 'mpif.h'

  ! input variables

  !   root   - id of root process
  !   procid - id of each process
  !   nproc  - total number of processes running

  integer, intent(in) :: root, procid, nproc

  !   imethd  - PHF method to use
  !   iwfnty  - type of wavefunction to use
  !   norb    - number of orbitals
  !   nbas    - number of basis functions
  !   nbct    - number of Cartesian basis functions
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons

  integer, intent(in) :: imethd, iwfnty, norb, nbas, nbct, nup, ndn

  !   ndspin - number of configs per bs determinant (spin)
  !   ndsptl - number of configs per bs determinant (spatial)
  !   ndcplx - number of configs per bs determinant (complex)
  !   ndet   - total number of determinants stored
  !   ngrdx  - total number of grid points

  integer, intent(in) :: ndspin, ndsptl, ndcplx, ndet, ngrdx

  !   iopt   - type of optimization (1 - FED, 2 - RES)
  !   ndtopt - number of determinants to update
  !   nstat  - number of states in wavefunction
  !   idetst - array with number of determinants per state

  integer, intent(in) :: iopt, ndtopt, nstat
  integer, dimension(nstat), intent(in) :: idetst

  !   gdim   - dimension of each gradient vector

  integer, intent(in) :: gdim


  ! other variables

  integer :: ik, j, k1, k2, ij, nb
  integer :: nosq, ntot, iwfnt1, gdim1, ndt1, ndt2, ndtt, nevl
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: istatus


  ! functions

  integer, external :: ilaenv


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


  ! set some important variables

  nosq = norb*norb
  ntot = nup + ndn

  iwfnt1 = iwfnty
  if ( lsuhf ) iwfnt1 = 3

  gdim1 = gdim
  if ( lsuhf ) gdim1 = (2*norb-ntot)*ntot

  ndt1 = ndspin*ndsptl*ndcplx*ndtopt
  ndt2 = ndspin*ndsptl*ndcplx*idetst(nstat)
  ndtt = ndspin*ndsptl*ndcplx*ndet
  nevl = ndet*ndtopt*ngrdx*ndcplx


  ! distribute tasks over processors

  ntask = nevl/nproc

  if ( procid < mod(nevl,nproc) ) then
    jstrt = procid*(ntask+1) + 1
    ntask = ntask+1
  else
    jstrt = procid*ntask + mod(nevl,nproc) + 1
  end if


  allocate (ovlpr(ndt1,ndtt), hmltr(ndt1,ndtt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'ovlpr, hmltr', 'setup_hf')

  allocate (garr1r(gdim,ndt1,ndtt), garr2r(gdim,ndt1,ndtt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'garr?r', 'setup_hf')

  if ( procid == root ) then
    allocate (ovmt(ndt2,ndt2), hmmt(ndt2,ndt2), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'ovmt, hmmt', 'setup_hf')
  end if


  !$omp  parallel default(private) &
  !$omp& shared(procid, root, lsptl, lspin, lsuhf, lsghf, &
  !$omp&        iwfnty, iwfnt1, nup, ndn, &
  !$omp&        nosq, ntot, ndt1, ndt2, ndtt, &
  !$omp&        ik, ndtopt, ndet, &
  !$omp&        gdim, gdim1, norb, nbas, nbct)

  ! prepare indexing array of states over which we need to loop

  allocate (indx(ndtopt*ndet,3), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'indx', 'setup_hf')

    j = 1;  ij = 1
  do k1 = ik+1, ik+ndtopt
    do k2 = 1, ndet
      indx(j,1) = k1
      indx(j,2) = k2
      indx(j,3) = ij    ! ij1, k1 represent the same det; different indexing

      j = j + 1
    end do
    ij = ij + 1
  end do


  ! scratch overlap and Hamiltonian arrays

  allocate (ovlp(ndt1,ndtt), hmlt(ndt1,ndtt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'ovlp, hmlt', 'setup_hf')


  ! scratch gradient array

  allocate (garr1(gdim,ndt1,ndtt), garr2(gdim,ndt1,ndtt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'garr?', 'setup_hf')


  ! allocate space for d1st, d2st

  if ( iwfnt1 == 1 ) then
    allocate (d1st(nosq), d2st(nosq), stat=istatus)
  else if ( iwfnt1 == 2 ) then
    allocate (d1st(2*nosq), d2st(2*nosq), stat=istatus)
  else if ( iwfnt1 == 3 ) then
    allocate (d1st(4*nosq), d2st(4*nosq), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'd?st', 'setup_hf')

  ! some scratch versions of d1st, d2st

  if ( lsptl ) then
    if ( iwfnty == 1 ) then
      allocate (d2stY(nosq), stat=istatus)
    else if ( iwfnty == 2 ) then
      allocate (d2stY(2*nosq), stat=istatus)
    else if ( iwfnty == 3 ) then
      allocate (d2stY(4*nosq), stat=istatus)
    end if

    if ( istatus /= 0 ) call error_alloc (1, 'd2stY', 'setup_hf')
  end if

  if ( lspin ) then
    allocate (d2stX(4*nosq), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'd2stX', 'setup_hf')
  end if

  ! allocate space for rho

  if ( iwfnt1 == 1 ) then
    allocate (rho(norb,norb), stat=istatus)
  else if ( iwfnt1 == 2 ) then
    allocate (rho(norb,2*norb), stat=istatus)
  else if ( iwfnt1 == 3 ) then
    allocate (rho(2*norb,2*norb), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'rho', 'setup_hf')

  ! allocate space for gam, tpgam

  if ( iwfnt1 == 1 ) then
    allocate (gam(norb,norb), tpgam(norb,norb), stat=istatus)
  else if ( iwfnt1 == 2 ) then
    allocate (gam(norb,2*norb), tpgam(norb,2*norb), stat=istatus)
  else if ( iwfnt1 == 3 ) then
    allocate (gam(2*norb,2*norb), tpgam(2*norb,2*norb), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'gam', 'setup_hf')


  ! setup trdnmt scratch

  call setup_trdnmt (iwfnt1, norb, nup, ndn)

  ! setup erictr scratch

  call setup_erictr (iwfnt1, norb, nbas, nbct)

  ! setup bldgrad scratch

  call setup_bldgrad (iwfnt1, norb, nup, ndn)


  ! allocate gr1sc, gr2sc

  allocate (gr1sc(gdim1), gr2sc(gdim1), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'gr?sc', 'setup_hf')


  ! allocate space for gr1t, gr2t

  allocate (gr1t(gdim), gr2t(gdim), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'gr?t', 'setup_hf')

  !$omp  end parallel


  ! scratch arrays to perform diagonalization

  if ( procid == root ) then
    allocate (xvec(ndt2,ndt2), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'xvec', 'setup_hf')

    allocate (xval(ndt2), xoval(ndt2), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'xval, xoval', 'setup_hf')

    allocate (scx1(ndt2*ndt2), scx2(ndt2*ndt2), &
            & scx3(ndt2*ndt2), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'scx?', 'setup_hf')

    ! nb = 2
    nb = ilaenv (1, 'zhetrd', 'u', ndt2, -1, -1, -1)
    lwrkx = max (1, (nb+1)*ndt2)

    allocate (wrkx(lwrkx), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'wrkx', 'setup_hf')

    ! allocate (rwrkx(3*ndt2-2), stat=istatus)
    allocate (rwrkx(7*ndt2), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'rwrkx', 'setup_hf')
  end if


  return
end subroutine setup_hf



subroutine shutdown_hf (root, procid, imethd, iwfnty, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_hf  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in hamilt_hf.    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   root   - id of root process
  !   procid - id of each process

  integer, intent(in) :: root, procid

  !   imethd  - PHF method to use
  !   iwfnty  - type of wavefunction to use
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons

  integer, intent(in) :: imethd, iwfnty, nup, ndn


  ! other variables

  integer :: iwfnt1
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: istatus


  lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
  lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
  lspin = mod(imethd-1,2) == 1

  lsuhf = .false.
  lsghf = .false.

  if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
  if ( lspin .and. iwfnty == 3 )  lsghf = .true.

  iwfnt1 = iwfnty
  if ( lsuhf ) iwfnt1 = 3


  deallocate (ovlpr, hmltr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ovlpr, hmltr', 'shutdown_hf')

  deallocate (garr1r, garr2r, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'garr?r', 'shutdown_hf')

  if ( procid == root ) then
    deallocate (ovmt, hmmt, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'ovmt, hmmt', 'shutdown_hf')
  end if


  !$omp  parallel default(private) &
  !$omp& shared(procid, root, lsptl, lspin, lsuhf, lsghf, &
  !$omp&        iwfnty, iwfnt1, nup, ndn)

  ! deallocate indexing array

  deallocate (indx, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'indx', 'shutdown_hf')


  ! deallocate some scratch arrays

  deallocate (ovlp, hmlt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ovlp, hmlt', 'shutdown_hf')

  deallocate (garr1, garr2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'garr?', 'shutdown_hf')


  ! deallocate some variables

  deallocate (d1st, d2st, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'd?st', 'shutdown_hf')

  if ( lsptl ) then
    deallocate (d2stY, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'd2stY', 'shutdown_hf')
  end if

  if ( lspin ) then
    deallocate (d2stX, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'd2stX', 'shutdown_hf')
  end if

  deallocate (rho, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rho', 'shutdown_hf')

  deallocate (gam, tpgam, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'gam', 'shutdown_hf')


  ! shutdown trdnmt scratch

  call shutdown_trdnmt (iwfnt1, nup, ndn)

  ! shutdown erictr scratch

  call shutdown_erictr

  ! shutdown bldgrad scratch

  call shutdown_bldgrad


  ! deallocate gr1sc, gr2sc

  deallocate (gr1sc, gr2sc, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'gr?sc', 'shutdown_hf')


  ! deallocate space for gr1t, gr2t

  deallocate (gr1t, gr2t, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'gr?t', 'shutdown_hf')

  !$omp  end parallel


  ! deallocate scratch arrays to perform diagonalization

  if ( procid == root ) then
    deallocate (xvec, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'xvec', 'shutdown_hf')

    deallocate (xval, xoval, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'xval, xoval', 'shutdown_hf')

    deallocate (scx1, scx2, scx3, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'scx?', 'shutdown_hf')

    deallocate (wrkx, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'wrkx', 'shutdown_hf')

    deallocate (rwrkx, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'rwrkx', 'shutdown_hf')
  end if


  return
end subroutine shutdown_hf



subroutine hamilt_hf (root, procid, iout, outfile, iprint, imethd, iwfnty, &
     & norb, nbas, nbct, nup, ndn, xmat, hmat, i2sv, ni2s, lsf, enr, ndspin, &
     & ndsptl, ndcplx, darr, xdim, ndet, ndtt, ovdet, hmdet, nstat, idetst, &
     & mxdtst, inovst, inhmst, mdvec, iopt, ndtopt, imult, ngrda, ngrdb, &
     & ngrdg, grda, grdb, grdg, wgta, wgtb, wgtg, nop, rotmat, wgtmat, &
     & energy, garr, gdim, thrsh, iflg)

  use i2sint
  use spin

! +----------------------------------------------------------------+
! |                                                                |
! | hamilt_hf  --  CAJH, 11.2012                                   |
! |                                                                |
! |                                                                |
! | Evaluate the energy of a linear combination (FED or RES) of    |
! | symmetry-projected Hartree-Fock states. The local gradient     |
! | wrt the underlying determinants is also computed.              |
! |                                                                |
! | This subroutine supports both spatial and spin projection with |
! | all types of HF determinants (RHF, UHF, GHF).                  |
! |                                                                |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful evaluation                           |
! |         =  1,  orthogonality between some dets found           |
! |         =  2,  diagonalization in current state failed         |
! |                                                                |
! +----------------------------------------------------------------+

  include 'mpif.h'

  ! input variables

  !   root    - id of root process
  !   procid  - id of each process
  !   iout    - output file unit number
  !   outfile - output file name
  !   iprint  - printing level

  integer, intent(in) :: root, procid, iout, iprint

  character(len=*), intent(in) :: outfile

  !   imethd  - PHF method to use
  !   iwfnty  - type of wavefunction to use
  !   norb    - number of orbitals
  !   nbas    - number of basis functions
  !   nbct    - number of Cartesian basis functions
  !   nup     - number of spin-up electrons
  !   ndn     - number of spin-dn electrons
  !   xmat    - transformation matrix [ = S^(-1/2) ]
  !   hmat    - core Hamiltonian matrix ( ort AO basis )

  integer, intent(in) :: imethd, iwfnty, norb, nbas, nbct, nup, ndn

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,norb), intent(in) :: hmat

  !   i2sv    - vector of 2-electron integrals
  !   ni2s    - number of 2-electron integrals stored in i2sv
  !   lsf     - whether stored integrals include symmetry factor
  !   enr     - nuclear repulsion energy

  integer, intent(in) :: ni2s
  logical, intent(in) :: lsf

  real(kind=dp) :: enr
  !jag20
  real(kind=dp) :: gr, gi

  type (i2s), dimension(ni2s), intent(in) :: i2sv

  !   ndspin  - number of configs per bs determinant (spin)
  !   ndsptl  - number of configs per bs determinant (spatial)
  !   ndcplx  - number of configs per bs determinant (complex)

  integer, intent(in) :: ndspin, ndsptl, ndcplx

  !   darr    - array of determinants
  !   xdim    - leading dimension of array darr
  !   ndet    - total number of determinants stored
  !   ndtt    - dimension of ovdet, hmdet
  !             (ndtt = ndet*ndspin*ndsptl*ndcplx)
  !   ovdet   - overlap matrix among dets [updated]
  !   hmdet   - hamiltonian matrix among dets [updated]

  integer, intent(in) :: xdim, ndet, ndtt

  complex(kind=dp), dimension(xdim,ndet), intent(in) :: darr
  complex(kind=dp), dimension(ndtt,ndtt), intent(inout) :: ovdet, hmdet

  !   nstat   - number of states in wavefunction
  !   idetst  - array with number of determinants per state
  !   mxdtst  - maximum number of determinants in a state
  !   inovst  - inverse overlap among nstat-1 previous states
  !   inhmst  - S-1 . H . S-1 for nstat-1 previous states
  !   mdvec   - CI expansion of states among its dets [updated]

  integer, intent(in) :: nstat, mxdtst
  integer, dimension(nstat), intent(in) :: idetst

  complex(kind=dp), dimension(:,:), intent(in) :: inovst, inhmst
  complex(kind=dp), dimension(ndspin*ndsptl*ndcplx*mxdtst,nstat), &
       & intent(inout) :: mdvec

  !   iopt    - type of optimization (1 - FED, 2 - RES)
  !   ndtopt  - number of determinants to update

  integer, intent(in) :: iopt, ndtopt

  !   imult   - multiplicity to use in spin projection
  !   ngrda   - number of grid points (alpha)
  !   ngrdb   - number of grid points (beta )
  !   ngrdg   - number of grid points (gamma)
  !   grda    - array of grid points (alpha)
  !   grdb    - array of grid points (beta )
  !   grdg    - array of grid points (gamma)
  !   wgta    - array of integration weights (alpha)
  !   wgtb    - array of integration weights (beta )
  !   wgtg    - array of integration weights (gamma)

  integer, intent(in) :: imult, ngrda, ngrdb, ngrdg

  real(kind=dp), dimension(ngrda), intent(in) :: grda, wgta
  real(kind=dp), dimension(ngrdb), intent(in) :: grdb, wgtb
  real(kind=dp), dimension(ngrdg), intent(in) :: grdg, wgtg

  !   nop     - number of symmetry operations for spatial projection
  !   rotmat  - rotation matrices for spatial projection ( ort AO basis )
  !   wgtmat  - weight matrices for spatial projection

  integer, intent(in) :: nop

  complex(kind=dp), dimension(ndsptl,ndsptl,nop), intent(in) :: wgtmat
  complex(kind=dp), dimension(norb,norb,nop), intent(in) :: rotmat

  !   thrsh   - threshold to use in determining linear independency

  real(kind=dp), intent(in) :: thrsh


  ! output variables

  !   energy  - energy of current state
  !   garr    - gradient of energy per determinant optimized
  !   gdim    - dimension of each gradient vector
  !   iflg    - return code

  integer, intent(in) :: gdim
  integer, intent(out) :: iflg

  real(kind=dp), intent(out) :: energy

  complex(kind=dp), dimension(gdim,ndtopt), intent(inout) :: garr


  ! other variables

  integer :: ierr, idum, ldim, iflg1
  integer :: j, ik, p1a, p1b, p2a, p2b
  integer :: nosq, ntot, iwfnt1, ndt1, ndt2, ngrds
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: ind1, ind2, ind3, ind4, k1, k2, k3
  integer :: ind1_old, ind2_old, ind3_old, k1_old, k2_old
  integer :: aind, bind, gind
  integer :: jval, mval, nval
  integer :: ih, ip
  integer :: j1, j2, j1a, j1b, j2a, j2b, ik1, ik2, ik3
  integer :: ik1b, ik2b
  logical :: lint

  integer :: ij1, ij2, ix1, ix2
  integer :: ist1, ist2, idet1, idet2
  integer :: iz1, iz2, izx1, izx2, izy2
  integer :: iz1a, iz1b, iz1r

  integer :: mint, mdc
  parameter ( mint = mpi_integer, mdc = mpi_double_complex )

  integer :: mworld
  parameter ( mworld = mpi_comm_world )

  real(kind=dp) :: alpha, beta, gamma, wgtY

  complex(kind=dp) :: ovpX, hmovpX, wgtX
  complex(kind=dp) :: ztr, ztr_up, ztr_dn


  iflg = 0

  lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
  lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
  lspin = mod(imethd-1,2) == 1

  lsuhf = .false.
  lsghf = .false.

  if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
  if ( lspin .and. iwfnty == 3 )  lsghf = .true.


  ! set some pointers

    ik = 0
  do j = 1, nstat-1
    ik = ik + idetst(j)
  end do

    p1a = ik+1;  p1b = 1    ! index of 1st det in current state
    p2a = ik+1;  p2b = 1    ! index of 1st det to optimize

  if ( iopt == 1 ) then
    p2a = p2a + idetst(nstat)-1   ! adjust for FED
    p2b = p2b + idetst(nstat)-1
  end if


  ! set some important variables

  nosq = norb*norb
  ntot = nup + ndn

  iwfnt1 = iwfnty
  if ( lsuhf ) iwfnt1 = 3

  ndt1  = ndspin*ndsptl*ndcplx*ndtopt
  ndt2  = ndspin*ndsptl*ndcplx*idetst(nstat)
  ngrds = ngrda*ngrdb*ngrdg

  ! decide spin to use in projection

  if ( lspin .and. mod(imult,2) == 0 ) then
    jval = imult-1
    lint = .false.
  else if ( lspin ) then
    jval = (imult-1)/2
    lint = .true.
  end if


  ovlpr(1:ndt1,1:ndtt) = z0
  hmltr(1:ndt1,1:ndtt) = z0

  garr1r(1:gdim,1:ndt1,1:ndtt) = z0
  garr2r(1:gdim,1:ndt1,1:ndtt) = z0


  ! set overlap and Hamiltonian arrays to 0

  !$omp  parallel default(private) &
  !$omp& shared(lcplx, lsptl, lspin, lsuhf, lsghf, jval, lint, &
  !$omp&        nosq, ntot, ndt1, ndt2, ndtt, ngrds, nop, &
  !$omp&        ndspin, ndsptl, ndcplx, &
  !$omp&        iwfnty, iwfnt1, nup, ndn, &
  !$omp&        gdim, norb, nbas, nbct, &
  !$omp&        ngrda, ngrdb, ngrdg, grda, grdb, grdg, &
  !$omp&        wgta, wgtb, wgtg, wgtmat, &
  !$omp&        darr, rotmat, i2sv, ni2s, lsf, xmat, hmat, &
  !$omp&        jstrt, ntask, ovlpr, hmltr, garr1r, garr2r, iflg1)

  ovlp(1:ndt1,1:ndtt) = z0
  hmlt(1:ndt1,1:ndtt) = z0


  ! set gradient to zero

  garr1(1:gdim,1:ndt1,1:ndtt) = z0
  garr2(1:gdim,1:ndt1,1:ndtt) = z0


  ! loop over grid points / determinants to evaluate overlap and
  ! Hamiltonian matrix elements, as well as the gradient wrt reference
  ! determinants

    iflg1 = 0

    aind = 0
    bind = 0
    gind = 0

    ind1_old = 0
    ind2_old = 0
    ind3_old = 0
    k1_old   = 0
    k2_old   = 0

  !$omp  do schedule(static,1) reduction(+:iflg1)

  do j = jstrt, jstrt+ntask-1

    ind1 = mod(j-1,ngrds) + 1                        ! spin grid index
    ind2 = mod((j-1)/ngrds,nop) + 1                  ! spatial symm index
    ind3 = mod((j-1)/(ngrds*nop),ndcplx) + 1         ! complex left index
    ind4 = (j-1)/(ngrds*nop*ndcplx) + 1              ! det index

    if ( lspin ) then
      gind = mod(ind1-1,ngrdg) + 1
      bind = mod((ind1-1)/ngrdg,ngrdb) + 1
      aind = (ind1-1)/(ngrdg*ngrdb) + 1

      alpha = grda(aind)
      beta  = grdb(bind)
      gamma = grdg(gind)
    end if

    k1 = indx(ind4,1)    ! index of left state wrt all dets
    k2 = indx(ind4,2)    ! index of right state wrt all dets
    k3 = indx(ind4,3)    ! index of left state wrt dets to optimize


    ! load determinants (we actually load D*)
    ! multiply D2 by spatial rotation matrices

    if ( k1 /= k1_old .or. ind3 /= ind3_old ) then
      if ( iwfnty == 1 ) then
        if ( ind3 == 1 ) then
          d1st(1:nosq) = conjg (darr(1:nosq,k1))
        else
          d1st(1:nosq) = darr(1:nosq,k1)
        end if

      else if ( iwfnty == 2 ) then
        if ( ind3 == 1 ) then
          d1st(1:2*nosq) = conjg (darr(1:2*nosq,k1))
        else
          d1st(1:2*nosq) = darr(1:2*nosq,k1)
        end if

      else if ( iwfnty == 3 ) then
        if ( ind3 == 1 ) then
          d1st(1:4*nosq) = conjg (darr(1:4*nosq,k1))
        else
          d1st(1:4*nosq) = darr(1:4*nosq,k1)
        end if
      end if
    end if

    if ( k2 /= k2_old ) then
      if ( iwfnty == 1 ) then
        if ( .not. lsptl ) then
          d2st(1:nosq) = conjg (darr(1:nosq,k2))
        else
          d2stY(1:nosq) = conjg (darr(1:nosq,k2))

          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY, norb, z0, d2st, norb)
        end if

      else if ( iwfnty == 2 ) then
        if ( .not. lsptl .and. .not. lspin ) then
          d2st(1:2*nosq) = conjg (darr(1:2*nosq,k2))
        else if ( .not. lsptl .and. lspin ) then
          d2stX(1:2*nosq) = conjg (darr(1:2*nosq,k2))
        else
          d2stY(1:2*nosq) = conjg (darr(1:2*nosq,k2))

          if ( lspin ) then
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(1), norb, z0, d2stX(1), norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(nosq+1), norb, z0, d2stX(nosq+1), norb)

          else
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(1), norb, z0, d2st(1), norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(nosq+1), norb, z0, d2st(nosq+1), norb)
          end if
        end if

      else if ( iwfnty == 3 ) then
        if ( .not. lsptl .and. .not. lspin ) then
          d2st(1:4*nosq) = conjg (darr(1:4*nosq,k2))
        else if ( .not. lsptl .and. lspin ) then
          d2stX(1:4*nosq) = conjg (darr(1:4*nosq,k2))
        else
          d2stY(1:4*nosq) = conjg (darr(1:4*nosq,k2))

          if ( lspin ) then
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(1), 2*norb, z0, d2stX(1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(norb+1), 2*norb, z0, d2stX(norb+1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(2*nosq+1), 2*norb, z0, d2stX(2*nosq+1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(2*nosq+norb+1), 2*norb, z0, d2stX(2*nosq+norb+1), 2*norb)

          else
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(1), 2*norb, z0, d2st(1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(norb+1), 2*norb, z0, d2st(norb+1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(2*nosq+1), 2*norb, z0, d2st(2*nosq+1), 2*norb)
          call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
               & d2stY(2*nosq+norb+1), 2*norb, z0, d2st(2*nosq+norb+1), 2*norb)
          end if
        end if
      end if

    else if ( k2 == k2_old .and. ind2 /= ind2_old ) then
      if ( iwfnty == 1 ) then
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY, norb, z0, d2st, norb)

      else if ( iwfnty == 2 ) then
        if ( lspin ) then
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(1), norb, z0, d2stX(1), norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(nosq+1), norb, z0, d2stX(nosq+1), norb)

        else
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(1), norb, z0, d2st(1), norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(nosq+1), norb, z0, d2st(nosq+1), norb)
        end if

      else if ( iwfnty == 3 ) then
        if ( lspin ) then
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(1), 2*norb, z0, d2stX(1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(norb+1), 2*norb, z0, d2stX(norb+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(2*nosq+1), 2*norb, z0, d2stX(2*nosq+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(2*nosq+norb+1), 2*norb, z0, d2stX(2*nosq+norb+1), 2*norb)

        else
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(1), 2*norb, z0, d2st(1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(norb+1), 2*norb, z0, d2st(norb+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(2*nosq+1), 2*norb, z0, d2st(2*nosq+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & d2stY(2*nosq+norb+1), 2*norb, z0, d2st(2*nosq+norb+1), 2*norb)
        end if
      end if
    end if


    ! form full spin-orbital matrices for SUHF

    if ( lsuhf .and. &
       & ( k1 /= k1_old .or. ind3 /= ind3_old ) ) then

      call uhf_to_ghf (norb, nup, ndn, d1st)
    end if

    if ( lsuhf .and. &
       & ( k2 /= k2_old .or. &
         & ( k2 == k2_old .and. ind2 /= ind2_old ) ) ) then

      call uhf_to_ghf (norb, nup, ndn, d2stX)
    end if


    ! multiply D2 by spin rotation matrix

    if ( lspin ) then
      call spinrot (norb, alpha, beta, gamma, d2stX, d2st, 1)
    end if


    ! compute overlap
    ! build transition density matrix

    if ( iwfnt1 == 1 ) then
      call trdnmt_rhf (norb, nup, d1st, d2st, rho, ovpX, iflg1)

    else if ( iwfnt1 == 2 ) then
      call trdnmt_uhf (norb, nup, ndn, d1st, d2st, rho, ovpX, iflg1)

    else if ( iwfnt1 == 3 ) then
      call trdnmt_ghf (norb, nup, ndn, d1st, d2st, rho, ovpX, iflg1)
    end if

    if ( iflg1 /= 0 ) cycle


    ! contract transition density matrix with two-particle integrals

    if ( iwfnt1 == 1 ) then
      call erictr_rhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

    else if ( iwfnt1 == 2 ) then
      call erictr_uhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

    else if ( iwfnt1 == 3 ) then
      call erictr_ghf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)
    end if


    ! compute Hamiltonian overlap

    hmovpX = z0

    !   .. first the single-particle contribution ..

    if ( iwfnt1 == 1 ) then
      ztr = trab_zmat (norb, hmat, norb, rho, norb, idum)
      hmovpX = hmovpX + d2*ztr

    else if ( iwfnt1 == 2 ) then
      ztr_up = trab_zmat (norb, hmat, norb, rho(1,1), norb, idum)
      ztr_dn = trab_zmat (norb, hmat, norb, rho(1,norb+1), norb, idum)
      hmovpX = hmovpX + ztr_up + ztr_dn

    else if ( iwfnt1 == 3 ) then
      ztr_up = trab_zmat (norb, hmat, norb, rho(1,1), 2*norb, idum)
      ztr_dn = trab_zmat (norb, hmat, norb, rho(norb+1,norb+1), 2*norb, idum)
      hmovpX = hmovpX + ztr_up + ztr_dn
    end if

    !   .. now the two-particle contribution ..

    if ( iwfnt1 == 1 ) then
      ztr = trab_zmat (norb, rho, norb, gam, norb, idum)
      hmovpX = hmovpX + ztr

    else if ( iwfnt1 == 2 ) then
      ztr_up = trab_zmat (norb, rho(1,1), norb, gam(1,1), norb, idum)
      ztr_dn = trab_zmat (norb, rho(1,norb+1), norb, gam(1,norb+1), norb, idum)
      hmovpX = hmovpX + 0.50e0_dp * (ztr_up + ztr_dn)

    else if ( iwfnt1 == 3 ) then
      ztr = trab_zmat (2*norb, rho, 2*norb, gam, 2*norb, idum)
      hmovpX = hmovpX + 0.50e0_dp * ztr
    end if


    ! build tpgam

    if ( iwfnt1 == 1 ) then
      tpgam(1:norb,1:norb) = hmat(1:norb,1:norb) &
                         & + gam(1:norb,1:norb)

    else if ( iwfnt1 == 2 ) then
      tpgam(1:norb,1:norb) = hmat(1:norb,1:norb) &
                         & + gam(1:norb,1:norb)

      tpgam(1:norb,norb+1:2*norb) = hmat(1:norb,1:norb) &
                                & + gam(1:norb,norb+1:2*norb)

    else if ( iwfnt1 == 3 ) then
      tpgam(1:2*norb,1:2*norb) = gam(1:2*norb,1:2*norb)

      tpgam(1:norb,1:norb) &
           & = tpgam(1:norb,1:norb) &
           & + hmat(1:norb,1:norb)
      tpgam(norb+1:2*norb,norb+1:2*norb) &
           & = tpgam(norb+1:2*norb,norb+1:2*norb) &
           & + hmat(1:norb,1:norb)
    end if


    ! build contributions to gradient

    if ( iwfnt1 == 1 ) then
      call bldgrad_rhf (norb, nup)
    else if ( iwfnt1 == 2 ) then
      call bldgrad_uhf (norb, nup, ndn)
    else if ( iwfnt1 == 3 ) then
      call bldgrad_ghf (norb, nup, ndn)
    end if


    ! reshape contributions to gradient

    !   gr1t = (hmovpX * gr1sc + gr2sc)
    !   gr2t = gr1sc

    if ( iwfnty == 2 .and. lsuhf ) then
        iz1 = 0; iz2 = 0

      do ih = 1, nup
        do ip = 1, norb-nup
          gr1sc(iz1+ip) = gr1sc(iz2+ip)
          gr2sc(iz1+ip) = gr2sc(iz2+ip)
        end do

        iz1 = iz1 + norb-nup
        iz2 = iz2 + 2*norb-ntot
      end do

        iz2 = iz2 + norb-nup

      do ih = 1, ndn
        do ip = 1, norb-ndn
          gr1sc(iz1+ip) = gr1sc(iz2+ip)
          gr2sc(iz1+ip) = gr2sc(iz2+ip)
        end do

        iz1 = iz1 + norb-ndn
        iz2 = iz2 + 2*norb-ntot
      end do
    end if

    if ( iwfnty == 1 ) then
      gr1t(1:gdim) = d2*hmovpX * gr1sc(1:gdim) + d2*gr2sc(1:gdim)
      gr2t(1:gdim) = d2*gr1sc(1:gdim)

    else if ( iwfnty == 2 .or. iwfnty == 3 ) then
      gr1t(1:gdim) = hmovpX * gr1sc(1:gdim) + gr2sc(1:gdim)
      gr2t(1:gdim) = gr1sc(1:gdim)
    end if


    ! loop over ndspin, ndsptl to finish accumulating quantities

    do j1 = 1, ndspin*ndsptl
      j1a = mod(j1-1,ndspin) + 1
      j1b = (j1-1)/ndspin + 1

      ik1  = (k3-1)*ndspin*ndsptl*ndcplx &
         & + (ind3-1)*ndspin*ndsptl + j1
      ik1b = ik1

      if ( lcplx .and. ind3 == 1 ) ik1b = ik1 + ndspin*ndsptl
      if ( lcplx .and. ind3 == 2 ) ik1b = ik1 - ndspin*ndsptl

      do j2 = 1, ndspin*ndsptl
        j2a = mod(j2-1,ndspin) + 1
        j2b = (j2-1)/ndspin + 1

        ik2  = (k2-1)*ndspin*ndsptl*ndcplx + j2
        ik2b = ik2

        if ( lcplx ) ik2b = ik2 + ndspin*ndsptl

        if ( lsghf .and. .not. lint ) then
          mval = jval-2*(j1a-1)
          nval = jval-2*(j2a-1)
        else if ( lsghf ) then
          mval = jval-(j1a-1)
          nval = jval-(j2a-1)
        end if

        if ( lsuhf .and. .not. lint ) then
          mval = nup-ndn
          nval = mval
        else if ( lsuhf ) then
          mval = (nup-ndn)/2
          nval = mval
        end if

        wgtX = z1

        if ( lspin ) then
          call wignerd (wgtY, lint, jval, mval, nval, beta, idum)

          if ( .not. lint ) then
            wgtX = wgtX * wgtY &
               & * exp(zi * real(mval,dp)/d2 * alpha) &
               & * exp(zi * real(nval,dp)/d2 * gamma)
          else
            wgtX = wgtX * wgtY &
               & * exp(zi * real(mval,dp) * alpha) &
               & * exp(zi * real(nval,dp) * gamma)
          end if

          wgtX = wgtX * wgta(aind) * wgtb(bind) * wgtg(gind)
        end if

        if ( lsptl ) then
          wgtX = wgtX * conjg(wgtmat(j1b,j2b,ind2))
        end if

        ovlp(ik1,ik2) = ovlp(ik1,ik2) &
                    & + wgtX * ovpX
        hmlt(ik1,ik2) = hmlt(ik1,ik2) &
                    & + wgtX * ovpX * hmovpX

        garr1(1:gdim,ik1,ik2) = garr1(1:gdim,ik1,ik2) &
                            & + wgtX * ovpX * gr1t(1:gdim)
        garr2(1:gdim,ik1,ik2) = garr2(1:gdim,ik1,ik2) &
                            & + wgtX * ovpX * gr2t(1:gdim)

        if ( lcplx ) then
          ovlp(ik1b,ik2b) = ovlp(ik1b,ik2b) &
                        & + conjg(wgtX * ovpX)
          hmlt(ik1b,ik2b) = hmlt(ik1b,ik2b) &
                        & + conjg(wgtX * ovpX * hmovpX)

          garr1(1:gdim,ik1b,ik2b) = garr1(1:gdim,ik1b,ik2b) &
                                & + conjg(wgtX * ovpX * gr1t(1:gdim))
          garr2(1:gdim,ik1b,ik2b) = garr2(1:gdim,ik1b,ik2b) &
                                & + conjg(wgtX * ovpX * gr2t(1:gdim))
        end if
      end do
    end do

    ind1_old = ind1
    ind2_old = ind2
    ind3_old = ind3
    k1_old   = k1
    k2_old   = k2
  end do

  !$omp  end do


  !$omp  critical

    ovlpr(1:ndt1,1:ndtt) = ovlpr(1:ndt1,1:ndtt) + ovlp(1:ndt1,1:ndtt)
    hmltr(1:ndt1,1:ndtt) = hmltr(1:ndt1,1:ndtt) + hmlt(1:ndt1,1:ndtt)

    garr1r(1:gdim,1:ndt1,1:ndtt) = garr1r(1:gdim,1:ndt1,1:ndtt) + &
         & garr1(1:gdim,1:ndt1,1:ndtt)
    garr2r(1:gdim,1:ndt1,1:ndtt) = garr2r(1:gdim,1:ndt1,1:ndtt) + &
         & garr2(1:gdim,1:ndt1,1:ndtt)

  !$omp  end critical
  !$omp  barrier

  !$omp  end parallel


  ! check if any processor had a non-zero flag

  ! call mpi_allreduce (iflg1, iflg2, 1, mint, mpi_sum, mworld, ierr)

  if ( iflg1 /= 0 ) then
    write (6, *) 'error: Energy evaluation failed.'
    stop
  end if
 

  ! collect all contributions to overlap, Hamiltonian matrices
  ! collect all contributions to gradient

  if ( procid == root ) then
    call mpi_reduce (mpi_in_place, ovlpr(1,1), ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (mpi_in_place, hmltr(1,1), ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (mpi_in_place, garr1r(1,1,1), gdim*ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (mpi_in_place, garr2r(1,1,1), gdim*ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
  else
    call mpi_reduce (ovlpr(1,1), idum, ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (hmltr(1,1), idum, ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (garr1r(1,1,1), idum, gdim*ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
    call mpi_reduce (garr2r(1,1,1), idum, gdim*ndt1*ndtt, mdc, &
         & mpi_sum, root, mworld, ierr)
  end if


    ! update ovdet, hmdet

  if ( procid == root ) then
    do j = 1, ndet*ndtopt
      k1 = indx(j,1)
      k2 = indx(j,2)
      k3 = indx(j,3)

      do j1 = 1, ndspin*ndsptl*ndcplx
        ik1 = (k1-1)*ndspin*ndsptl*ndcplx + j1
        ik3 = (k3-1)*ndspin*ndsptl*ndcplx + j1

        do j2 = 1, ndspin*ndsptl*ndcplx
          ik2 = (k2-1)*ndspin*ndsptl*ndcplx + j2

          ovdet(ik1,ik2) = ovlpr(ik3,ik2)
          hmdet(ik1,ik2) = hmltr(ik3,ik2)

          if ( k2 < p2a ) then
            ovdet(ik2,ik1) = conjg (ovdet(ik1,ik2))
            hmdet(ik2,ik1) = conjg (hmdet(ik1,ik2))
          end if
        end do
      end do
    end do
  end if


  ! form Hamiltonian and overlap matrices to diagonalize

  if ( procid == root ) then
    ovmt(1:ndt2,1:ndt2) = z0
    hmmt(1:ndt2,1:ndt2) = z0

    ! contributions from current state dets

      izx1 = (p1a-1)*ndspin*ndsptl*ndcplx
    do iz1 = 1, ndt2
      izx1 = izx1 + 1
        izx2 = (p1a-1)*ndspin*ndsptl*ndcplx
      do iz2 = 1, ndt2
        izx2 = izx2 + 1

        ovmt(iz1,iz2) = ovmt(iz1,iz2) + ovdet(izx1,izx2)
        hmmt(iz1,iz2) = hmmt(iz1,iz2) + hmdet(izx1,izx2)
      end do
    end do

    ! contributions from other state dets

    if ( nstat > 1 ) then
      ij1 = 0
    do ist1 = 1, nstat-1
      idet1 = idetst(ist1)
      do ix1 = 1, idet1*ndspin*ndsptl*ndcplx
        ij1 = ij1 + 1

          ij2 = 0
        do ist2 = 1, nstat-1
          idet2 = idetst(ist2)
          do ix2 = 1, idet2*ndspin*ndsptl*ndcplx
            ij2 = ij2 + 1

              izx1 = (p1a-1)*ndspin*ndsptl*ndcplx
            do iz1 = 1, ndt2
              izx1 = izx1 + 1
                izx2 = (p1a-1)*ndspin*ndsptl*ndcplx
              do iz2 = 1, ndt2
                izx2 = izx2 + 1

          ovmt(iz1,iz2) = ovmt(iz1,iz2) &
             & - mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & inovst(ist1,ist2) * &
               & ovdet(izx1,ij1) * ovdet(ij2,izx2)

          hmmt(iz1,iz2) = hmmt(iz1,iz2) &
             & - mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & inovst(ist1,ist2) * &
               & ( ovdet(izx1,ij1) * hmdet(ij2,izx2) + &
                 & hmdet(izx1,ij1) * ovdet(ij2,izx2) ) &
             & + mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & inhmst(ist1,ist2) * &
               & ovdet(izx1,ij1) * ovdet(ij2,izx2)
              end do
            end do
          end do
        end do
      end do
    end do
    end if

    do iz1 = 1, ndt2
      ovmt(iz1,iz1) = real (ovmt(iz1,iz1))
      hmmt(iz1,iz1) = real (hmmt(iz1,iz1))
    end do
  end if


  ! solve CI problem

  if ( procid == root ) then
    call zgenevs ('u', ndt2, hmmt, ovmt, ldim, xvec, xval, xoval, &
       & scx1, scx2, scx3, wrkx, lwrkx, rwrkx, thrsh, iflg1)

    if ( iflg1 /= 0 ) then
      iflg = 2
    end if

    ! save energy and CI vector

    energy = xval(1)
    mdvec(1:ndt2,nstat) = xvec(1:ndt2,1)

    ! add nuclear repulsion energy

    xval(1:ldim) = xval(1:ldim) + enr


    if ( iprint > 1 ) then
      open (unit = iout, file = outfile, status = 'old', &
          & form = 'formatted', position = 'append')

      if ( iprint >= 2 ) then

        ! print eigenvalues

        write (iout, *)
        call print_dmat (iout, 1, ndt2, 1, 'hmmt eigenvalues', xval)
        call print_dmat (iout, 1, ndt2, 1, 'ovmt eigenvalues', xoval)
      end if

      if ( iprint >= 3 ) then

        ! print Hamiltonian and overlap matrix

        call print_zmat (iout, 1, ndt2, ndt2, 'hmmt matrix', hmmt)
        call print_zmat (iout, 1, ndt2, ndt2, 'ovmt matrix', ovmt)
      end if

      close (unit = iout)
    end if

  else
    energy = d0
  end if

  ! broadcast iflg

  call mpi_bcast (iflg, 1, mint, root, mworld, ierr)

  if ( iflg /= 0 ) return


  ! finish building gradient

  if ( procid == root ) then
    garr(1:gdim,1:ndtopt) = z0

    ! contributions from current state dets

      izx1 = (p2b-1)*ndspin*ndsptl*ndcplx
    do iz1a = 1, ndtopt
    do iz1r = 1, ndcplx
    do iz1b = 1, ndspin*ndsptl
      izx1 = izx1 + 1
      iz1  = (iz1a-1)*ndspin*ndsptl*ndcplx &
         & + (iz1r-1)*ndspin*ndsptl + iz1b

        izx2 = (p1b-1)*ndspin*ndsptl*ndcplx
        izy2 = (p1a-1)*ndspin*ndsptl*ndcplx
      do iz2 = 1, ndt2
        izx2 = izx2 + 1
        izy2 = izy2 + 1

        if ( iz1r == 1 ) then
          garr(1:gdim,iz1a) = garr(1:gdim,iz1a) &
               & + conjg(mdvec(izx1,nstat)) * mdvec(izx2,nstat) * &
               &   ( garr1r(1:gdim,iz1,izy2) - &
               &     energy * garr2r(1:gdim,iz1,izy2) )

        else
          garr(1:gdim,iz1a) = garr(1:gdim,iz1a) &
               & + mdvec(izx1,nstat) * conjg(mdvec(izx2,nstat)) * &
               &   ( conjg(garr1r(1:gdim,iz1,izy2)) - &
               &     energy * conjg(garr2r(1:gdim,iz1,izy2)) )
        end if
      end do
    end do
    end do
    end do

    ! contributions from other state dets

    if ( nstat > 1 ) then
      ij1 = 0
    do ist1 = 1, nstat-1
      idet1 = idetst(ist1)
      do ix1 = 1, idet1*ndspin*ndsptl*ndcplx
        ij1 = ij1 + 1

          ij2 = 0
        do ist2 = 1, nstat-1
          idet2 = idetst(ist2)
          do ix2 = 1, idet2*ndspin*ndsptl*ndcplx
            ij2 = ij2 + 1

              izx1 = (p2b-1)*ndspin*ndsptl*ndcplx
            do iz1a = 1, ndtopt
            do iz1r = 1, ndcplx
            do iz1b = 1, ndspin*ndsptl
              izx1 = izx1 + 1
              iz1  = (iz1a-1)*ndspin*ndsptl*ndcplx &
                   & + (iz1r-1)*ndspin*ndsptl + iz1b

                izx2 = (p1b-1)*ndspin*ndsptl*ndcplx
                izy2 = (p1a-1)*ndspin*ndsptl*ndcplx
              do iz2 = 1, ndt2
                izx2 = izx2 + 1
                izy2 = izy2 + 1

            if ( iz1r == 1 ) then
        garr(1:gdim,iz1a) = garr(1:gdim,iz1a) &
             & - conjg(mdvec(izx1,nstat)) * mdvec(izx2,nstat) * &
             &   mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & inovst(ist1,ist2) * &
               & ( garr2r(1:gdim,iz1,ij1) * hmdet(ij2,izy2) - &
               &   energy * garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2) ) &
             & - conjg(mdvec(izx1,nstat)) * mdvec(izx2,nstat) * &
             &   mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & inovst(ist1,ist2) * &
               & ( garr1r(1:gdim,iz1,ij1) * ovdet(ij2,izy2) - &
               &   energy * garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2) ) &
             & + conjg(mdvec(izx1,nstat)) * mdvec(izx2,nstat) * &
             &   mdvec(ix1,ist1) * conjg(mdvec(ix2,ist2)) * &
               & ( inhmst(ist1,ist2) * &
               &   garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2) - &
               &   inovst(ist1,ist2) * &
               &   energy * garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2) )

            else
        garr(1:gdim,iz1a) = garr(1:gdim,iz1a) &
             & - mdvec(izx1,nstat) * conjg(mdvec(izx2,nstat)) * &
             &   conjg(mdvec(ix1,ist1)) * mdvec(ix2,ist2) * &
               & conjg(inovst(ist1,ist2)) * &
               & ( conjg(garr2r(1:gdim,iz1,ij1) * hmdet(ij2,izy2)) - &
               &   energy * conjg(garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2)) ) &
             & - mdvec(izx1,nstat) * conjg(mdvec(izx2,nstat)) * &
             &   conjg(mdvec(ix1,ist1)) * mdvec(ix2,ist2) * &
               & conjg(inovst(ist1,ist2)) * &
               & ( conjg(garr1r(1:gdim,iz1,ij1) * ovdet(ij2,izy2)) - &
               &   energy * conjg(garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2)) ) &
             & + mdvec(izx1,nstat) * conjg(mdvec(izx2,nstat)) * &
             &   conjg(mdvec(ix1,ist1)) * mdvec(ix2,ist2) * &
               & ( conjg(inhmst(ist1,ist2)) * &
               &   conjg(garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2)) - &
               &   conjg(inovst(ist1,ist2)) * &
               &   energy * conjg(garr2r(1:gdim,iz1,ij1) * ovdet(ij2,izy2)) )
            end if
              end do
            end do
            end do
            end do
          end do
        end do
      end do
    end do
    end if
  end if

  !jag20
if ( procid == root ) then
    do ik1 = 1, ndtopt
      do ik2 = 1, gdim
        if ( abs(real(garr(ik2,ik1))) > 1.0e-12_dp ) then
          gr = real(garr(ik2,ik1))
        else
          gr = d0
        end if
        if ( abs(aimag(garr(ik2,ik1))) > 1.0e-12_dp ) then
          gi = aimag(garr(ik2,ik1))
        else
          gi = d0
        end if
        garr(ik2,ik1) = cmplx(gr,gi,dp)
      end do
    end do
end if
 

  return
end subroutine hamilt_hf



subroutine setup_bldgrad (iwfnty, norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_bldgrad  --  CAJH, 11.2012                               |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in bldgrad.        |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons

  integer, intent(in) :: iwfnty, norb, nup, ndn


  ! other variables

  integer :: istatus, ntot


  istatus = 0
  ntot = nup + ndn

  ! allocate space for scr1, scr2

  if ( iwfnty == 1 ) then
    allocate (scr2(norb,nup), scr1(norb,nup), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (scr2(norb,ntot), scr1(norb,ntot), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scr2(2*norb,ntot), scr1(2*norb,ntot), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'setup_bldgrad')

  return
end subroutine setup_bldgrad



subroutine shutdown_bldgrad

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_bldgrad  --  CAJH, 11.2012                            |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in bldgrad.      |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: istatus


  ! deallocate space for scr?

  deallocate (scr1, scr2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'shutdown_bldgrad')


  return
end subroutine shutdown_bldgrad



subroutine bldgrad_ghf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | bldgrad_ghf  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Build contributions to the gradient of a Hamiltonian           |
! | expectation value wrt |phi1>. We here build                    |
! |                                                                |
! |   gr1sc  =  D1(:,N+1:end)^T . rho . D1(:,1:N)*                 |
! |   gr2sc  =  D1(:,N+1:end)^T . (I - rho) .                      |
! |                                (T + Gam) . rho . D1(:,1:N)*    |
! |                                                                |
! | On input, T + Gam should be loaded in tpgam.                   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: k, ntot


  ntot = nup + ndn

  ! compute:  scr1 = rho . D1(:,1:N)*

  call zgemm ('n', 'n', 2*norb, ntot, 2*norb, z1, rho, 2*norb, &
       & d1st, 2*norb, z0, scr1, 2*norb)

  ! compute:  gr1sc = D1(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', 2*norb-ntot, ntot, 2*norb, z1, d1st(2*norb*ntot+1), &
       & 2*norb, scr1, 2*norb, z0, gr1sc, 2*norb-ntot)


  ! build rho-I; use rho as scratch space

  do k = 1, 2*norb
    rho(k,k) = rho(k,k) - z1
  end do

  ! compute:  scr2 = (T + Gam) . scr1

  call zgemm ('n', 'n', 2*norb, ntot, 2*norb, z1, tpgam, 2*norb, &
       & scr1, 2*norb, z0, scr2, 2*norb)

  ! compute:  scr1 = (I - rho) . scr2

  call zgemm ('n', 'n', 2*norb, ntot, 2*norb, -z1, rho, 2*norb, &
       & scr2, 2*norb, z0, scr1, 2*norb)

  ! compute:  gr2sc = D1(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', 2*norb-ntot, ntot, 2*norb, z1, d1st(2*norb*ntot+1), &
       & 2*norb, scr1, 2*norb, z0, gr2sc, 2*norb-ntot)


  return
end subroutine bldgrad_ghf



subroutine bldgrad_uhf (norb, nup, ndn)

! +----------------------------------------------------------------+
! |                                                                |
! | bldgrad_uhf  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | ( UHF version of bldgrad_ghf. )                                |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons
  !   ndn  - number of spin-dn electrons

  integer, intent(in) :: norb, nup, ndn


  ! other variables

  integer :: k, nosq


  nosq = norb*norb

  ! compute:  scr1 = rho . D1(:,1:N)*

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nup, norb, z1, rho(1,1), norb, &
       & d1st, norb, z0, scr1(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, ndn, norb, z1, rho(1,norb+1), norb, &
       & d1st(nosq+1), norb, z0, scr1(1,nup+1), norb)
  end if

  ! compute:  gr1sc = D1(:,N+1:end)^T . scr1

  if ( nup > 0 ) then
  call zgemm ('c', 'n', norb-nup, nup, norb, z1, d1st(norb*nup+1), &
       & norb, scr1(1,1), norb, z0, gr1sc(1), norb-nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('c', 'n', norb-ndn, ndn, norb, z1, d1st(nosq+norb*ndn+1), &
       & norb, scr1(1,nup+1), norb, z0, gr1sc(1+(norb-nup)*nup), norb-ndn)
  end if


  ! build rho-I; use rho as scratch space

  do k = 1, norb
    rho(k,k) = rho(k,k) - z1
  end do

  do k = 1, norb
    rho(k,k+norb) = rho(k,k+norb) - z1
  end do

  ! compute:  scr2 = (T + Gam) . scr1

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nup, norb, z1, tpgam(1,1), norb, &
       & scr1(1,1), norb, z0, scr2(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, ndn, norb, z1, tpgam(1,norb+1), norb, &
       & scr1(1,nup+1), norb, z0, scr2(1,nup+1), norb)
  end if

  ! compute:  scr1 = (I - rho) . scr2

  if ( nup > 0 ) then
  call zgemm ('n', 'n', norb, nup, norb, -z1, rho(1,1), norb, &
       & scr2(1,1), norb, z0, scr1(1,1), norb)
  end if

  if ( ndn > 0 ) then
  call zgemm ('n', 'n', norb, ndn, norb, -z1, rho(1,norb+1), norb, &
       & scr2(1,nup+1), norb, z0, scr1(1,nup+1), norb)
  end if

  ! compute:  gr2sc = D1(:,N+1:end)^T . scr1

  if ( nup > 0 ) then
  call zgemm ('c', 'n', norb-nup, nup, norb, z1, d1st(norb*nup+1), &
       & norb, scr1(1,1), norb, z0, gr2sc(1), norb-nup)
  end if

  if ( ndn > 0 ) then
  call zgemm ('c', 'n', norb-ndn, ndn, norb, z1, d1st(nosq+norb*ndn+1), &
       & norb, scr1(1,nup+1), norb, z0, gr2sc(1+(norb-nup)*nup), norb-ndn)
  end if


  return
end subroutine bldgrad_uhf



subroutine bldgrad_rhf (norb, nup)

! +----------------------------------------------------------------+
! |                                                                |
! | bldgrad_rhf  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | ( RHF version of bldgrad_ghf. )                                |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   norb - number of orbitals
  !   nup  - number of spin-up electrons

  integer, intent(in) :: norb, nup


  ! other variables

  integer :: k


  ! compute:  scr1 = rho . D1(:,1:N)*

  call zgemm ('n', 'n', norb, nup, norb, z1, rho, norb, &
       & d1st, norb, z0, scr1, norb)

  ! compute:  gr1sc = D1(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', norb-nup, nup, norb, z1, d1st(norb*nup+1), &
       & norb, scr1, norb, z0, gr1sc, norb-nup)


  ! build rho-I; use rho as scratch space

  do k = 1, norb
    rho(k,k) = rho(k,k) - z1
  end do

  ! compute:  scr2 = (T + Gam) . scr1

  call zgemm ('n', 'n', norb, nup, norb, z1, tpgam, norb, &
       & scr1, norb, z0, scr2, norb)

  ! compute:  scr1 = (I - rho) . scr2

  call zgemm ('n', 'n', norb, nup, norb, -z1, rho, norb, &
       & scr2, norb, z0, scr1, norb)

  ! compute:  gr2sc = D1(:,N+1:end)^T . scr1

  call zgemm ('c', 'n', norb-nup, nup, norb, z1, d1st(norb*nup+1), &
       & norb, scr1, norb, z0, gr2sc, norb-nup)


  return
end subroutine bldgrad_rhf



subroutine uhf_to_ghf (nbs, nup, ndn, dst)

! +----------------------------------------------------------------+
! |                                                                |
! | uhf_to_ghf  --  CAJH, 06.2013                                  |
! |                                                                |
! |                                                                |
! | Transform a matrix of orbital coefficients from UHF storage    |
! | to GHF storage. We assume that the allocated space in dst is   |
! | sufficient to accomodate the full matrix.                      |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   nbs - dimension of spin-blocks
  !   nup - number of spin-up electrons
  !   ndn - number of spin-dn electrons

  integer, intent(in) :: nbs, nup, ndn

  ! input / output variables

  complex(kind=dp), dimension(4*nbs*nbs), intent(inout) :: dst

  ! other variables

  integer :: nbssq, ntot, iz1, iz2
  integer :: k, ik_up, ik_dn, ikx


  nbssq = nbs*nbs
  ntot = nup + ndn

  ! first swap down-occupied with up-virtual

  dst(2*nbssq+1:2*nbssq+ndn*nbs) = dst(nbssq+1:nbssq+ndn*nbs)

    iz1 = nbssq-nbs;  iz2 = nbssq-nbs + ndn*nbs

  do k = nbs, nup+1, -1
    dst(iz2+1:iz2+nbs) = dst(iz1+1:iz1+nbs)

    iz1 = iz1 - nbs
    iz2 = iz2 - nbs
  end do

  dst(nbs*nup+1:nbs*ntot) = dst(2*nbssq+1:2*nbssq+ndn*nbs)

  ! now loop over all vectors and fill with zeros

  do k = 2*nbs, 1, -1
    ik_up = (k-1)*2*nbs
    ik_dn = (k-1)*2*nbs + nbs
    ikx   = (k-1)*nbs

    if ( k <= nup .or. ( k > nup+ndn .and. k <= nbs+ndn ) ) then
      dst(ik_up+1:ik_up+nbs)  = dst(ikx+1:ikx+nbs)
      dst(ik_dn+1:ik_dn+nbs)  = z0
    else
      dst(ik_up+1:ik_up+nbs)  = z0
      dst(ik_dn+1:ik_dn+nbs)  = dst(ikx+1:ikx+nbs)
    end if
  end do


  return
end subroutine uhf_to_ghf


end module hamilt


