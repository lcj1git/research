

program phfmol

  use constants
  use i2sint
  use util
  use spatpr
  use iguess
  use purcart
  use phfscf
  ! use ifport

! +----------------------------------------------------------------+
! |                                                                |
! | phfmol  --  CAJH, 01.2013                                      |
! |                                                                |
! |                                                                |
! | Program to perform calculations based on projected Hartree-    |
! | Fock (HF) or Hartree-Fock-Bogoliubov (HFB) configurations on   |
! | molecular systems. The following methods are supported:        |
! |                                                                |
! | - projected HF / HFB                                           |
! |     ( equiv to VAMP, Variation After Mean-field Projection )   |
! |                                                                |
! | - excited VAMP                                                 |
! |                                                                |
! | - FED (few-determinant) VAMP                                   |
! |                                                                |
! | - excited FED VAMP                                             |
! |                                                                |
! | - resonating HF / HFB                                          |
! |                                                                |
! |                                                                |
! | DISCLAIMER : phfmol currently knows nothing about molecular    |
! |   structure. It needs to be provided with one- and two-        |
! |   electron integrals (and even the nuclear repulsion energy).  |
! |   This program concerns only with PHF / PHFB optimizations.    |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | Some useful references discussing the methods used in this     |
! | program are the following:                                     |
! |                                                                |
! | - K. W. Schmid, 'On the use of general symmetry-projected      |
! |     Hartree-Fock-Bogoliubov configurations in variational      |
! |     approaches to the nuclear many-body problem', Prog. Part.  |
! |     Nucl. Phys. 52, 565 (2004).                                |
! |                                                                |
! | - H. Fukutome, 'Theory of Resonating Quantum Fluctuations in   |
! |     a Fermion System. Resonating Hartree-Fock Approximation',  |
! |     Prog. Theor. Phys. 80, 417 (1988).                         |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | phfmol assumes a molecular Hamiltonian of the form             |
! |                                                                |
! |   H  =  ENR  +  sum_{ik}  h_{ik}  c!_i c_k                     |
! |      +  1/4 * sum_{ijkl}  < ij |v| kl > c!_i c!_j c_l c_k,     |
! |                                                                |
! | where h_{ik} are one-electron integrals (typically the sum of  |
! | a kinetic energy part and a nuclear-electron attraction) and   |
! | < ij | v | kl > are antisymmetrized electron-repulsion ints.   |
! | ENR is a constant nuclear-nuclear repulsion term.              |
! |                                                                |
! | The molecular Hamiltonian is assumed to be written in a non-   |
! | orthogonal basis. One typically uses Gaussian or other atom-   |
! | centered basis functions in molecular calculations. Thus, an   |
! | overlap matrix among basis states is expected by the program.  |
! |                                                                |
! | The following properties are assumed about the molecular       |
! | Hamiltonian (these should be true in most cases):              |
! |                                                                |
! | - It commutes with the electron number operator. This allows   |
! |   calculations with particle-number projection.                |
! |                                                                |
! | - It commutes with the total spin operator (S^2) as well as    |
! |   with the z-projection of spin (S_z). This permits us to      |
! |   use spin projection.                                         |
! |                                                                |
! | - It commutes with certain spatial symmetry operations. The    |
! |   molecule defines a point group. Exact solutions to the       |
! |   Schrodinger equation must therefore have spatial symmetry    |
! |   properties. The use of spatial symmetry projection relies    |
! |   on this fact.                                                |
! |                                                                |
! | - It commutes with the complex conjugation operator. Hence we  |
! |   can perform complex conjugation restoration, that is,        |
! |   diagonalize the Hamiltonian in the basis of |phi>, |phi*>.   |
! |                                                                |
! +----------------------------------------------------------------+

! +----------------------------------------------------------------+
! |                                                                |
! | The following variables and concepts are used throughout:      |
! |                                                                |
! |                                                                |
! | std AO basis - standard (non-orthogonal) Atomic Orbital basis  |
! |                                                                |
! | ort AO basis - orthonormal AO basis, obtained by a             |
! |   'canonical orthogonalization' of the std AO basis            |
! |                                                                |
! |   The transformation matrix from the std AO basis to the ort   |
! |   AO basis is stored in xmat, and referred to as the           |
! |   'transformation matrix' [ = S^(-1/2) ].                      |
! |                                                                |
! |                                                                |
! | nbas - number of AO basis functions                            |
! |                                                                |
! | nbct - number of Cartesian AO basis functions                  |
! |                                                                |
! |   2el-integrals should be provided on the Cartesian AO basis.  |
! |   1el-integrals, on the other hand, are expected in the        |
! |   pure function basis if nbas < nbct. **                       |
! |                                                                |
! | norb - number of orbitals (basis functions in ort AO basis)    |
! |                                                                |
! |   This may be smaller than nbas if the set of AO basis         |
! |   functions is not linearly independent.                       |
! |                                                                |
! | ** Note that if pure functions are used (nbas < nbct), then    |
! |    the transformation matrices between pure and Cartesian      |
! |    functions should be provided in 'purcart.dat'.              |
! |      ( see load_purcart for more details )                     |
! |                                                                |
! +----------------------------------------------------------------+

! +----------------------------------------------------------------+
! |                                                                |
! | description of files required by phfmol :                      |
! |                                                                |
! |                                                                |
! | - All integer quantities should be stored as 4-byte, unless    |
! |   otherwise noted.                                             |
! |                                                                |
! | - Each line denoted by => denotes a Fortran record.            |
! |                                                                |
! |                                                                |
! | fi1s - [unformatted] 1el-integrals                             |
! |                                                                |
! |   =>  nbas                                                     |
! |   =>  overlap matrix (dp real, UT matrix)                      |
! |   =>  core Hamiltonian matrix (dp real, UT matrix)             |
! |                                                                |
! |                                                                |
! | fi2s - [unformatted] 2el-integrals                             |
! |                                                                |
! |   =>  nbct, ni2s, lenrec, isymf **                             |
! |   =>  2el-ints, record 1                                       |
! |   =>  2el-ints, record 2                                       |
! |   =>  2el-ints, record ...                                     |
! |                                                                |
! |   'lenrec' denotes the record length, that is, the (maximum)   |
! |   number of two-electron integrals stored in each record.      |
! |                                                                |
! |   Each record should be stored as                              |
! |     (ii(k), ij(k), ik(k), il(k), int(k), k = 1, lenrec)        |
! |   where ii, ij, ik, il are indices of the integral and 'int'   |
! |   is the actual value of the integral:                         |
! |     (ii ij | ik il)  =  int.                                   |
! |                                                                |
! |   Indices should be stored as 2-byte integers, whereas the     |
! |   integral should be stored as a double precision real number. |
! |                                                                |
! |   Two-el integrals SHOULD BE stored in Mulliken notation.      |
! |   Integral symmetry MUST BE used to reduce the number of       |
! |     integrals stored in the file.                              |
! |   Molecular symmetry SHOULD NOT BE used.                       |
! |                                                                |
! |   ** isymf = 1 lets phfmol know that the integrals are stored  |
! |      with a symmetry factor (see fdet)                         |
! |                                                                |
! |                                                                |
! | fdet - [unformatted] MO coefficients                           |
! |                                                                |
! |   =>  iwfnty, nbas, norb, ndet                                 |
! |   =>  MO coeffs (dp complex, full matrix)                      |
! |                                                                |
! |   The form of the MO matrix depends on iwfnty (RHF, UHF, GHF). |
! |   MO coeffs should be provided wrt std AO basis.               |
! |                                                                |
! |                                                                |
! | frot - [unformatted] rotation matrices                         |
! |                                                                |
! |   =>  name of point group (character, len=4)                   |
! |   =>  nop, nbas                                                |
! |   =>  rotmat(:,:,1) (dp real, full nbas x nbas matrix)         |
! |   =>  rotmat(:,:,2)                                            |
! |   =>  ...                                                      |
! |                                                                |
! |   Please check the spatpr module for more details.             |
! |                                                                |
! +----------------------------------------------------------------+

  implicit none

  include 'mpif.h'
  include 'omp_lib.h'


  ! --  CALCULATION PARAMETERS  -- (start) -------------------------

  ! +--------------------------+
  ! |  Hamiltonian parameters  |
  ! +--------------------------+

  !  nbas   - number of spatial symmetry basis functions
  !  ni2s   - number of 2-el integrals

  !  NOTE: nbct is NOT used; we assume integrals are provided in the
  !        desired format
  !  xx nbct   - number of Cartesian basis functions
  !           (used for 2el-integral contraction)

  !   ** Note that if pure functions are used, the transformation
  !      matrices x2pure and x2cart should be provided on file
  !      purcart.dat.
  !        ( see load_purcart for more details )

  !  NOTE: enr is read from fints
  !  xx enr    - nuclear repulsion energy

  !  NOTE: 1-, 2-el integrals are read from fints (below)
  !  xx fi1s   - name of file with 1el-integrals
  !  xx fi2s   - name of file with 2el-integrals

  !  fints  - FCIdump file with integrals

  !  fdet   - name of file with guess determinant to read

  !  frot   - name of file with rotation matrices for spatial
  !           symmetry projection

  ! +-----------------------+
  ! |  number of electrons  |
  ! +-----------------------+

  !  nup    - number of spin-up electrons

  !  ndn    - number of spin-dn electrons

  !   NOTE:  nup, ndn  set the value of m_s (Sz-projection) to use.
  !          The determinants constructed will have exactly nup or
  !          ndn electrons for UHF-type ones, or approximately that
  !          for GHF determinants.

  ! +-------------------------------+
  ! |  method selection parameters  |
  ! +-------------------------------+

  !  iwfnty - type of wavefunction to use **
  !             = defaults to 3

  !             iwfnty  =  1,  RHF type determinant
  !                     =  2,  UHF type determinant
  !                     =  3,  GHF type determinant

  !   ** set by variable wfntyp in input file
  !        wfntyp  =  'rhf'  =>   iwfnty = 1
  !                =  'uhf'  =>          = 2
  !                =  'ghf'  =>          = 3

  !  imethd - PHF method to use **
  !             = defaults to 1

  !             imethd  =  1,  HF
  !                     =  2,  spin-projected HF
  !                     =  3,  PG-projected HF
  !                     =  4,  PG + spin-projected HF
  !                     =  5,  complex-projected HF
  !                     =  6,  complex + spin-projected HF
  !                     =  7,  complex + PG-projected HF
  !                     =  8,  complex + PG + spin-projected HF

  !  Here, PG = point group (spatial symmetry projection).

  !   ** set by variable method in input file
  !        method  =  'hf'     =>   imethd = 1
  !                =  'shf'    =>          = 2
  !                =  'pghf'   =>          = 3
  !                =  'pgshf'  =>          = 4
  !                =  'khf'    =>          = 5
  !                =  'kshf'   =>          = 6
  !                =  'kpghf'  =>          = 7
  !                =  'kpgshf' =>          = 8

  ! +---------------------------+
  ! |  optimization parameters  |
  ! +---------------------------+

  !  iprint - printing level
  !             = defaults to 1

  !             iprint  =  0, print energy and gradient norm per iter
  !                     =  1, also print determinants to be optimized
  !                           before and after the optimization
  !                     =  2, also print Hamiltonian and overlap
  !                           eigenvalues at each function eval
  !                     =  3, also print Hamiltonian and overlap
  !                           matrix at each function eval

  !  iopt   - type of optimization to perform **
  !             = defaults to 1

  !             iopt  =  1,  FED-type (only last det is state is optimized) 
  !                   =  2,  RES-type (all dets in state are optimized)

  !   ** set by variable opttyp in input file
  !        opttyp  =  'fed'   =>   iopt = 1
  !                =  'res'   =>        = 2

  !  igstyp - type of initial guess to use
  !             = defaults to 1

  !             igstyp  =  0,  read from file [ fdet ]
  !                     =  1,  diagonalize core Hamiltonian

  !  igsnbf - number of determinant to load from [ fdet ]
  !             ( only useful when igstyp = 0 )
  !             = defaults to 0 **

  !   ** igsnbf = 0 will load the last det in [ fdet ]

  !  igsmix - mixing of initial guess orbital coefficients
  !             = defaults to 0

  !             igsmix  =  0,  no mixing
  !                     =  1,  very mild mixing
  !                     =  2,  mild mixing
  !                     =  3,  normal mixing
  !                     =  4,  aggressive mixing
  !                     =  5,  very aggressive mixing

  !  igsnmx - number of occ-virtual pairs of orbitals (of each spin)
  !             to mix in igsmix ( only useful when igsmix > 0 )
  !             = defaults to 0 **

  !   ** igsnmx = 0 mixes all the orbitals

  !  maxit  - maximum number of SCF iterations
  !             = defaults to 512

  !  mxfunc - maximum number of energy + gradient evaluations
  !             = defaults to 4*maxit

  !  nbuf   - size of LBFGS buffer (number of vectors to use)
  !             = defaults to 20

  !  gtol   - convergence criterion (norm of the global gradient)
  !             = defaults to 1.0e-4

  ! +-------------------------+
  ! |  projection parameters  |
  ! +-------------------------+

  !  imult  - multiplicity to recover by spin projection **
  !             = defaults to 2*(nup-ndn) + 1

  !   **  s  =    0   =>  imult  =  1 
  !          =  1/2              =  2
  !          =    1              =  3

  !  ngrda  - number of grid points for alpha integration
  !             = defaults to 11 for SGHF
  !                 ( set to 1 otherwise )

  !  ngrdb  - number of grid points for beta  integration
  !             = defaults to 20 for SGHF / SUHF
  !                 ( set to 1 otherwise )

  !  ngrdg  - number of grid points for gamma integration
  !             = defaults to 11 for SGHF
  !                 ( set to 1 otherwise )

  !  igrda  - integration quadrature to use for alpha
  !             = defaults to 1

  !  igrdb  - integration quadrature to use for beta
  !             = defaults to 2

  !  igrdg  - integration quadrature to use for gamma
  !             = defaults to 1

  !             igrdX  =  1,  use trapezoidal grid
  !                    =  2,  use Gauss-Legendre quadrature

  !  sptirr - irrep to recover in spatial symmetry restoration

  !  sptgr  - point group to use in spatial symmetry restoration

  !   ** currently supported combinations sptgr / sptirr:

  !     sptgr  =  'cs  ',    sptirr  =  1,  A'
  !                                  =  2,  A''

  !     sptgr  =  'c2  ',    sptirr  =  1,  A
  !                                  =  2,  B

  !     sptgr  =  'c2h ',    sptirr  =  1,  A_g
  !                                  =  2,  B_g
  !                                  =  3,  A_u
  !                                  =  4,  B_u

  !     sptgr  =  'c2v ',    sptirr  =  1,  A_1
  !                                  =  2,  A_2
  !                                  =  3,  B_1
  !                                  =  4,  B_2

  !     sptgr  =  'c4v ',    sptirr  =  1,  A_1
  !                                  =  2,  A_2
  !                                  =  3,  B_1
  !                                  =  4,  B_2
  !                                  =  5,  E

  !     sptgr  =  'd2  ',    sptirr  =  1,  A
  !                                  =  2,  B_1
  !                                  =  3,  B_2
  !                                  =  4,  B_3

  !     sptgr  =  'd2h ',    sptirr  =  1,  A_g
  !                                  =  2,  B_1g
  !                                  =  3,  B_2g
  !                                  =  4,  B_3g
  !                                  =  5,  A_u
  !                                  =  6,  B_1u
  !                                  =  7,  B_2u
  !                                  =  8,  B_3u

  !     sptgr  =  'd4h ',    sptirr  =  1,  A_1g
  !                                  =  2,  A_2g
  !                                  =  3,  B_1g
  !                                  =  4,  B_2g
  !                                  =  5,  E_g
  !                                  =  6,  A_1u
  !                                  =  7,  A_2u
  !                                  =  8,  B_1u
  !                                  =  9,  B_2u
  !                                  = 10,  E_u

  ! +-----------------+
  ! |  state control  |
  ! +-----------------+

  !  lread  - whether to read old wfn from file 'phfrun.dat'
  !             = defaults to F

  !  lsave  - whether to save optimized wfn in file 'phfrun.dat'
  !             = defaults to T

  !  lnewdt - whether to add a determinant to the expansion
  !             = defaults to F

  !  lnewst - whether to treat a new determinant as a new state
  !           for EXCITED methods
  !             = defaults to F


  integer :: nbas, nbct
  integer :: nup, ndn
  integer :: iwfnty, imethd
  integer :: iprint, iopt, igstyp, igsnbf, igsmix, igsnmx
  integer :: maxit, mxfunc, nbuf
  integer :: ngrda, ngrdb, ngrdg, imult
  integer :: igrda, igrdb, igrdg, sptirr
  logical :: lread, lsave, lnewdt, lnewst

  character(len=4) :: sptgr
  character(len=16) :: fints, fdet, frot
  ! character(len=16) :: fi1s, fi2s

  real(kind=dp) :: enr, gtol


  ! --  CALCULATION PARAMETERS  -- (end) ---------------------------


  ! input, output, data files

  character(len=*) :: inpfile, outfile, datfile
  parameter ( inpfile = 'phfrun.inp', outfile = 'phfrun.out' )
  parameter ( datfile = 'phfrun.dat' )

  integer :: iinp, iout, idat, idat2 = 20
  parameter ( iinp = 11, iout = 12, idat = 13 )

  logical :: ftest


  ! some parameters

  !   mxi2s - maximum number of 2el-ints allowed

  integer :: mxi2s  ! ( 512 MB )
  parameter ( mxi2s = 33554432 )

  !   mxbas - maximum number of basis functions allowed

  integer :: mxbas
  parameter ( mxbas = 500 )

  !   mxgrd - maximum grid size

  integer :: mxgrd
  parameter ( mxgrd = 50 )

  !   mxdet - maximum number of dets allowed
  !   mxstt - maximum number of states allowed

  integer :: mxdet, mxstt
  parameter ( mxdet = 256, mxstt = 20 )


  ! other variables

  real(kind=dp) :: tsrt, tend

  character(len=24) :: systime
  character(len=40) :: hname
  ! character(len=MAX_HOSTNAM_LENGTH+1) :: hname

  ! mpi variables

  integer :: root, procid, nproc, nthrd, ierr
  parameter ( root = 0 )

  ! general variables

  integer :: istatus, nflag
  integer :: ni2s, norb, nosq, xdim, edim
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf, lsf

  integer :: ndspin, ndsptl, ndcplx
  integer :: ndet, ndtt, nstat, mxdtst, ndtopt
  integer, dimension(:), allocatable :: idetst

  real(kind=dp), dimension(:), allocatable :: stval
  real(kind=dp), dimension(:,:), allocatable :: enarr

  complex(kind=dp), dimension(:,:), allocatable :: ovdet, hmdet
  complex(kind=dp), dimension(:,:), allocatable :: ovstat, hmstat, stvec
  complex(kind=dp), dimension(:,:), allocatable :: mdvec
  complex(kind=dp), dimension(:,:), allocatable :: darr

  complex(kind=dp), dimension(:), allocatable :: dmat
  complex(kind=dp), dimension(:,:), allocatable :: smat, hmat, xmat

  type (i2s), dimension(:), allocatable :: i2sv

  integer :: ib


  call ctime (time (), systime)
  call hostnm (hname)
  ! systime = ctime (time ())
  ! istatus = hostnm (hname)


  ! initialize mpi

  call mpi_init (ierr)


  ! determine number of processes, id of each process

  call mpi_comm_rank (mpi_comm_world, procid, ierr)
  call mpi_comm_size (mpi_comm_world, nproc, ierr)

  ! find number of threads
  if ( procid == root ) then
    nthrd = omp_get_max_threads ()
  end if


  if ( procid == root ) then

    ! open output file
    
    inquire (file = outfile, exist = ftest)
    
    if ( ftest ) then
      write (6, *) 'error: File ', outfile, ' exists in current directory.'
      stop
    end if
    
    open (unit = iout, file = outfile, status = 'new', form = 'formatted')

    ! print header information to output file

    write (iout, *)
    write (iout, *) '+------------------------------------------------------+'
    write (iout, *) '|                                                      |'
    write (iout, *) '| phfmol - program to perform projected HF / HFB based |'
    write (iout, *) '|          calculations on molecular systems           |'
    write (iout, *) '|                                                      |'
    write (iout, *) '|                              CAJH, v0.1  --  01.2013 |'
    write (iout, *) '+------------------------------------------------------+'
    write (iout, *)
    write (iout, *)

    ! print date, hostname

    write (iout, *) ' phfmol started on ', trim (hname), ' at ', trim (systime), '.'
    write (iout, '(2X,A,I6,A)') 'phfmol running on up to ', nproc*nthrd, ' thread(s).'

    write (iout, *)
    write (iout, *)

    ! close output file

    close (unit = iout)

  end if


  ! read input file
  ! read 1el-ints, 2el-ints

  if ( procid == root ) then
    call read_input

    ! call read_i1s
    ! call read_i2s
    lsf = .false.
    call read_fints

    norb = nbas
    ! build S, X as identity matrices
    allocate (smat(nbas,nbas), xmat(nbas,nbas))
    smat = d0
    xmat = d0

    do ib = 1, nbas
      smat(ib,ib) = d1
      xmat(ib,ib) = d1
    end do    
  end if


  ! read data file if requested
  ! set defaults otherwise

  if ( procid == root ) then
    nosq = norb*norb
    xdim = 0
    edim = 0

    if ( iwfnty == 1 ) then
      xdim = nosq
      edim = norb
    else if ( iwfnty == 2 ) then
      xdim = 2*nosq
      edim = 2*norb
    else if ( iwfnty == 3 ) then
      xdim = 4*nosq
      edim = 2*norb
    end if

    if ( lread ) then
      call read_data

    else
      ndet   = 1
      ndtt   = ndspin*ndsptl*ndcplx
      nstat  = 1
      mxdtst = 1
      ndtopt = 1

      allocate (darr(xdim,1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'darr', 'phfmol')

      allocate (enarr(edim,1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'enarr', 'phfmol')

      allocate (idetst(1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'idetst', 'phfmol')

      allocate (ovdet(ndtt,ndtt), &
              & hmdet(ndtt,ndtt), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'ovdet, hmdet', 'phfmol')

      allocate (ovstat(1,1), hmstat(1,1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'ovstat, hmstat', 'phfmol')

      allocate (stvec(1,1), stval(1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'stvec, stval', 'phfmol')

      allocate (mdvec(ndtt,1), stat=istatus)
      if ( istatus /= 0 ) call error_alloc (1, 'mdvec', 'phfmol')

      idetst(1) = 1

      ovstat(1,1) = z1
      hmstat(1,1) = z1
      stvec(1,1)  = z1
      stval(1)    = d1

      ovdet(1:ndtt,1:ndtt) = z0
      hmdet(1:ndtt,1:ndtt) = z0

      mdvec(1:ndtt,1) = z0
    end if
  end if


  ! prepare new determinant if needed
  ! finish preparing darr

  if ( procid == root ) then
    if ( ( .not. lread ) .or. ( lread .and. lnewdt ) ) then
      call build_new_dmat
    end if

    if ( lread .and. lnewdt ) then
      darr(1:xdim,ndet) = dmat(1:xdim)
      enarr(1:edim,ndet) = d0
    end if

    if ( .not. lread ) then
      darr(1:xdim,1) = dmat(1:xdim)
      enarr(1:edim,1) = d0
    end if
  end if


  ! broadcast calculation parameters

  call phfmol_send_par


  ! load pure to Cartesian transformation matrices

  ! if ( nbct > nbas )  call load_purcart


  ! perform optimization

    if ( procid == root ) then
      open (unit = iout, file = outfile, status = 'old', &
          & form = 'formatted', position = 'append')

      write (iout, *)
      write (iout, '(2X,A)') 'SCF start.'
      write (iout, *)

      close (unit = iout)
    end if

    tsrt = mpi_wtime ()
    ! call cpu_time (tsrt)

  call scfhf_lbfgs (root, procid, nproc, iout, outfile, iprint, imethd, &
       & iwfnty, norb, nbas, nbct, nup, ndn, xmat, hmat, smat, i2sv, ni2s, &
       & lsf, enr, imult, ngrda, ngrdb, ngrdg, igrda, igrdb, igrdg, sptgr, &
       & sptirr, frot, ndspin, ndsptl, ndcplx, darr, enarr, xdim, edim, &
       & ndet, ndtt, ovdet, hmdet, nstat, idetst, mxdtst, ovstat, hmstat, &
       & stvec, stval, mdvec, iopt, ndtopt, maxit, mxfunc, nbuf, gtol, nflag)

    tend = mpi_wtime ()
    ! call cpu_time (tend)

    if ( procid == root ) then
      open (unit = iout, file = outfile, status = 'old', &
          & form = 'formatted', position = 'append')

      write (iout, *)

      select case (nflag)
        case (0)
          write (iout, 21)
        case (1)
          write (iout, 22)
        case (2)
          write (iout, 23)
        case (3)
          write (iout, 24)
      end select

      21 format (2X, 'SCF done. Successful completion.')
      22 format (2X, 'SCF done. Optimization failed.')
      23 format (2X, 'SCF done. Number of iterations exhausted. Re-run.')
      24 format (2X, 'SCF done. Number of fcn evals exhausted. Re-run.')

      write (iout, '(2X,A,F12.2,A)') 'Elapsed time: ', tend-tsrt, ' sec.'
      write (iout, *)

      close (unit = iout)
    end if


  ! save data file

  if ( procid == root ) then
    if ( lsave )  call save_data
    call save_det
  end if


  ! deallocate memory

  !$omp  parallel

  if ( nbct > nbas )  call shutdown_purcart

  !$omp  end parallel

  call shutdown_phfmol


  ! finalize program

  call mpi_finalize (ierr)

  stop

contains


subroutine read_input

! +----------------------------------------------------------------+
! |                                                                |
! | read_input  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Read input file 'phfrun.inp' for phfmol calculations and       |
! | interpret it. Verify the consistency of all parameters.        |
! | Lastly, the full list of parameters is printed to the output   |
! | file 'phfrun.out'.                                             |
! |                                                                |
! | Standard Fortran namelists are used to retrieve information    |
! | from the input file.                                           |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  character(len=10) :: wfntyp, method, opttyp


  ! namelist declarations

  ! namelist /file_data/ fi1s, fi2s, fdet, frot
  ! namelist /ham_data/ nbas, nbct, enr
  namelist /ham_data/ nbas, ni2s, fints, fdet, frot

  namelist /nel_data/ nup, ndn

  namelist /mtd_data/ wfntyp, method

  namelist /opt_data/ iprint, opttyp, igstyp, igsnbf, igsmix, &
                    & igsnmx, maxit, mxfunc, nbuf, gtol

  namelist /prj_data/ imult, ngrda, ngrdb, ngrdg, &
                    & igrda, igrdb, igrdg, sptgr, sptirr

  namelist /stt_data/ lread, lsave, lnewdt, lnewst
  

  ! default values in ham_data
  !   ( Hamiltonian parameters )

  nbas = 0
  ni2s = 0
  ! nbct = 0
  ! enr  = 0.0e0_dp

  ! default values in nel_data

  nup = 0
  ndn = 0

  ! default values in mtd_data
  !   ( method selection )

  wfntyp = 'ghf'
  method = 'hf'

  ! default values in opt_data
  !   ( optimization parameters )

  iprint = 1
  opttyp = 'fed'
  igstyp = -1
  igsnbf = 0
  igsmix = 0
  igsnmx = 0
  maxit  = 512
  mxfunc = -1
  nbuf   = 20
  gtol   = 1.0e-4_dp

  ! default values in prj_data
  !   ( projection parameters )

  imult = -1
  igrda = 1
  igrdb = 2
  igrdg = 1
  ngrda = -1
  ngrdb = -1
  ngrdg = -1
  sptgr = ' '
  sptirr = 0

  ! default values in stt_data
  !   ( state control )

  lread  = .false.
  lsave  = .true.
  lnewdt = .false.
  lnewst = .false.

  ! default values in file_data

  ! fi1s = ' '
  ! fi2s = ' '
  fints = ''
  fdet = ''
  frot = ''


  ! open input file

  inquire (file = inpfile, exist = ftest)

  if ( .not. ftest ) then
    write (6, *) 'error: Input file ', inpfile, ' does not exist.'
    stop
  end if

  open (unit = iinp, file = inpfile, status = 'old', &
      & form = 'formatted', delim = 'apostrophe')

  read (iinp, ham_data)  ! read ham_data
  nbct = nbas

  read (iinp, nel_data)  ! read nel_data

  read (iinp, mtd_data)  ! read mtd_data

    select case (wfntyp)
      case ('rhf')
        iwfnty = 1
      case ('uhf')
        iwfnty = 2
      case ('ghf')
        iwfnty = 3
      case default
        iwfnty = 0
    end select

    select case (method)
      case ('hf')
        imethd = 1
      case ('shf')
        imethd = 2
      case ('pghf')
        imethd = 3
      case ('pgshf')
        imethd = 4
      case ('khf')
        imethd = 5
      case ('kshf')
        imethd = 6
      case ('kpghf')
        imethd = 7
      case ('kpgshf')
        imethd = 8
      case default
        imethd = 0
    end select

  read (iinp, opt_data)  ! read opt_data

    select case (opttyp)
      case ('fed')
        iopt = 1
      case ('res')
        iopt = 2
      case default
        iopt = 0
    end select

    ! set default values  (igstyp, mxfunc)

    if ( igstyp == -1 ) then
      igstyp = 1
    end if

    if ( mxfunc == -1 ) then
      mxfunc = 4*maxit
    end if

  read (iinp, prj_data)  ! read prj_data

    ! some logical variables

    lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
    lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
    lspin = mod(imethd-1,2) == 1

    lsuhf = .false.
    lsghf = .false.

    if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
    if ( lspin .and. iwfnty == 3 )  lsghf = .true.

    ! set default values for spin projection

    if ( imult == -1 .and. lspin ) then
      imult = 2*(nup-ndn) + 1
    else if ( .not. lspin ) then
      imult = 1
    end if

    if ( ngrda == -1 .and. lsghf ) then
      ngrda = 11
    else if ( .not. lsghf ) then
      ngrda = 1
    end if

    if ( ngrdb == -1 .and. lspin ) then
      ngrdb = 20
    else if ( .not. lspin ) then
      ngrdb = 1
    end if

    if ( ngrdg == -1 .and. lsghf ) then
      ngrdg = 11
    else if ( .not. lsghf ) then
      ngrdg = 1
    end if

  read (iinp, stt_data)  ! read state control data

  ! read (iinp, file_data) ! read file description data


  ! close input file

  close (unit = iinp)


  ! verify input variables

  call check_input


  ! print input selection to output file

  open (unit = iout, file = outfile, status = 'old', &
      & form = 'formatted', position = 'append')


  write (iout, *) ' Hamiltonian parameters :'
  write (iout, *)

  write (iout, 10) 'nbas  ', nbas
  write (iout, 10) 'ni2s  ', ni2s
  ! write (iout, 10) 'nbct  ', nbct
  ! write (iout, 11) 'enr   ', enr

  write (iout, 15) 'fints ', fints
  if ( igstyp == 0 ) then
    if ( fdet == '' ) then
      write (6, *) 'error: igstyp = 0 but fdet not provided'
      stop
    end if
    write (iout, 15) 'fdet  ', fdet
  end if

  if ( lsptl ) then
    if ( frot == '' ) then
      write (6, *) 'error: lsptl = T but frot not provided'
      stop
    end if
    write (iout, 15) 'frot  ', frot
  end if

  write (iout, *)
  write (iout, *) ' Number of electrons :'
  write (iout, *)

  write (iout, 10) 'nup   ', nup
  write (iout, 10) 'ndn   ', ndn

  write (iout, *)
  write (iout, *) ' Method selection :'
  write (iout, *)

  select case (iwfnty)
    case (1)
      write (iout, 120) iwfnty, 'RHF'
    case (2)
      write (iout, 120) iwfnty, 'UHF'
    case (3)
      write (iout, 120) iwfnty, 'GHF'
  end select

  select case (imethd)
    case (1)
      write (iout, 130) imethd
    case (2)
      write (iout, 131) imethd, 'spin'
    case (3)
      write (iout, 131) imethd, 'PG'
    case (4)
      write (iout, 131) imethd, 'PG + spin'
    case (5)
      write (iout, 131) imethd, 'complex'
    case (6)
      write (iout, 131) imethd, 'complex + spin'
    case (7)
      write (iout, 131) imethd, 'complex + PG'
    case (8)
      write (iout, 131) imethd, 'complex + PG + spin'
  end select

  write (iout, *)
  write (iout, *) ' Optimization parameters :'
  write (iout, *)

  write (iout, 10) 'iprint', iprint
  write (iout, 10) 'igstyp', igstyp

  if ( igstyp == 0 ) then
    write (iout, 10) 'igsnbf', igsnbf
  end if

  write (iout, 10) 'igsmix', igsmix

  if ( igsmix > 0 ) then
    write (iout, 10) 'igsnmx', igsnmx
  end if

  write (iout, 10) 'maxit ', maxit
  write (iout, 10) 'mxfunc', mxfunc
  write (iout, 10) 'nbuf  ', nbuf
  write (iout, 13) 'gtol  ', gtol

  select case (iopt)
    case (1)
      write (iout, 140) iopt, 'FED'
    case (2)
      write (iout, 140) iopt, 'RES'
  end select

  write (iout, *)
  write (iout, *) ' Projection parameters :'
  write (iout, *)

  if ( lspin ) then
    write (iout, 10) 'imult ', imult

    if ( lsuhf ) then
      write (iout, 10) 'ngrdb ', ngrdb
      write (iout, 10) 'igrdb ', igrdb

    else if ( lsghf ) then
      write (iout, 10) 'ngrda ', ngrda
      write (iout, 10) 'ngrdb ', ngrdb
      write (iout, 10) 'ngrdg ', ngrdg
      write (iout, 10) 'igrda ', igrda
      write (iout, 10) 'igrdb ', igrdb
      write (iout, 10) 'igrdg ', igrdg
    end if
  end if

  if ( lsptl ) then
    write (iout, 14) 'sptgr ', sptgr
    write (iout, 10) 'sptirr', sptirr
  end if

  write (iout, *)
  write (iout, *) ' State control :'
  write (iout, *)

  write (iout, 12) 'lread ', lread
  write (iout, 12) 'lsave ', lsave
  write (iout, 12) 'lnewdt', lnewdt
  write (iout, 12) 'lnewst', lnewst

  ! write (iout, *)
  ! write (iout, *) ' File names :'
  ! write (iout, *)

  ! write (iout, 15) 'fi1s  ', fi1s
  ! write (iout, 15) 'fi2s  ', fi2s

  ! if ( igstyp == 0 ) then
  !   write (iout, 15) 'fdet  ', fdet
  ! end if

  ! if ( lsptl ) then
  !   write (iout, 15) 'frot  ', frot
  ! end if


  write (iout, *)

  close (unit = iout)


  ! format statements

  10  format (4X, A, ' = ', I12)
  11  format (4X, A, ' = ', F12.4)
  12  format (4X, A, ' = ', L12)
  13  format (4X, A, ' = ', 1P, E12.3)
  14  format (4X, A, ' = ', 8X, A4)
  15  format (4X, A, ' = ', A16)

  120 format (4X, 'iwfnty = ', I12, '  --  ', A, ' determinant')
  130 format (4X, 'imethd = ', I12, '  --  ', 'HF')
  131 format (4X, 'imethd = ', I12, '  --  ', A, '-projected HF')
  140 format (4X, 'iopt   = ', I12, '  --  ', A, ' type optimization')


  return
end subroutine read_input



subroutine check_input

! +----------------------------------------------------------------+
! |                                                                |
! | check_input  --  CAJH, 01.2013                                 |
! |                                                                |
! |                                                                |
! | Perform consistency checks among input parameters.             |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: ntot


  ! check that some parameters have been recognized

  if ( iwfnty == 0 ) then
    write (6, *) 'error: The specified wfn type is not supported.'
    stop
  end if

  if ( imethd == 0 ) then
    write (6, *) 'error: The specified method is not supported.'
    stop
  end if

  if ( iopt == 0 ) then
    write (6, *) 'error: The specified opt type is not supported.'
    stop
  end if


  ! +--------------------------+
  ! |  Hamiltonian parameters  |
  ! +--------------------------+

  ! check number of basis functions

  if ( nbas <= 0 ) then
    write (6, *) 'error: Number of basis must be positive.'
    stop
  else if ( nbas > mxbas ) then
    write (6, *) 'error: Number of basis functions is too large.'
    stop
  end if


  ! check that nbas <= nbct

  if ( nbct < nbas ) then
    write (6, *) 'error: nbct must be >= to nbas.'
    stop
  end if

  
  ! check that number of 2el-ints is within bounds

  if ( ni2s <= 0 ) then
    write (6, *) 'error: ni2s <= 0 in fi2s file.'
    stop

  else if ( ni2s > mxi2s ) then
    write (6, *) 'error: ni2s is too large.'
    stop
  end if


  ! +-----------------------+
  ! |  number of electrons  |
  ! +-----------------------+

  ! ntot - total number of electrons

  ntot = nup + ndn

  ! check number of electrons

  if ( nup < 0 .or. ndn < 0 ) then
    write (6, *) 'error: Number of electrons cannot be negative.'
    stop

  else if ( nup == 0 .and. ndn == 0 ) then
    write (6, *) 'error: There are no electons in the system.'
    stop

  else if ( nup > nbas .or. ndn > nbas ) then
    write (6, *) 'error: Too many up or dn electrons.'
    stop
  end if


  ! +--------------------+
  ! |  method selection  |
  ! +--------------------+

  ! methods involving spin projection require UHF or GHF

  if ( mod(imethd-1,2) == 1 ) then
    if ( iwfnty == 1 ) then
      write (6, *) 'error: No spin projection for RHF wfns.'
      stop
    end if
  end if

  ! RHF requires equal number of up and dn electrons

  if ( iwfnty == 1 .and. nup /= ndn ) then
    write (6, *) 'error: RHF wfn but nup /= ndn.'
    stop
  end if


  ! +---------------------------+
  ! |  optimization parameters  |
  ! +---------------------------+

  ! printing level

  if ( iprint < 0 .or. iprint > 3 ) then
    iprint = 0
  end if

  ! initial guess selection

  !if ( igstyp < 0 .or. igstyp > 1 ) then
!gomez59
  if ( igstyp < 0 .or. igstyp > 2 ) then
    write (6, *) 'error: Selected igstyp is not available.'
    stop
  end if

  if ( igstyp == 0 .and. igsnbf < 0 ) then
    write (6, *) 'error: igsnbf must be non-negative.'
    stop
  end if

  if ( igsmix < 0 .or. igsmix > 5 ) then
    write (6, *) 'error: Selected igsmix is not available.'
    stop
  end if

  if ( igsmix > 0 ) then
    if ( igsnmx < 0 ) then
      write (6, *) 'error: igsnmx must be non-negative.'
      stop
    end if

    if ( iwfnty == 1 ) then
      if ( igsnmx > nup .or. igsnmx > nbas-nup ) then
        write (6, *) 'error: igsnmx is too large.'
        stop
      end if

    else if ( iwfnty == 2 ) then
      if ( igsnmx > nup .or. igsnmx > ndn .or. &
           igsnmx > nbas-nup .or. igsnmx > nbas-ndn ) then
        write (6, *) 'error: igsnmx is too large.'
        stop
      end if

    else if ( iwfnty == 3 ) then
      if ( igsnmx > ntot .or. igsnmx > 2*nbas-ntot ) then
        write (6, *) 'error: igsnmx is too large.'
        stop
      end if
    end if
  end if


  ! adjust other parameters

  if ( maxit < 0 ) then
    write (6, *) 'error: maxit < 0 is not allowed.'
    stop
  end if

  if ( mxfunc < maxit ) then
    mxfunc = maxit
  end if

  if ( nbuf < 1 ) then
    write (6, *) 'error: Size of LBFGS buffer must be >= 1.'
    stop
  else if ( nbuf > 100 ) then
    write (6, *) 'error: Size of LBFGS buffer is too large.'
    stop
  end if

  if ( gtol > 1.0e1_dp ) then
    write (6, *) 'error: Unreasonable selection for gtol.'
    stop
  else if ( gtol < 1.0e-12_dp ) then
    gtol = 1.0e-12_dp
  end if


  ! +-------------------------+
  ! |  projection parameters  |
  ! +-------------------------+

    ! determine ndcplx

    if ( lcplx ) then
      ndcplx = 2
    else
      ndcplx = 1
    end if

  ! check multiplicity

  if ( imult < 0 .or. imult > 2*ntot + 1 ) then
    write (6, *) 'error: Unphysical multiplicity selected.'
    stop
  end if

    ! determine ndspin

    if ( lsghf ) then
      ndspin = imult
    else
      ndspin = 1
    end if

  ! check grid dimensions

  if ( ngrda < 1 .or. ngrda > mxgrd ) then
    write (6, *) 'error: Incorrect alpha integration dimension.'
    stop
  end if

  if ( ngrdb < 1 .or. ngrdb > mxgrd ) then
    write (6, *) 'error: Incorrect beta integration dimension.'
    stop
  end if

  if ( ngrdg < 1 .or. ngrdg > mxgrd ) then
    write (6, *) 'error: Incorrect gamma integration dimension.'
    stop
  end if

  ! check grid type

  if ( ngrda > 1 .and. igrda /= 1 .and. igrda /= 2 ) then
    write (6, *) 'error: Unrecognized quadrature in alpha integration.'
    stop
  end if

  if ( ngrdb > 1 .and. igrdb /= 1 .and. igrdb /= 2 ) then
    write (6, *) 'error: Unrecognized quadrature in alpha integration.'
    stop
  end if

  if ( ngrdg > 1 .and. igrdg /= 1 .and. igrdg /= 2 ) then
    write (6, *) 'error: Unrecognized quadrature in alpha integration.'
    stop
  end if

  ! check spatial projection
  ! determine ndsptl

  if ( lsptl ) then

    select case (sptgr)
      case ('cs  ')
        call chk_irrep_cs  (sptirr, ndsptl)
      case ('c2  ')
        call chk_irrep_c2  (sptirr, ndsptl)
      case ('c2h ')
        call chk_irrep_c2h (sptirr, ndsptl)
      case ('c2v ')
        call chk_irrep_c2v (sptirr, ndsptl)
      case ('c4v ')
        call chk_irrep_c4v (sptirr, ndsptl)
      case ('c6v ')
        call chk_irrep_c6v (sptirr, ndsptl)
      case ('c8v ')
        call chk_irrep_c8v (sptirr, ndsptl)
      case ('c12v')
        call chk_irrep_c12v(sptirr, ndsptl)
      case ('c16v')
        call chk_irrep_c16v(sptirr, ndsptl)
      case ('d2  ')
        call chk_irrep_d2  (sptirr, ndsptl)
      case ('d2h ')
        call chk_irrep_d2h (sptirr, ndsptl)
      case ('d4h ')
        call chk_irrep_d4h (sptirr, ndsptl)
      case ('d6h ')
        call chk_irrep_d6h (sptirr, ndsptl)
      case ('d8h ')
        call chk_irrep_d8h (sptirr, ndsptl)
      case ('d12h')
        call chk_irrep_d12h(sptirr, ndsptl)
      case ('d16h')
        call chk_irrep_d16h(sptirr, ndsptl)

      case default
        write (6, *) 'error: Point group not supported in phfmol.'
        stop
    end select

  else if ( .not. lsptl ) then
    ndsptl = 1
  end if


  ! +-----------------+
  ! |  state control  |
  ! +-----------------+

  if ( .not. lread ) then   ! cannot add new det unless reading
    lnewdt = .false.
    lnewst = .false.
  end if

  if ( lread .and. lnewst ) then   ! if adding new state,
    lnewdt = .true.                ! add new determinant
  end if


  return
end subroutine check_input


subroutine read_fints ()

  implicit none

  integer :: ist, iz
  integer :: ji, jj, jk, jl
  real(kind=dp) :: val

  integer :: norb, nelec, ms2, isym
!  integer, dimension(:), allocatable :: orbsym
  integer, dimension(1000) :: orbsym
  namelist /FCI/ norb, nelec, ms2, orbsym, isym

  inquire (file=fints, exist=ftest)
  if ( .not. ftest) then
    write (6,*) 'failed ftest'
    stop
  end if

!  allocate (orbsym(nbas))
  open (unit=idat, file=fints, status='old', form='formatted')
  read (idat,fci)

  if ( norb /= nbas ) then
    write (6,*) 'inconsistent norb in fints'
    stop
  end if

  ist = 0
  iz  = 1
  enr = d0

  allocate (hmat(nbas,nbas))
  hmat = d0
  allocate (i2sv(ni2s))

  do while ( ist == 0 )
    read (idat,*,iostat=ist) val, ji, jj, jk, jl

    if ( ist == 0 ) then
      if ( ji /= 0 .and. jj /= 0 .and. &
         & jk /= 0 .and. jl /= 0 ) then
        if ( iz > ni2s ) then
          iz = iz+1
          exit
        end if

        i2sv(iz)%ii = int(ji,kind=2)
        i2sv(iz)%ij = int(jj,kind=2)
        i2sv(iz)%ik = int(jk,kind=2)
        i2sv(iz)%il = int(jl,kind=2)
        i2sv(iz)%int = val
        iz = iz+1

      else if ( ji == 0 .and. jj == 0 .and. &
              & jk == 0 .and. jl == 0 ) then
        enr = enr + val

      else if ( jk == 0 .and. jl == 0 ) then
        if ( ji /= jj ) then
          hmat(ji,jj) = hmat(ji,jj) + val
          hmat(jj,ji) = hmat(jj,ji) + val
        else
          hmat(ji,ji) = hmat(ji,ji) + val
        end if
      else
        write (6,*) 'confusing line in file fints'
        stop
      end if
    end if
  end do

  if ( .not. is_iostat_end(ist) ) then
    write (6,*) 'read error in file fints'
    stop
  end if
  if ( iz-1 /= ni2s ) then
    write (6,*) 'number of 2-el integrals read is not ni2s'
    stop
  end if

!  deallocate (orbsym)
  close (unit=idat)

end subroutine read_fints


! subroutine read_i1s

! ! +----------------------------------------------------------------+
! ! |                                                                |
! ! | read_i1s  --  CAJH, 01.2013                                    |
! ! |                                                                |
! ! |                                                                |
! ! | Read one-electron integrals (overlap, core Hamiltonian) from   |
! ! | fi1s file. Working versions of the overlap and core Ham        |
! ! | matrices are prepared.                                         |
! ! |                                                                |
! ! | Additionally, the transformation matrix [ = S^(-1/2) ] is      |
! ! | prepared by diagonalization of the overlap matrix.             |
! ! |                                                                |
! ! +----------------------------------------------------------------+


!   implicit none

!   ! other variables

!   integer :: nbas1, ntt, j, k, icnt, idim
!   integer :: istatus, info

!   real(kind=dp), dimension(:), allocatable :: smtpk, hmtpk
!   real(kind=dp), dimension(:), allocatable :: sval, wrk
!   real(kind=dp), dimension(:,:), allocatable :: svec


!   ! some parameters
!   !   thresh - threshold used for linear dependencies in overlap matrix

!   real(kind=dp) :: thresh
!   parameter ( thresh = 1.0e-6_dp )


!   ! open fi1s file

!   inquire (file = fi1s, exist = ftest)

!   if ( .not. ftest ) then
!     write (6, *) 'error: 1el-ints file ', fi1s, ' does not exist.'
!     stop
!   end if

!   open (unit = idat, file = fi1s, status = 'old', form = 'unformatted')

!   read (idat) nbas1

!   if ( nbas1 /= nbas ) then
!     write (6, *) 'error: Inconsistent nbas in fi1s file.'
!     stop
!   end if

!   ! read overlap and core Hamiltonian matrices

!   ntt = nbas*(nbas+1)/2

!   allocate (smtpk(ntt), hmtpk(ntt), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'smtpk, hmtpk', 'read_i1s')

!   read (idat) (smtpk(k), k = 1, ntt)
!   read (idat) (hmtpk(k), k = 1, ntt)

!   ! close fi1s file

!   close (unit = idat)


!   ! prepare working versions of overlap, core Hamiltonian matrices

!   allocate (smat(nbas,nbas), hmat(nbas,nbas), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'smat, hmat', 'read_i1s')

!   icnt = 1

!   do j = 1, nbas
!     do k = 1, j
!       smat(k,j) = smtpk(icnt)
!       hmat(k,j) = hmtpk(icnt)

!       if ( k /= j ) then
!         smat(j,k) = smtpk(icnt)
!         hmat(j,k) = hmtpk(icnt)
!       end if

!       icnt = icnt + 1
!     end do
!   end do


!   ! +------------------+
!   ! |  construct xmat  |
!   ! +------------------+

!   ! allocate memory for several arrays

!   allocate (svec(nbas,nbas), sval(nbas), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'svec, sval', 'read_i1s')
!   do k = 2, nbas
!     if ( abs(sval(k)) < thresh ) then
!       idim = k-1
!       exit
!     end if
!   end do

!   allocate (wrk(3*nbas), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'wrk', 'read_i1s')


!   ! diagonalize smtpk
!   ! set to negative to recover largest eigenvalues on top

!   smtpk(1:ntt) = -smtpk(1:ntt)

!   call dspev ('v', 'u', nbas, smtpk, sval, svec, nbas, wrk, info)

!   if ( info /= 0 ) then
!     write (6, *) 'error: Problems diagonalizing overlap matrix.'
!     stop
!   end if

!   ! determine norb

!   idim = nbas

!   do k = 2, nbas
!     if ( abs(sval(k)) / abs(sval(1)) < thresh ) then
!       idim = k-1
!       exit
!     end if
!   end do

!   norb = idim

!   ! prepare xmat

!   allocate (xmat(nbas,norb), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'xmat', 'read_i1s')

!   do j = 1, norb
!     do k = 1, nbas
!       xmat(k,j) = svec(k,j) / sqrt(-sval(j))
!     end do
!   end do


!   ! deallocate memory

!   deallocate (smtpk, hmtpk, stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (2, 'smtpk, hmtpk', 'read_i1s')

!   deallocate (svec, sval, stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (2, 'svec, sval', 'read_i1s')

!   deallocate (wrk, stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (2, 'wrk', 'read_i1s')

!   return
! end subroutine read_i1s



! subroutine read_i2s

! ! +----------------------------------------------------------------+
! ! |                                                                |
! ! | read_i2s  --  CAJH, 01.2013                                    |
! ! |                                                                |
! ! |                                                                |
! ! | Read two-electron integrals from fi2s file.                    |
! ! |                                                                |
! ! +----------------------------------------------------------------+


!   implicit none

!   ! other variables

!   integer :: nbct1, isymf
!   integer :: nrec, lenrec
!   integer :: i, z, k, istatus


!   ! open fi2s file

!   inquire (file = fi2s, exist = ftest)

!   if ( .not. ftest ) then
!     write (6, *) 'error: 2el-ints file ', fi2s, ' does not exist.'
!     stop
!   end if

!   open (unit = idat, file = fi2s, status = 'old', form = 'unformatted')

!   read (idat) nbct1, ni2s, lenrec, isymf

!   if ( nbct1 /= nbct ) then
!     write (6, *) 'error: Inconsistent nbct in fi2s file.'
!     stop
!   end if

!   lsf = .false.
!   if ( isymf == 1 )  lsf = .true.

!   ! check that number of 2el-ints is within bounds

!   if ( ni2s <= 0 ) then
!     write (6, *) 'error: ni2s <= 0 in fi2s file.'
!     stop

!   else if ( ni2s > mxi2s ) then
!     write (6, *) 'error: ni2s is too large.'
!     stop
!   end if


!   ! allocate memory for 2el-ints

!   allocate (i2sv(ni2s), stat=istatus)
!   if ( istatus /= 0 ) call error_alloc (1, 'i2sv', 'read_i2s')


!   ! read 2el-ints

!   if ( mod (ni2s,lenrec) == 0 ) then
!     nrec = ni2s/lenrec
!   else
!     nrec = ni2s/lenrec + 1
!   end if

!   i = 1

!   do z = 1, nrec
!     if ( i+lenrec-1 <= ni2s ) then
!       read (idat) (i2sv(k), k = i, i+lenrec-1)
!     else
!       read (idat) (i2sv(k), k = i, ni2s)
!     end if

!     i = i + lenrec
!   end do


!   ! close fi2s file

!   close (unit = idat)


!   return
! end subroutine read_i2s



subroutine read_data

! +----------------------------------------------------------------+
! |                                                                |
! | read_data  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Read all important variables from previous data file.          |
! | Make sure that there are no inconsistencies wrt current        |
! | calculation.                                                   |
! |                                                                |
! | If all tests succeed, then load all previous information and   |
! | prepare some important arrays for current calculation.         |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: j, k
  integer :: nbas1, nbct1, norb1, nup1, ndn1
  integer :: iwfnty1, imethd1
  integer :: imult1, sptirr1
  integer :: ndet1, nstat1, ndtt1, mxdtst1, nflg1
  integer :: istatus
  logical :: lchk

  character(len=4) :: sptgr1

  complex(kind=dp), dimension(:,:), allocatable :: xmat1


  ! open data file

  inquire (file = datfile, exist = ftest)

  if ( .not. ftest ) then
    write (6, *) 'error: Data file ', datfile, ' does not exist.'
    stop
  end if

  open (unit = idat, file = datfile, status = 'old', form = 'unformatted')


  10 format ('error: Parameters in data file do not coincide ', &
          & 'with current ones: ', A, '.')


  ! read several integer scalars

  read (idat) nbas1, nbct1, norb1, nup1, ndn1, iwfnty1, imethd1

  if ( nbas1 /= nbas ) then
    write (6, 10) 'nbas'
    stop
  end if

  if ( nbct1 /= nbct ) then
    write (6, 10) 'nbct'
    stop
  end if

  if ( norb1 /= norb ) then
    write (6, 10) 'norb'
    stop
  end if

  if ( nup1 /= nup .or. ndn1 /= ndn ) then
    write (6, 10) 'nup, ndn'
    stop
  end if

  if ( iwfnty1 /= iwfnty ) then
    write (6, 10) 'iwfnty'
    stop
  end if

  if ( imethd1 /= imethd ) then
    write (6, 10) 'imethd'
    stop
  end if


  ! read spin projection parameters

  if ( mod(imethd-1,2) == 1 ) then
    read (idat) imult1
    read (idat)
    read (idat)

    if ( imult1 /= imult ) then
      write (6, 10) 'imult'
      stop
    end if
  end if


  ! read spatial projection parameters

  if ( mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3 ) then
    read (idat) sptgr1
    read (idat) sptirr1

    if ( sptgr1 /= sptgr ) then
      write (6, 10) 'sptgr'
      stop
    end if

    if ( sptirr1 /= sptirr ) then
      write (6, 10) 'sptirr'
      stop
    end if
  end if


  ! read xmat

  allocate (xmat1(nbas,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'xmat1', 'read_data')

  read (idat) ((xmat1(j,k), j = 1, nbas), k = 1, norb)

  ! test whether X! . S . X = I
  ! if so, use previous xmat; otherwise, abort calculation

  call check_xmat (xmat1, lchk)

  if ( lchk ) then
    xmat(1:nbas,1:norb) = xmat1(1:nbas,1:norb)
  else
    write (6, *) 'error: Test X!.S.X = I failed in read_data.'
    stop
  end if

  deallocate (xmat1, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'xmat1', 'read_data')


  ! read previously stored number of determinants

  read (idat) ndet1, nstat1

  ! fix ndet and nstat

  nstat = 0
  ndet  = 0

  if ( (.not. lnewst) .and. (.not. lnewdt) ) then
    nstat = nstat1
    ndet  = ndet1

  else if ( (.not. lnewst) .and. lnewdt ) then
    nstat = nstat1
    ndet  = ndet1 + 1

  else if ( lnewst .and. lnewdt ) then
    nstat = nstat1 + 1
    ndet  = ndet1 + 1
  end if

  ! decide whether ndet and nstat are acceptable

  if ( ndet > mxdet ) then
    write (6, *) 'error: Number of determinants is too large.'
    stop

  else if ( nstat > mxstt ) then
    write (6, *) 'error: Number of states is too large.'
    stop
  end if

  ! allocate idetst and fill

  allocate (idetst(nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'idetst', 'read_data')

  if ( (.not. lnewst) .and. (.not. lnewdt) ) then
    read (idat) (idetst(k), k = 1, nstat)

  else if ( (.not. lnewst) .and. lnewdt ) then
    read (idat) (idetst(k), k = 1, nstat)
    idetst(nstat) = idetst(nstat) + 1

  else if ( lnewst .and. lnewdt ) then
    read (idat) (idetst(k), k = 1, nstat-1)
    idetst(nstat) = 1
  end if

  ! determine several important variables

  ndtt   = ndet*ndspin*ndsptl*ndcplx
  ndtt1  = ndet1*ndspin*ndsptl*ndcplx
  mxdtst = maxval (idetst(1:nstat))

  if ( iopt == 1 ) then
    ndtopt = 1
  else
    ndtopt = idetst(nstat)
  end if

  if ( (.not. lnewst) .and. lnewdt ) then
    mxdtst1 = maxval (idetst(1:nstat-1))
    mxdtst1 = max (mxdtst1, idetst(nstat)-1)

  else if ( lnewst .and. lnewdt ) then
    mxdtst1 = maxval (idetst(1:nstat-1))

  else
    mxdtst1 = mxdtst
  end if

  ! allocate important arrays

  allocate (ovdet(ndtt,ndtt), hmdet(ndtt,ndtt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'ovdet, hmdet', 'read_data')

  allocate (ovstat(nstat,nstat), hmstat(nstat,nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'ovstat, hmstat', 'read_data')

  allocate (stvec(nstat,nstat), stval(nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'stvec, stval', 'read_data')

  allocate (mdvec(ndspin*ndsptl*ndcplx*mxdtst,nstat), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'mdvec', 'read_data')

  ! fill in arrays

  read (idat) ((ovdet(j,k), j = 1, ndtt1), k = 1, ndtt1)
  read (idat) ((hmdet(j,k), j = 1, ndtt1), k = 1, ndtt1)
  read (idat) ((ovstat(j,k), j = 1, nstat1), k = 1, nstat1)
  read (idat) ((hmstat(j,k), j = 1, nstat1), k = 1, nstat1)
  read (idat) ((stvec(j,k), j = 1, nstat1), k = 1, nstat1)
  read (idat) (stval(k), k = 1, nstat1)

  ! pad appropriately if ndet > ndet1

  if ( ndet1 < ndet ) then
    ovdet(ndtt1+1:ndtt,1:ndtt) = z0
    ovdet(1:ndtt,ndtt1+1:ndtt) = z0

    hmdet(ndtt1+1:ndtt,1:ndtt) = z0
    hmdet(1:ndtt,ndtt1+1:ndtt) = z0
  end if

  if ( nstat1 < nstat ) then
    ovstat(nstat1+1:nstat,1:nstat) = z0
    ovstat(1:nstat,nstat1+1:nstat) = z0

    hmstat(nstat1+1:nstat,1:nstat) = z0
    hmstat(1:nstat,nstat1+1:nstat) = z0

    stvec(nstat1+1:nstat,1:nstat) = z0
    stvec(1:nstat,nstat1+1:nstat) = z0

    stval(nstat1+1:nstat) = d0
  end if

  ! CI expansion of states among its dets

  read (idat) ((mdvec(j,k), j = 1, ndspin*ndsptl*ndcplx*mxdtst1), &
             & k = 1, nstat1)

  if ( mxdtst1 < mxdtst ) then
    mdvec(ndspin*ndsptl*ndcplx*mxdtst1+1:&
         &ndspin*ndsptl*ndcplx*mxdtst,1:nstat) = z0
  end if

  if ( nstat1 < nstat ) then
    mdvec(1:ndspin*ndsptl*ndcplx*mxdtst,nstat1+1:nstat) = z0
  end if

  ! read array of determinants

  allocate (darr(xdim,ndet), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'darr', 'read_data')

  do k = 1, ndet1
    read (idat) (darr(j,k), j = 1, xdim)
  end do

  if ( ndet1 < ndet ) then
    darr(1:xdim,ndet1+1:ndet) = z0
  end if

  allocate (enarr(edim,ndet), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'enarr', 'read_data')

  do k = 1, ndet1
    read (idat) (enarr(j,k), j = 1, edim)
  end do

  if ( ndet1 < ndet ) then
    enarr(1:edim,ndet1+1:ndet) = d0
  end if


  ! read nflag and decide whether to proceed

  read (idat) nflg1

  if ( ( nstat > nstat1 .and. nflg1 /= 0 ) .or. &
     & ( ndet > ndet1 .and. nflg1 /= 0 .and. iopt == 1 ) ) then
    write (6, *) 'error: Aborting calculation.'
    write (6, *) 'It appears like the previous calculation did not converge.'
    stop
  end if


  ! close data file

  close (unit = idat)


  return
end subroutine read_data



subroutine check_xmat (zmat, ltst)

! +----------------------------------------------------------------+
! |                                                                |
! | check_xmat  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Small subroutine to test whether the input matrix zmat can be  |
! | used as a transformation matrix [ = S^(-1/2) ].                |
! |                                                                |
! | In particular, this routine tests whether                      |
! |                                                                |
! |   Z! . S . Z  -  I  =  0,                                      |
! |                                                                |
! | where Z (zmat) is the input matrix (dim nbas x norb), and S    |
! | (smat) is the overlap matrix.                                  |
! |                                                                |
! | The output variable ltst determines whether the test was       |
! | successful (.true.) or not (.false.).                          |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   zmat - matrix to test as Z!.S.Z
  !   ltst - whether Z!.S.Z - I = 0

  logical, intent(out) :: ltst

  complex(kind=dp), dimension(nbas,norb), intent(in) :: zmat


  ! other variables

  integer :: j, k

  complex(kind=dp), dimension(:,:), allocatable :: scr1, scr2


  ! some parameters

  !   thresh - threshold to use in testing

  real(kind=dp) :: thresh
  parameter ( thresh = 1.0e-8_dp )


  allocate (scr1(nbas,norb), scr2(norb,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scratch', 'check_xmat')


  ! compute scr2 = Z! . S .Z

  call zgemm ('n', 'n', nbas, norb, nbas, z1, smat, nbas, &
       & zmat, nbas, z0, scr1, nbas)
  call zgemm ('c', 'n', norb, norb, nbas, z1, zmat, nbas, &
       & scr1, nbas, z0, scr2, norb)

  ! test whether scr2 - I = 0

  do k = 1, norb
    scr2(k,k) = scr2(k,k) - z1
  end do

  ltst = .true.

  do j = 1, norb
    do k = 1, norb
      if ( abs(scr2(k,j)) > thresh ) then
        ltst = .false.
        exit
      end if
    end do
  end do


  deallocate (scr1, scr2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch', 'check_xmat')


  return
end subroutine check_xmat



subroutine save_data

! +----------------------------------------------------------------+
! |                                                                |
! | save_data  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Saves all important variables from current run into a data     |
! | file so that the current calculation can be restarted or a     |
! | different calculation can be launched with additional dets.    |
! |                                                                |
! +----------------------------------------------------------------+

  integer :: j, k


  ! open data file

  inquire (file = datfile, exist = ftest)

  if ( ftest ) then
    open (unit = idat, file = datfile, status = 'old', form = 'unformatted')
    close (unit = idat, status = 'delete')
  end if

  open (unit = idat, file = datfile, status = 'new', form = 'unformatted')

  inquire (file = "detfile", exist = ftest)

  if ( ftest ) then
    open (unit = idat2, file = "detfile", status = 'old', form = 'formatted')
    close (unit = idat2, status = 'delete')
  end if

  open (unit = idat2, file = "detfile", status = 'new', form = 'formatted')

  ! save several integer scalars

  write (idat) nbas, nbct, norb, nup, ndn, iwfnty, imethd
  write (idat2,*) nbas, nbct, norb, nup, ndn, iwfnty, imethd

  ! save spin projection parameters

  if ( mod(imethd-1,2) == 1 ) then
    write (idat) imult
    write (idat) ngrda, ngrdb, ngrdg
    write (idat) igrda, igrdb, igrdg
  end if

  ! save spatial projection parameters

  if ( mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3 ) then
    write (idat) sptgr
    write (idat) sptirr
  end if

  ! save xmat

  write (idat) ((xmat(j,k), j = 1, nbas), k = 1, norb)
  write (idat2,*) ((xmat(j,k), j = 1, nbas), k = 1, norb)


  ! save determinant information

  write (idat) ndet, nstat
  write (idat) (idetst(k), k = 1, nstat)
  write (idat2,*) ndet, nstat
  write (idat2,*) (idetst(k), k = 1, nstat)

  ! save overlap, Hamiltonian matrices
  ! save eigenvectors and eigenvalues between states

  write (idat) ((ovdet(j,k), j = 1, ndtt), k = 1, ndtt)
  write (idat) ((hmdet(j,k), j = 1, ndtt), k = 1, ndtt)
  write (idat) ((ovstat(j,k), j = 1, nstat), k = 1, nstat)
  write (idat) ((hmstat(j,k), j = 1, nstat), k = 1, nstat)
  write (idat) ((stvec(j,k), j = 1, nstat), k = 1, nstat)
  write (idat) (stval(k), k = 1, nstat)

  write (idat2,*) ((ovdet(j,k), j = 1, ndtt), k = 1, ndtt)
  write (idat2,*) ((hmdet(j,k), j = 1, ndtt), k = 1, ndtt)
  write (idat2,*) ((ovstat(j,k), j = 1, nstat), k = 1, nstat)
  write (idat2,*) ((hmstat(j,k), j = 1, nstat), k = 1, nstat)
  write (idat2,*) ((stvec(j,k), j = 1, nstat), k = 1, nstat)
  write (idat2,*) (stval(k), k = 1, nstat)


  ! save CI expansion of states among its dets

  write (idat) ((mdvec(j,k), j = 1, ndspin*ndsptl*ndcplx*mxdtst), &
              & k = 1, nstat)
  write (idat2,*) ((mdvec(j,k), j = 1, ndspin*ndsptl*ndcplx*mxdtst), &
              & k = 1, nstat)
  
!gomez59
  inquire (file = 'dets.info', exist = ftest)
  if ( ftest ) then
    open (unit = 7, file = 'dets.info', status = 'old')
    close (unit = 7, status = 'delete')
  end if
  open (unit = 7, file = 'dets.info', status='new')
  write (7,*) "CI COEFFICIENTS"
  write (7,*) mdvec


  ! save array of determinants

  do k = 1, ndet
    write (idat) (darr(j,k), j = 1, xdim)
    write (idat2,*) (darr(j,k), j = 1, xdim)
  end do

  do k = 1, ndet
    write (idat2,*) (enarr(j,k), j = 1, edim)
  end do

  ! save nflag

  write (idat) nflag
  write (idat2,*) nflag


  ! close data file
  close(idat)
  close(idat2)

!gomez59


  return
end subroutine save_data



subroutine save_det

! +----------------------------------------------------------------+
! |                                                                |
! | save_det  --  CAJH, 01.2013                                    |
! |                                                                |
! |                                                                |
! | Save the list determinants in darr into file 'phfrun.det'.     |
! |                                                                |
! | Note that the determinants are saved in the std AO basis for   |
! | proper interaction with guess_read or other programs. This is  |
! | accomplished by letting                                        |
! |                                                                |
! |    D  <-  X . D,                                               |
! |                                                                |
! | where X is the transformation matrix [ = S^(-1/2) ].           |
! |                                                                |
! +----------------------------------------------------------------+

  integer :: j, k, xdimx, nosq
  integer :: istatus

  complex(kind=dp), dimension(:), allocatable :: dscr

  complex(kind=dp), dimension(:,:), allocatable :: scr_s1, scr_s2
  complex(kind=dp), dimension(:,:), allocatable :: dmat_s1, dmat_s2

  character(len=*) :: detfile
  parameter ( detfile = 'phfrun.det' )


  nosq = norb*norb

  ! open detfile file

  inquire (file = detfile, exist = ftest)

  if ( ftest ) then
    open (unit = idat, file = detfile, status = 'old', form = 'formatted')
    open (unit = 7, file = "dets.info", status = 'old', form = 'formatted')
    close (unit = idat, status = 'delete')
    close (unit = 7, status = 'delete')
  end if

  open (unit = idat, file = detfile, status = 'new', form = 'formatted')
  open (unit = 7, file = "dets.info", status = 'new', form = 'formatted')

  ! save header variables
  write (idat,*) iwfnty, nbas, norb, ndet

!gomez59
  write (7,*) "INFO"
  write (7,*) iwfnty, nbas, norb, ndet 

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
  if ( istatus /= 0 ) call error_alloc (1, 'dscr', 'save_det')

  ! scratch arrays for GHF wfns

  if ( iwfnty == 3 ) then
    allocate (dmat_s1(2*nbas,2*norb), dmat_s2(2*norb,2*norb), &
            & scr_s1(nbas,norb), scr_s2(norb,norb), stat=istatus)
  else
    allocate (dmat_s1(1,1), dmat_s2(1,1), &
            & scr_s1(1,1), scr_s2(1,1), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scratch', 'save_det')


  ! loop over determinants in darr

  do j = 1, ndet

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

    write (idat,*) (dscr(k), k = 1, xdimx)
    write (7,*) 'Determinant: ', j
    write (7,*) (dscr(k), k = 1, xdimx)
  end do

  ! save determinant energies

  do k = 1, ndet
    write (idat,*) (enarr(j,k), j = 1, edim)
  end do


  ! deallocate space for D in std AO basis

  deallocate (dscr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dscr', 'save_det')

  ! deallocate scratch arrays

  deallocate (dmat_s1, dmat_s2, scr_s1, scr_s2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch', 'save_det')


  ! close det file

  close (unit = idat)
  close(unit=7)


  return
end subroutine save_det



subroutine build_new_dmat

! +----------------------------------------------------------------+
! |                                                                |
! | build_new_dmat  --  CAJH, 01.2013                              |
! |                                                                |
! |                                                                |
! | Construct an initial guess for the orbital coefficients (dmat) |
! | of an additional determinant.                                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: istatus, i, ind


  ! allocate memory for matrix dmat

  allocate (dmat(xdim), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'dmat', 'build_new_dmat')

  ! generate a guess for dmat

  if ( igstyp == 0 ) then
    call guess_read (fdet, norb, nbas, nup, ndn, xdim, dmat, iwfnty, &
         & igsnbf, xmat, smat)

  else if ( igstyp == 1 ) then
    call guess_hmat (norb, nbas, nup, ndn, xdim, dmat, iwfnty, xmat, hmat)
!  end if
!gomez59
!initalize density matrix as identity
   else if ( igstyp == 2 ) then
     if ( iwfnty /= 1 ) then
       write (6,*) 'error: not RHF wavefunction'
       stop
     end if
     dmat(:) = z0
     do i = 1, norb
       ind = (i-1)*norb + i
       dmat(ind) = z1
     end do
   end if

  ! mix orbitals if requested

  if ( igsmix > 0 ) then
    call guess_mix (norb, nup, ndn, xdim, dmat, iwfnty, igsmix, igsnmx)
  end if


  return
end subroutine build_new_dmat



subroutine phfmol_send_par

  integer :: ierr, ipos
  integer :: istatus

  integer, dimension(100) :: ibuf

  integer, dimension(5) :: lens, locs, typs
  integer :: mpi_i2s, baseadd

  integer :: mint, mdp, mdc, mchr, mlog, mwrld

  parameter ( mint = mpi_integer )
  parameter ( mdp = mpi_double_precision )
  parameter ( mdc = mpi_double_complex )
  parameter ( mchr = mpi_character )
  parameter ( mlog = mpi_logical )
  parameter ( mwrld = mpi_comm_world )


  ! broadcast fundamental variables

  if ( procid == root ) then
    ipos = 0

    call mpi_pack (nthrd,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (iprint, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (nup,    1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ndn,    1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (imethd, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (iwfnty, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (nbas,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (nbct,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (norb,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ni2s,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (lsf,    1, mlog, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (enr,    1, mdp,  ibuf, 100, ipos, mwrld, ierr)
  end if

  call mpi_bcast (ibuf, 100, mpi_packed, root, mwrld, ierr)

  if ( procid /= root ) then
    ipos = 0

    call mpi_unpack (ibuf, 100, ipos, nthrd,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, iprint, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, nup,    1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ndn,    1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, imethd, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, iwfnty, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, nbas,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, nbct,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, norb,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ni2s,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, lsf,    1, mlog, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, enr,    1, mdp,  mwrld, ierr)

    call omp_set_num_threads (nthrd)
  end if

  ! broadcast one-electron integrals
  ! broadcast xmat

  if ( procid /= root ) then
    allocate (smat(nbas,nbas), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'smat', 'phfmol_send_par')

    allocate (xmat(nbas,norb), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'xmat', 'phfmol_send_par')

    allocate (hmat(nbas,nbas), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'hmat', 'phfmol_send_par')
  end if

  call mpi_bcast (smat, nbas*nbas, mdc, root, mwrld, ierr)
  call mpi_bcast (xmat, nbas*norb, mdc, root, mwrld, ierr)
  call mpi_bcast (hmat, nbas*nbas, mdc, root, mwrld, ierr)


  ! broadcast two-electron integrals

  if ( procid /= root ) then
    allocate (i2sv(ni2s), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'i2sv', 'phfmol_send_par')
  end if

    ! prepare structure

    call mpi_address (i2sv(1), baseadd, ierr)

    call mpi_address(i2sv(1)%ii, locs(1), ierr)

    lens(1) = 1
    locs(1) = locs(1) - baseadd
    typs(1) = mpi_integer2

    call mpi_address(i2sv(1)%ij, locs(2), ierr)

    lens(2) = 1
    locs(2) = locs(2) - baseadd
    typs(2) = mpi_integer2

    call mpi_address(i2sv(1)%ik, locs(3), ierr)

    lens(3) = 1
    locs(3) = locs(3) - baseadd
    typs(3) = mpi_integer2

    call mpi_address(i2sv(1)%il, locs(4), ierr)

    lens(4) = 1
    locs(4) = locs(4) - baseadd
    typs(4) = mpi_integer2

    call mpi_address(i2sv(1)%int, locs(5), ierr)

    lens(5) = 1
    locs(5) = locs(5) - baseadd
    typs(5) = mpi_double_precision

    call mpi_type_struct (5, lens, locs, typs, mpi_i2s, ierr)
    call mpi_type_commit (mpi_i2s, ierr)

  call mpi_bcast (i2sv, ni2s, mpi_i2s, root, mwrld, ierr)


  ! broadcast projection variables

  if ( procid == root ) then
    ipos = 0

    call mpi_pack (imult,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ngrda,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ngrdb,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ngrdg,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (igrda,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (igrdb,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (igrdg,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (sptirr, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (sptgr,  4, mchr, ibuf, 100, ipos, mwrld, ierr)
  end if

  call mpi_bcast (ibuf, 100, mpi_packed, root, mwrld, ierr)

  if ( procid /= root ) then
    ipos = 0

    call mpi_unpack (ibuf, 100, ipos, imult,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ngrda,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ngrdb,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ngrdg,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, igrda,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, igrdb,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, igrdg,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, sptirr, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, sptgr,  4, mchr, mwrld, ierr)
  end if


  ! broadcast variables related to number of determinants

  if ( procid == root ) then
    ipos = 0

    call mpi_pack (ndspin, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ndsptl, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ndcplx, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (xdim,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (edim,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (ndet,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (nstat,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (iopt,   1, mint, ibuf, 100, ipos, mwrld, ierr)
  end if

  call mpi_bcast (ibuf, 100, mpi_packed, root, mwrld, ierr)

  if ( procid /= root ) then
    ipos = 0

    call mpi_unpack (ibuf, 100, ipos, ndspin, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ndsptl, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ndcplx, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, xdim,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, edim,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, ndet,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, nstat,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, iopt,   1, mint, mwrld, ierr)
  end if

  ! broadcast idetst array

  if ( procid /= root ) then
    allocate (idetst(nstat), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'idetst', 'phfmol_send_par')
  end if

  call mpi_bcast (idetst(1), nstat, mint, root, mwrld, ierr)

  ! allocate other arrays

  if ( procid /= root ) then
    ndtt   = ndet*ndspin*ndsptl*ndcplx
    mxdtst = maxval (idetst(1:nstat))

    if ( iopt == 1 ) then
      ndtopt = 1
    else
      ndtopt = idetst(nstat)
    end if

    allocate (ovdet(ndtt,ndtt), hmdet(ndtt,ndtt), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'ovdet, hmdet', 'phfmol_send_par')

    allocate (ovstat(nstat,nstat), hmstat(nstat,nstat), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'ovstat, hmstat', 'phfmol_send_par')

    allocate (stvec(nstat,nstat), stval(nstat), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'stvec, stval', 'phfmol_send_par')

    allocate (mdvec(ndspin*ndsptl*ndcplx*mxdtst,nstat), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'mdvec', 'phfmol_send_par')
  end if

  ! .. no need to broadcast overlap and Hamiltonian matrices ..


  ! broadcast array of D matrices

  if ( procid /= root ) then
    allocate (darr(xdim,ndet), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'darr', 'phfmol_send_par')

    allocate (enarr(edim,ndet), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'enarr', 'phfmol_send_par')
  end if

  call mpi_bcast (darr(1,1),  xdim*ndet, mdc, root, mwrld, ierr)
  call mpi_bcast (enarr(1,1), edim*ndet, mdp, root, mwrld, ierr)


  ! broadcast optimization variables

  if ( procid == root ) then
    ipos = 0

    call mpi_pack (maxit,  1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (mxfunc, 1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (nbuf,   1, mint, ibuf, 100, ipos, mwrld, ierr)
    call mpi_pack (gtol,   1, mdp,  ibuf, 100, ipos, mwrld, ierr)
  end if

  call mpi_bcast (ibuf, 100, mpi_packed, root, mwrld, ierr)

  if ( procid /= root ) then
    ipos = 0

    call mpi_unpack (ibuf, 100, ipos, maxit,  1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, mxfunc, 1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, nbuf,   1, mint, mwrld, ierr)
    call mpi_unpack (ibuf, 100, ipos, gtol,   1, mdp,  mwrld, ierr)
  end if


  return
end subroutine phfmol_send_par



subroutine shutdown_phfmol

  integer :: istatus

  deallocate (idetst, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'idetst', 'shutdown_phfmol')

  deallocate (ovdet, hmdet, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ovdet, hmdet', 'shutdown_phfmol')

  deallocate (ovstat, hmstat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'ovstat, hmstat', 'shutdown_phfmol')

  deallocate (stvec, stval, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'stvec, stval', 'shutdown_phfmol')

  deallocate (mdvec, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'mdvec', 'shutdown_phfmol') 
  deallocate (darr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'darr', 'shutdown_phfmol')

  if ( allocated (dmat) ) then
    deallocate (dmat, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'dmat', 'shutdown_phfmol')
  end if

  deallocate (smat, xmat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'smat, xmat', 'shutdown_phfmol')

  deallocate (hmat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'hmat', 'shutdown_phfmol')

  deallocate (i2sv, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'i2sv', 'shutdown_phfmol')


  return
end subroutine shutdown_phfmol



subroutine load_purcart

! +----------------------------------------------------------------+
! |                                                                |
! | load_purcart  --  CAJH, 02.2013                                |
! |                                                                |
! |                                                                |
! | Read, from file 'purcart.dat', the transformation matrices     |
! | required for pure <-> Cartesian transformations.               |
! |                                                                |
! | The matrices are loaded into the module purcart to be          |
! | accessible throughout the rest of the program.                 |
! |                                                                |
! | The file 'purcart.dat' should have the following structure:    |
! |                                                                |
! |   =>  nbas, nbct                                               |
! |   =>  x2pure (real, nbct x nbas matrix)                        |
! |   =>  x2cart (real, nbas x nbct matrix)                        |
! |                                                                |
! +----------------------------------------------------------------+

  integer :: j, k, ierr, istatus
  integer :: nbas_old, nbct_old

  real(kind=dp), dimension(:,:), allocatable :: x2pr, x2cr
  complex(kind=dp), dimension(:,:), allocatable :: x2pure, x2cart


  ! allocate memory for pure <-> Cartesian transformation matrices

  if ( procid == root ) then
    allocate (x2pr(nbct,nbas), x2cr(nbas,nbct), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'x2r', 'load_purcart')
  end if

  allocate (x2pure(nbct,nbas), x2cart(nbas,nbct), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'x2c', 'load_purcart')


  ! read transformation matrices

  if ( procid == root ) then
    inquire (file = 'purcart.dat', exist = ftest)
    
    if ( .not. ftest ) then
      write (6, *) 'error: File purcart.dat cannot be found.'
      stop
    end if

    open (unit = idat, file = 'purcart.dat', status = 'old', &
        & form = 'unformatted')

    read (idat) nbas_old, nbct_old

    if ( nbas_old /= nbas ) then
      write (6, *) 'error: Inconsistent nbas and nbas_old in load_purcart.'
      stop
    end if

    if ( nbct_old /= nbct ) then
      write (6, *) 'error: Inconsistent nbct and nbct_old in load_purcart.'
      stop
    end if

    read (idat) ((x2pr(j,k), j = 1, nbct), k = 1, nbas)
    read (idat) ((x2cr(j,k), j = 1, nbas), k = 1, nbct)

    close (unit = idat)
  end if


  ! transform matrices to complex

  if ( procid == root ) then
    x2pure(1:nbct,1:nbas) = x2pr(1:nbct,1:nbas)
    x2cart(1:nbas,1:nbct) = x2cr(1:nbas,1:nbct)
  end if


  ! broadcast transformation matrices

  call mpi_bcast (x2pure, nbct*nbas, mpi_double_complex, root, &
       & mpi_comm_world, ierr)
  call mpi_bcast (x2cart, nbct*nbas, mpi_double_complex, root, &
       & mpi_comm_world, ierr)


  ! load module purcart to store transformation matrices

  !$omp  parallel default(shared)

  call setup_purcart (nbas, nbct, x2pure, x2cart)

  !$omp  end parallel


  ! deallocate memory

  deallocate (x2pure, x2cart, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'x2c', 'load_purcart')

  if ( procid == root ) then
    deallocate (x2pr, x2cr, stat=istatus)
    if ( istatus /= 0 ) call error_alloc (2, 'x2r', 'load_purcart')
  end if


  return
end subroutine load_purcart


end program phfmol


