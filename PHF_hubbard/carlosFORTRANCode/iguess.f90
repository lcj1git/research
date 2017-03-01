

module iguess

  use constants
  use util

  implicit none
  private

  public :: guess_read, guess_hmat, guess_mix

! +----------------------------------------------------------------+
! |                                                                |
! | iguess                                                         |
! |                                                                |
! |                                                                |
! | A collection of subroutines related to the preparation of an   |
! | initial guess for an additional determinant in phfmol.         |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The matrix of orbital coefficients, D, characterizes the       |
! | linear canonical transformation from the basis states to the   |
! | HF states. In this code, we take our basis states to be the    |
! | ort AO basis ones.                                             |
! |                                                                |
! | Explicitly, we use the following convention                    |
! |                                                                |
! |   b!_a  =  sum_{j}  D*_{ja}  c!_j,                             |
! |                                                                |
! | where b are HF operators and c are basis states operators.     |
! | The matrix D is unitary ( D . D! = I ). As it is commonly      |
! | done, the first N columns in D correspond to the occupied      |
! | (hole) states, while the last 2*norbs-N columns correspond to  |
! | the unoccupied (particle) states.                              |
! |                                                                |
! | We can understand D as the transformation matrix               |
! |                                                                |
! |   D*_{ja}  =  <a|j>   ==>   D_{ja}  =  <j|a>                   |
! |                                                                |
! +----------------------------------------------------------------+

contains


subroutine guess_read (fdet, norb, nbas, nup, ndn, xdim, dmat, iwfnty, &
     & igsnbf, xmat, smat)

! +----------------------------------------------------------------+
! |                                                                |
! | guess_read  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Prepare an initial guess of the matrix D of orbital            |
! | coefficients by reading an old one stored in file [ fdet ].    |
! |                                                                |
! | We assume that the wavefunctions stored in [ fdet ] are        |
! | in the std AO basis. Two actions are performed here:           |
! |                                                                |
! |   1. Make sure that, in the std AO basis, we have D!.S.D = I,  |
! |      where S is the overlap matrix. If this is not the case,   |
! |      we project the orbitals in order to restore               |
! |      orthonormality wrt S:                                     |
! |                                                                |
! |        =>   D!. S . D  =  A                                    |
! |        =>   D'! . S . D'  =  I,                                |
! |                                                                |
! |      where D' = D . A^(-1/2).                                  |
! |                                                                |
! |   2. Transform the stored wavefunction to the ort AO basis     |
! |      using the transformation matrix xmat.                     |
! |                                                                |
! |        =>   D! . S . D  =  I        [ std AO ]                 |
! |        =>   D'! . D'  =  I,         [ ort AO ]                 |
! |                                                                |
! |      where D' = X! . S . D.                                    |
! |                                                                |
! | Switching of wavefunction types (iwfnty) is allowed as long as |
! | no information is lost. We do require, however, that the       |
! | dimension of the wavefunctions (nbas, norbs) is consistent.    |
! |                                                                |
! | In general, there are ndet determinants stored in [ fdet ].    |
! | The input variable igsnbf can be used to determine which one   |
! | to load. If set to 0, the last one will be loaded.             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   fdet   - name of file to load
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   xdim   - leading dimension of dmat
  !   dmat   - matrix of orbital coefficients [ updated here ]
  !   iwfnty - type of dmat wavefunction
  !   igsnbf - number of determinant to read
  !              ( if 0, the last one will be loaded )
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )

  character(len=*), intent(in) :: fdet

  integer, intent(in) :: norb, nbas, nup, ndn
  integer, intent(in) :: xdim, iwfnty, igsnbf

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat
  complex(kind=dp), dimension(xdim), intent(inout) :: dmat


  ! other variables

  integer :: j, k, xdimx, xdimt, nosq
  integer :: iwfn_old, nbas_old, norb_old
  integer :: iup, idn, ik_up, ik_dn, ik1_up, ik1_dn
  integer :: inum, ndet, istatus

  integer :: iflg, iflg_up, iflg_dn

  complex(kind=dp), dimension(:), allocatable :: dmat_old, dmatt
  complex(kind=dp), dimension(:,:), allocatable :: s12

  complex(kind=dp), dimension(:,:), allocatable :: scr_s1, scr_s2
  complex(kind=dp), dimension(:,:), allocatable :: dmat_s1, dmat_s2


  ! data file for old wavefunction

  integer :: idat
  parameter ( idat = 13 )

  logical :: ftest


  ! error checking

  if ( iwfnty < 1 .or. iwfnty > 3 ) then
    write (6, *) 'error: Wfn type not supported in guess_read.'
    stop
  end if


  ! open guess file
    
  inquire (file = fdet, exist = ftest)
    
  if ( .not. ftest ) then
    write (6, *) 'error: File ', fdet, ' cannot be found.'
    stop
  end if

  open (unit = idat, file = fdet, status = 'old', form = 'formatted')


  ! read old wfn parameters

  read (idat,*) iwfn_old, nbas_old, norb_old, ndet


  ! error checking:
  !   - nbas must be equal to nbas_old
  !   - norb must be equal to norb_old
  !   - iwfnty must be consistent with iwfn_old
  !   - igsnbf must be <= ndet

  if ( nbas_old /= nbas ) then
    write (6, *) 'error: Inconsistent nbas and nbas_old in guess_read.'
    stop
  end if

  if ( norb_old /= norb ) then
    write (6, *) 'error: Inconsistent norb and norb_old in guess_read.'
    stop
  end if

  if ( iwfnty >= 1 .and. iwfnty <= 3 ) then
    if ( iwfn_old > iwfnty ) then
      write (6, *) 'error: Inconsistent wfn and wfn_old in guess_read.'
      stop
    end if
  end if

  if ( igsnbf > ndet ) then
    write (6, *) 'error: No. of dets in [ fdet ] is smaller than igsnbf.'
    stop
  end if

  inum = igsnbf
  if ( inum == 0 ) inum = ndet


  ! allocate memory for old determinant

  xdimx = 0

  if ( iwfn_old == 1 ) then
    xdimx = nbas*norb
  else if ( iwfn_old == 2 ) then
    xdimx = 2*nbas*norb
  else if ( iwfn_old == 3 ) then
    xdimx = 4*nbas*norb
  end if

  allocate (dmat_old(xdimx), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'dmat_old', 'guess_read')


  ! read old wfn

  do j = 1, inum-1
    read (idat,*)
  end do

  read (idat,*) (dmat_old(k), k = 1, xdimx)

  close (unit = idat)  ! close guess file


  ! project orbitals into current std AO basis

    iflg = 0
    iflg_up = 0;  iflg_dn = 0

  if ( iwfn_old == 1 ) then
    call gsproj1 (norb, nbas, xdimx, dmat_old, smat, iflg)

  else if ( iwfn_old == 2 ) then
    call gsproj1 (norb, nbas, xdimx/2, &
                      & dmat_old(1:norb*nbas), smat, iflg_up)
    call gsproj1 (norb, nbas, xdimx/2, &
                      & dmat_old(norb*nbas+1:2*norb*nbas), smat, iflg_dn)

  else if ( iwfn_old == 3 ) then
    call gsproj2 (norb, nbas, xdimx, dmat_old, smat, iflg)
  end if

  if ( iwfn_old == 1 .or. iwfn_old == 3 ) then
    if ( iflg /= 0 ) then
      write (6, *) 'error: Projection of guess orbitals failed.'
      stop
    end if

  else if ( iwfn_old == 2 ) then
    if ( iflg_up /= 0 .or. iflg_dn /= 0 ) then
      write (6, *) 'error: Projection of guess orbitals failed.'
      stop
    end if
  end if


  ! +-----------------------------+
  ! |  transform to ort AO basis  |
  ! +-----------------------------+

  nosq  = norb*norb
  xdimt = 0

  if ( iwfn_old == 1 ) then
    xdimt = nosq
  else if ( iwfn_old == 2 ) then
    xdimt = 2*nosq
  else if ( iwfn_old == 3 ) then
    xdimt = 4*nosq
  end if

  allocate (dmatt(xdimt), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'dmatt', 'guess_read')


  ! scratch arrays for GHF wfns

  if ( iwfn_old == 3 ) then
    allocate (dmat_s1(2*nbas,2*norb), dmat_s2(2*norb,2*norb), &
            & scr_s1(nbas,norb), scr_s2(norb,norb), stat=istatus)
  else
    allocate (dmat_s1(1,1), dmat_s2(1,1), &
            & scr_s1(1,1), scr_s2(1,1), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scratch', 'guess_read')


  ! build X! . S;  store in s12

  allocate (s12(norb,nbas), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 's12', 'guess_read')

  call zgemm ('c', 'n', norb, nbas, nbas, z1, xmat, nbas, smat, &
       & nbas, z0, s12, norb)


  ! form the new orbitals according to
  !   D'  =  X! . S . D  =  s12 . D

  ! store new orbitals in dmatt

  if ( iwfn_old == 1 ) then
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & dmat_old, nbas, z0, dmatt, norb)

  else if ( iwfn_old == 2 ) then
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & dmat_old(1), nbas, z0, dmatt(1), norb)
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & dmat_old(norb*nbas+1), nbas, z0, dmatt(nosq+1), norb)

  else if ( iwfn_old == 3 ) then
    dmat_s1(1:2*nbas,1:2*norb) = &
         & reshape (dmat_old(1:4*nbas*norb), (/ 2*nbas, 2*norb /))

      scr_s1(1:nbas,1:norb) = dmat_s1(1:nbas,1:norb)
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & scr_s1, nbas, z0, scr_s2, norb)
      dmat_s2(1:norb,1:norb) = scr_s2(1:norb,1:norb)

      scr_s1(1:nbas,1:norb) = dmat_s1(1:nbas,norb+1:2*norb)
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & scr_s1, nbas, z0, scr_s2, norb)
      dmat_s2(1:norb,norb+1:2*norb) = scr_s2(1:norb,1:norb)

      scr_s1(1:nbas,1:norb) = dmat_s1(nbas+1:2*nbas,1:norb)
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & scr_s1, nbas, z0, scr_s2, norb)
      dmat_s2(norb+1:2*norb,1:norb) = scr_s2(1:norb,1:norb)

      scr_s1(1:nbas,1:norb) = dmat_s1(nbas+1:2*nbas,norb+1:2*norb)
    call zgemm ('n', 'n', norb, norb, nbas, z1, s12, norb, &
         & scr_s1, nbas, z0, scr_s2, norb)
      dmat_s2(norb+1:2*norb,norb+1:2*norb) = scr_s2(1:norb,1:norb)

    dmatt(1:4*nosq) = &
         & reshape (dmat_s2(1:2*norb,1:2*norb), (/ 4*nosq /))
  end if


  ! +-------------------------------+
  ! |  fill new dmat using old one  |
  ! +-------------------------------+

  if ( iwfn_old == 1 ) then

    if ( iwfnty == 1 ) then
      dmat(1:nosq) = dmatt(1:nosq)

    else if ( iwfnty == 2 ) then
      dmat(1:nosq)        = dmatt(1:nosq)
      dmat(nosq+1:2*nosq) = dmatt(1:nosq)

    else if ( iwfnty == 3 ) then
      iup = 1
      idn = 1

      do k = 1, 2*norb
        ik_up = (k-1)*2*norb
        ik_dn = (k-1)*2*norb + norb

        ik1_up = (iup-1)*norb
        ik1_dn = (idn-1)*norb

        if ( k <= nup ) then
          dmat(ik_up+1:ik_up+norb) = dmatt(ik1_up+1:ik1_up+norb)
          dmat(ik_dn+1:ik_dn+norb) = z0
          iup = iup + 1

        else if ( k > nup .and. k <= nup+ndn ) then
          dmat(ik_up+1:ik_up+norb) = z0
          dmat(ik_dn+1:ik_dn+norb) = dmatt(ik1_dn+1:ik1_dn+norb)
          idn = idn + 1

        else if ( k > nup+ndn .and. k <= norb+ndn ) then
          dmat(ik_up+1:ik_up+norb) = dmatt(ik1_up+1:ik1_up+norb)
          dmat(ik_dn+1:ik_dn+norb) = z0
          iup = iup + 1

        else if ( k > norb+ndn .and. k <= 2*norb ) then
          dmat(ik_up+1:ik_up+norb) = z0
          dmat(ik_dn+1:ik_dn+norb) = dmatt(ik1_dn+1:ik1_dn+norb)
          idn = idn + 1
        end if
      end do

    end if

  else if ( iwfn_old == 2 ) then

    if ( iwfnty == 2 ) then
      dmat(1:2*nosq) = dmatt(1:2*nosq)

    else if ( iwfnty == 3 ) then
      iup = 1
      idn = 1

      do k = 1, 2*norb
        ik_up = (k-1)*2*norb
        ik_dn = (k-1)*2*norb + norb

        ik1_up = (iup-1)*norb
        ik1_dn = (idn-1)*norb + nosq

        if ( k <= nup ) then
          dmat(ik_up+1:ik_up+norb) = dmatt(ik1_up+1:ik1_up+norb)
          dmat(ik_dn+1:ik_dn+norb) = z0
          iup = iup + 1

        else if ( k > nup .and. k <= nup+ndn ) then
          dmat(ik_up+1:ik_up+norb) = z0
          dmat(ik_dn+1:ik_dn+norb) = dmatt(ik1_dn+1:ik1_dn+norb)
          idn = idn + 1

        else if ( k > nup+ndn .and. k <= norb+ndn ) then
          dmat(ik_up+1:ik_up+norb) = dmatt(ik1_up+1:ik1_up+norb)
          dmat(ik_dn+1:ik_dn+norb) = z0
          iup = iup + 1

        else if ( k > norb+ndn .and. k <= 2*norb ) then
          dmat(ik_up+1:ik_up+norb) = z0
          dmat(ik_dn+1:ik_dn+norb) = dmatt(ik1_dn+1:ik1_dn+norb)
          idn = idn + 1
        end if
      end do

    end if

  else if ( iwfn_old == 3 ) then

    if ( iwfnty == 3 ) then
      dmat(1:4*nosq) = dmatt(1:4*nosq)
    end if

  end if


  ! deallocate memory

  deallocate (dmatt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dmatt', 'guess_read')

  deallocate (s12, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 's12', 'guess_read')

  deallocate (dmat_old, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'dmat_old', 'guess_read')

  deallocate (dmat_s1, dmat_s2, scr_s1, scr_s2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch', 'guess_read')



  return
end subroutine guess_read



subroutine guess_hmat (norb, nbas, nup, ndn, xdim, dmat, iwfnty, &
     & xmat, hmat)

! +----------------------------------------------------------------+
! |                                                                |
! | guess_hmat  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Prepare an initial guess of the matrix D of orbital            |
! | coefficients by diagonalizing the matrix hmat (core            |
! | Hamiltonian) and occupying the lowest-energy eigenvectors.     |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   xdim   - leading dimension of dmat
  !   dmat   - matrix of orbital coefficients [ updated here ]
  !   iwfnty - type of dmat wavefunction
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   hmat   - core Hamiltonian matrix ( std AO basis )

  integer, intent(in) :: norb, nbas, nup, ndn
  integer, intent(in) :: xdim, iwfnty

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: hmat
  complex(kind=dp), dimension(xdim), intent(inout) :: dmat


  ! other variables

  integer :: nosq
  integer :: nb, lwork
  integer :: k, iup, idn, ik_up, ik_dn
  integer :: istatus, info

  integer :: idum, inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  integer, dimension(:), allocatable :: iwork, ifail

  real(kind=dp), dimension(:), allocatable :: rwrk, val

  complex(kind=dp), dimension(:), allocatable :: work
  complex(kind=dp), dimension(:,:), allocatable :: hmatt, scrt, rvec


  ! external functions

  integer, external :: ilaenv
  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  ! error checking

  if ( iwfnty < 1 .or. iwfnty > 3 ) then
    write (6, *) 'error: Wfn type not supported in guess_hmat.'
    stop
  end if


  ! transform core Hamiltonian to ort AO basis

  allocate (hmatt(norb,norb), scrt(nbas,norb), rvec(norb,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'hmatt', 'guess_hmat')

  call zgemm ('n', 'n', nbas, norb, nbas, z1, hmat, nbas, &
       & xmat, nbas, z0, scrt, nbas)
  call zgemm ('c', 'n', norb, norb, nbas, z1, xmat, nbas, &
       & scrt, nbas, z0, hmatt, norb)


  ! +---------------------+
  ! |  diagonalize hmatt  |
  ! +---------------------+

  ! allocate memory for eigenvalues of hmatt (val)

  allocate (val(norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'val', 'guess_hmat')

  ! determine optimal lwork for zheev
  !   ( this code was extracted directly from zheev )

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', norb, -1, -1, -1)
  lwork = max (1, (nb+1)*norb)

  ! allocate memory for scratch arrays (zheev)

  allocate (iwork(5*norb), ifail(norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'guess_hmat')

  allocate (work(lwork), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'work', 'guess_hmat')

  ! allocate (rwrk(3*norb-2), stat=istatus)
  allocate (rwrk(7*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rwrk', 'guess_hmat')


  ! diagonalize hmatt:
  !   eigenvalues are stored in val
  !   eigenvectors are stored in hmatt

  call zheevx ('v', 'a', 'u', norb, hmatt, norb, ddum, ddum, &
       & idum, idum, abstol, inum, val, rvec, norb, work, &
       & lwork, rwrk, iwork, ifail, info)

  hmatt(1:norb,1:norb) = rvec(1:norb,1:norb)

  ! call zheev ('v', 'u', norb, hmatt, norb, val, work, lwork, &
  !      & rwrk, info)

  if ( info /= 0 ) then
    write (6, *) 'error: zheev error in guess_hmat.'
    stop
  end if


  ! fill dmat with eigenvectors of hmatt

  nosq = norb*norb

  if ( iwfnty == 1 ) then

    dmat(1:nosq) = &
         & reshape (conjg (hmatt(1:norb,1:norb)), (/ nosq /))

  else if ( iwfnty == 2 ) then

    dmat(1:nosq)        = &
         & reshape (conjg (hmatt(1:norb,1:norb)), (/ nosq /))
    dmat(nosq+1:2*nosq) = &
         & reshape (conjg (hmatt(1:norb,1:norb)), (/ nosq /))

  else if ( iwfnty == 3 ) then

    iup = 0
    idn = 0

    do k = 1, 2*norb
      ik_up = (k-1)*2*norb
      ik_dn = (k-1)*2*norb + norb

      if ( k <= nup ) then
        iup = iup + 1
        dmat(ik_up+1:ik_up+norb) = conjg (hmatt(1:norb,iup))
        dmat(ik_dn+1:ik_dn+norb) = z0

      else if ( k > nup .and. k <= nup+ndn ) then
        idn = idn + 1
        dmat(ik_up+1:ik_up+norb) = z0
        dmat(ik_dn+1:ik_dn+norb) = conjg (hmatt(1:norb,idn))

      else if ( k > nup+ndn .and. k <= norb+ndn ) then
        iup = iup + 1
        dmat(ik_up+1:ik_up+norb) = conjg (hmatt(1:norb,iup))
        dmat(ik_dn+1:ik_dn+norb) = z0

      else if ( k > norb+ndn .and. k <= 2*norb ) then
        idn = idn + 1
        dmat(ik_up+1:ik_up+norb) = z0
        dmat(ik_dn+1:ik_dn+norb) = conjg (hmatt(1:norb,idn))
      end if
    end do

  end if


  ! deallocate memory

  deallocate (hmatt, scrt, rvec, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'hmatt', 'guess_hmat')

  deallocate (val, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'val', 'guess_hmat')

  deallocate (iwork, ifail, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'iwork', 'guess_hmat')

  deallocate (work, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'work', 'guess_hmat')

  deallocate (rwrk, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rwrk', 'guess_hmat')


  return
end subroutine guess_hmat



subroutine guess_mix (norb, nup, ndn, xdim, dmat, iwfnty, itype, igsnmx)

  use linalg, only : zmat_symm
  ! use ifport

! +----------------------------------------------------------------+
! |                                                                |
! | guess_mix  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Mix the orbitals in a previously prepared initial guess of the |
! | matrix D of orbital coefficients. That is, this subroutine     |
! | updates D according to                                         |
! |                                                                |
! |   D  <--  D . X,                                               |
! |                                                                |
! | where X is a unitary matrix built as                           |
! |                                                                |
! |   X  =  exp(i*fc*A).                                           |
! |                                                                |
! | Here, A is a random (complex) Hermitian matrix. fc controls    |
! | the amount of mixing introduced by X; for small fc the mixing  |
! | is small. itype is used to control the factor fc used:         |
! |                                                                |
! |   itype  =  1,  fc = 5.0e-2,  very mild mixing                 |
! |          =  2,     = 1.0e-1,  mild mixing                      |
! |          =  3,     = 2.0e-1,  normal mixing                    |
! |          =  4,     = 4.0e-1,  aggressive mixing                |
! |          =  5,     = 8.0e-1,  very aggressive mixing           |
! |                                                                |
! | The number of occ-virtual orbital pairs to mix (of each spin)  |
! | is controlled by the input variable igsnmx. This must satisfy  |
! |                                                                |
! |   igsnmx <= nup,  igsnmx <= norb-nup                           |
! |   igsnmx <= ndn,  igsnmx <= norb-ndn                           |
! |                                                                |
! | If igsnmx = 0, then all orbitals are mixed, regardless of      |
! | whether they are occupied or virtual.                          |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb   - number of orbitals
  !   nup    - number of spin-up electrons
  !   ndn    - number of spin-dn electrons
  !   xdim   - leading dimension of dmat
  !   dmat   - matrix of orbital coefficients [ updated here ]
  !   iwfnty - type of dmat wavefunction
  !   itype  - type of mixing to perform (see above)
  !   igsnmx - number of orbitals (of each spin) to mix

  integer, intent(in) :: norb, nup, ndn
  integer, intent(in) :: xdim, iwfnty, itype, igsnmx

  complex(kind=dp), dimension(xdim), intent(inout) :: dmat


  ! some parameters

  !   arrfac - stored values for different mixings (see above)

  real(kind=dp), dimension(5) :: arrfac
  parameter ( arrfac = (/ 5.0e-2_dp, 1.0e-1_dp, 2.0e-1_dp, &
                        & 4.0e-1_dp, 8.0e-1_dp /) )


  ! other variables

  integer :: nosq, dima, adim
  integer :: j, k, ind
  integer :: nb, lwork
  integer :: istatus, info

  integer :: szseed, pid, s
  integer, dimension(8) :: dt
  integer, dimension(:), allocatable :: seed

  real(kind=dp) :: fac, fc_r, fc_i

  integer :: idum, inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  integer, dimension(:), allocatable :: iwork, ifail

  real(kind=dp), dimension(:), allocatable :: rwrk, val

  complex(kind=dp), dimension(:), allocatable :: work, scr
  complex(kind=dp), dimension(:), allocatable :: amat, xmat, avec


  ! external functions

  integer, external :: ilaenv
  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  nosq = norb*norb

  ! error checking

  if ( iwfnty < 1 .or. iwfnty > 3 ) then
    write (6, *) 'error: Wfn type not supported in guess_mix.'
    stop
  end if

  if ( itype < 1 .or. itype > 5 ) then
    write (6, *) 'error: itype not recognized in guess_mix.'
    stop
  end if

  fac = arrfac(itype)   ! determine mixing factor


  ! check whether igsnmx is within bounds

  if ( igsnmx < 0 ) then
    write (6, *) 'error: igsnmx < 0 in guess_mix.'
    stop
  end if

  if ( igsnmx > 0 ) then
    if ( iwfnty == 1 ) then
      if ( igsnmx > nup .or. igsnmx > norb-nup ) then
        write (6, *) 'error: igsnmx too large in guess_mix.'
        stop
      end if

    else if ( iwfnty == 2 ) then
      if ( igsnmx > nup .or. igsnmx > norb-nup .or. &
         & igsnmx > ndn .or. igsnmx > norb-ndn ) then
        write (6, *) 'error: igsnmx too large in guess_mix.'
        stop
      end if

    else if ( iwfnty == 3 ) then
      if ( 2*igsnmx > nup+ndn .or. 2*igsnmx > 2*norb-nup-ndn ) then
        write (6, *) 'error: igsnmx too large in guess_mix.'
        stop
      end if
    end if
  end if

  if ( igsnmx == 0 ) then
    dima = norb
  else
    dima = 2*igsnmx
  end if


  ! prepare calls to rand ()

  pid = getpid ()

  call random_seed (size=szseed)
  allocate (seed(szseed))
  seed(:) = 0

  call date_and_time (values=dt)

  seed(1) = dt(8) + 1000*dt(7) + 60*1000*dt(6) + 60*60*1000*dt(5)
  s = ieor(seed(1), pid)

  seed(1) = seed(1) + 32771
  if ( szseed > 1 ) seed(2) = pid
  if ( szseed > 2 ) seed(3:) = s + 521 * (/ (j, j = 1, szseed-2) /)

  call random_seed (put=seed)


  ! +-----------------------------------+
  ! |  build random Hermitian matrix A  |
  ! +-----------------------------------+

  ! determine dimension of matrix A

  adim = 0

  if ( iwfnty == 1 ) then
    adim = dima*dima
  else if ( iwfnty == 2 ) then
    adim = 2*dima*dima
  else if ( iwfnty == 3 ) then
    adim = 4*dima*dima
  end if


  ! allocate memory for amat, xmat
  ! allocate memory for scr (scratch space for matrix operations)

  allocate (amat(adim), avec(adim), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'amat', 'guess_mix')

  allocate (xmat(adim), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'xmat', 'guess_mix')

  if ( iwfnty == 1 .or. iwfnty == 2 ) then
    allocate (scr(norb*dima), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (scr(4*norb*dima), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'scr', 'guess_mix')


  ! fill amat with random numbers

  do k = 1, adim
    call random_number (harvest=fc_r)
    call random_number (harvest=fc_i)

    fc_r = -1.0e0_dp + 2.0e0_dp * fc_r
    fc_i = -1.0e0_dp + 2.0e0_dp * fc_i

    !jag20
    !amat(k) = cmplx (fc_r, fc_i, dp)
     amat(k) = cmplx (d0, fc_i, dp)
  end do

  deallocate (seed)


  ! hermitize amat

  if ( iwfnty == 1 ) then
    call zmat_symm (1, 'u', dima, amat)

  else if ( iwfnty == 2 ) then
    call zmat_symm (1, 'u', dima, amat)
    call zmat_symm (1, 'u', dima, amat(dima*dima+1))

  else if ( iwfnty == 3 ) then
    call zmat_symm (1, 'u', 2*dima, amat)
  end if


  ! +------------------------+
  ! |  diagonalize matrix A  |
  ! +------------------------+

  ! allocate memory for eigenvalues (val) of amat

  if ( iwfnty == 1 ) then
    allocate (val(dima), stat=istatus)
  else if ( iwfnty == 2 .or. iwfnty == 3 ) then
    allocate (val(2*dima), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'val', 'guess_mix')

  ! determine optimal lwork for zheev
  !   ( this code was extracted directly from zheev )

  if ( iwfnty == 1 .or. iwfnty == 2 ) then
    ! nb = 2
    nb = ilaenv (1, 'zhetrd', 'u', dima, -1, -1, -1)
    lwork = max (1, (nb+1)*dima)

  else if ( iwfnty == 3 ) then
    ! nb = 2
    nb = ilaenv (1, 'zhetrd', 'u', 2*dima, -1, -1, -1)
    lwork = max (1, (nb+1)*2*dima)
  end if


  ! allocate memory for scratch arrays (zheev)

  if ( iwfnty == 1 .or. iwfnty == 2 ) then
    allocate (iwork(5*dima), ifail(dima), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'guess_mix')
  else if ( iwfnty == 3 ) then
    allocate (iwork(10*dima), ifail(2*dima), stat=istatus)
    if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'guess_mix')
  end if

  allocate (work(lwork), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'work', 'guess_mix')

  if ( iwfnty == 1 .or. iwfnty == 2 ) then
    ! allocate (rwrk(3*dima-2), stat=istatus)
    allocate (rwrk(7*dima), stat=istatus)
  else if ( iwfnty == 3 ) then
    ! allocate (rwrk(6*dima-2), stat=istatus)
    allocate (rwrk(14*dima), stat=istatus)
  end if

  if ( istatus /= 0 ) call error_alloc (1, 'rwrk', 'guess_mix')


  ! diagonalize amat:
  !   eigenvalues are stored in val
  !   eigenvectors are stored in amat

  if ( iwfnty == 1 ) then
    call zheevx ('v', 'a', 'u', dima, amat, dima, ddum, ddum, &
         & idum, idum, abstol, inum, val, avec, dima, work, &
         & lwork, rwrk, iwork, ifail, info)

    amat(1:adim) = avec(1:adim)

    ! call zheev ('v', 'u', dima, amat, dima, val, work, lwork, &
    !      & rwrk, info)

    if ( info /= 0 ) then
      write (6, *) 'error: zheev error in guess_mix.'
      stop
    end if

  else if ( iwfnty == 2 ) then
    call zheevx ('v', 'a', 'u', dima, amat(1), dima, ddum, ddum, &
         & idum, idum, abstol, inum, val(1), avec(1), dima, work, &
         & lwork, rwrk, iwork, ifail, info)

    ! call zheev ('v', 'u', dima, amat, dima, val, work, lwork, &
    !      & rwrk, info)

    if ( info /= 0 ) then
      write (6, *) 'error: zheev error in guess_mix.'
      stop
    end if

    call zheevx ('v', 'a', 'u', dima, amat(dima*dima+1), dima, ddum, ddum, &
         & idum, idum, abstol, inum, val(dima+1), avec(dima*dima+1), dima, &
         & work, lwork, rwrk, iwork, ifail, info)

    amat(1:adim) = avec(1:adim)

    ! call zheev ('v', 'u', dima, amat(dima*dima+1), dima, val(dima+1), &
    !      & work, lwork, rwrk, info)

    if ( info /= 0 ) then
      write (6, *) 'error: zheev error in guess_mix.'
      stop
    end if 

  else if ( iwfnty == 3 ) then
    call zheevx ('v', 'a', 'u', 2*dima, amat, 2*dima, ddum, ddum, &
         & idum, idum, abstol, inum, val, avec, 2*dima, work, &
         & lwork, rwrk, iwork, ifail, info)

    amat(1:adim) = avec(1:adim)

    ! call zheev ('v', 'u', 2*dima, amat, 2*dima, val, work, lwork, &
    !      & rwrk, info)

    if ( info /= 0 ) then
      write (6, *) 'error: zheev error in guess_mix.'
      stop
    end if
  end if


  ! +---------------------------------------------+
  ! |  construct unitary matrix X = exp(i*fac*A)  |
  ! +---------------------------------------------+

  ! build diag(xmat):
  !   diag(xmat) = exp(i * fac * diag(amat))

  xmat = z0

  if ( iwfnty == 1 ) then
    do k = 1, dima
      ind = (k-1)*dima + k
      xmat(ind) = exp (zi * fac * val(k))
    end do

  else if ( iwfnty == 2 ) then
    do k = 1, dima
      ind = (k-1)*dima + k
      xmat(ind) = exp (zi * fac * val(k))
    end do

    do k = 1, dima
      ind = (k-1)*dima + k + dima*dima
      xmat(ind) = exp (zi * fac * val(k))
    end do

  else if ( iwfnty == 3 ) then
    do k = 1, 2*dima
      ind = (k-1)*2*dima + k
      xmat(ind) = exp (zi * fac * val(k))
    end do
  end if
 

  ! finish building xmat:
  !   xmat = vec . diag(xmat) . vec!

  if ( iwfnty == 1 ) then
    call zgemm ('n', 'n', dima, dima, dima, z1, amat, dima, &
         & xmat, dima, z0, scr, dima)
    call zgemm ('n', 'c', dima, dima, dima, z1, scr, dima, &
         & amat, dima, z0, xmat, dima)

  else if ( iwfnty == 2 ) then
    call zgemm ('n', 'n', dima, dima, dima, z1, amat, dima, &
         & xmat, dima, z0, scr, dima)
    call zgemm ('n', 'c', dima, dima, dima, z1, scr, dima, &
         & amat, dima, z0, xmat, dima)

    call zgemm ('n', 'n', dima, dima, dima, z1, amat(dima*dima+1), dima, &
         & xmat(dima*dima+1), dima, z0, scr, dima)
    call zgemm ('n', 'c', dima, dima, dima, z1, scr, dima, &
         & amat(dima*dima+1), dima, z0, xmat(dima*dima+1), dima)

  else if ( iwfnty == 3 ) then
    call zgemm ('n', 'n', 2*dima, 2*dima, 2*dima, z1, amat, 2*dima, &
         & xmat, 2*dima, z0, scr, 2*dima)
    call zgemm ('n', 'c', 2*dima, 2*dima, 2*dima, z1, scr, 2*dima, &
         & amat, 2*dima, z0, xmat, 2*dima)
  end if


  ! compute  dmat = dmat . xmat

  if ( iwfnty == 1 ) then
    if ( igsnmx > 0 ) then
      ind = (nup-igsnmx)*norb + 1
    else
      ind = 1
    end if

    call zgemm ('n', 'n', norb, dima, dima, z1, dmat(ind), norb, &
         & xmat, dima, z0, scr, norb)

    dmat(ind:ind+norb*dima-1) = scr(1:norb*dima)

  else if ( iwfnty == 2 ) then
    if ( igsnmx > 0 ) then
      ind = (nup-igsnmx)*norb + 1
    else
      ind = 1
    end if

    call zgemm ('n', 'n', norb, dima, dima, z1, dmat(ind), norb, &
         & xmat, dima, z0, scr, norb)

    dmat(ind:ind+norb*dima-1) = scr(1:norb*dima)

    if ( igsnmx > 0 ) then
      ind = nosq + (ndn-igsnmx)*norb + 1
    else
      ind = nosq + 1
    end if

    call zgemm ('n', 'n', norb, dima, dima, z1, dmat(ind), norb, &
         & xmat(dima*dima+1), dima, z0, scr, norb)

    dmat(ind:ind+norb*dima-1) = scr(1:norb*dima)

  else if ( iwfnty == 3 ) then
    if ( igsnmx > 0 ) then
      ind = (nup+ndn-2*igsnmx)*2*norb + 1
    else
      ind = 1
    end if

    call zgemm ('n', 'n', 2*norb, 2*dima, 2*dima, z1, dmat(ind), 2*norb, &
         & xmat, 2*dima, z0, scr, 2*norb)

    dmat(ind:ind+4*norb*dima-1) = scr(1:4*norb*dima)
  end if


  ! deallocate memory

  deallocate (amat, avec, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'amat', 'guess_mix')

  deallocate (xmat, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'xmat', 'guess_mix')

  deallocate (scr, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scr', 'guess_mix')

  deallocate (iwork, ifail, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'iwork', 'guess_mix')

  deallocate (work, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'work', 'guess_mix')

  deallocate (rwrk, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rwrk', 'guess_mix')

  deallocate (val, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'val', 'guess_mix')


  return
end subroutine guess_mix



subroutine gsproj2 (norb, nbas, xdim, dmat, smat, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | gsproj2  --  CAJH, 01.2013                                     |
! |                                                                |
! |                                                                |
! | Project a given matrix D of orbital coefficients into the      |
! | current std AO basis [for GHF type determinants].              |
! |                                                                |
! | That is, if D is the matrix of orbital coefficients provided   |
! | on input, this subroutine returns D' such that                 |
! |                                                                |
! |        =>   D!. S . D  =  A                                    |
! |        =>   D'! . S . D'  =  I,                                |
! |                                                                |
! | where D' = D . A^(-1/2).                                       |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  amat is not positive definite                   |
! |         =  2,  diagonalization of amat failed                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb - number of orbitals
  !   nbas - number of basis functions
  !   xdim - leading dimension of dmat
  !   dmat - matrix of orbital coefficients [ updated here ]
  !   smat - overlap matrix ( std AO basis )

  integer, intent(in) :: norb, nbas, xdim
  integer, intent(out) :: iflg

  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat
  complex(kind=dp), dimension(xdim), intent(inout) :: dmat


  ! other variables

  integer :: lwrk, nb, info, k
  integer :: istatus

  integer :: idum, inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  integer, dimension(:), allocatable :: iwork, ifail

  real(kind=dp), dimension(:), allocatable :: aval, rwrk

  complex(kind=dp), dimension(:), allocatable :: wrk
  complex(kind=dp), dimension(:,:), allocatable :: smatX

  complex(kind=dp), dimension(:,:), allocatable :: amat, avec
  complex(kind=dp), dimension(:,:), allocatable :: scra


  ! constants

  !   eps - estimate of machine precision

  real(kind=dp) :: eps
  parameter ( eps = 1.0-15_dp )


  ! external functions

  integer, external :: ilaenv
  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  iflg = 0

  ! build smatX (overlap matrix in spin-orbital basis)

  allocate (smatX(2*nbas,2*nbas), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'smatX', 'gsproj2')

  smatX(1:nbas,1:nbas)               = smat(1:nbas,1:nbas)
  smatX(nbas+1:2*nbas,1:nbas)        = z0
  smatX(1:nbas,nbas+1:2*nbas)        = z0
  smatX(nbas+1:2*nbas,nbas+1:2*nbas) = smat(1:nbas,1:nbas)

  ! allocate amat, avec, aval

  allocate (amat(2*norb,2*norb), avec(2*norb,2*norb), &
          & aval(2*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'amat', 'gsproj2')

  ! allocate space for scratch arrays

  allocate (scra(2*nbas,2*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'gsproj2')

  ! determine lwrk
  ! allocate wrk, rwrk

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', 2*norb, -1, -1, -1)
  lwrk = max (1, (nb+1)*2*norb)

  allocate (iwork(10*norb), ifail(2*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'gsproj2')

  ! allocate (wrk(lwrk), rwrk(6*norb-2), stat=istatus)
  allocate (wrk(lwrk), rwrk(14*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wrk', 'gsproj2')


  ! compute D! . S . D = A
  ! store in amat

  call zgemm ('n', 'n', 2*nbas, 2*norb, 2*nbas, z1, smatX, &
       & 2*nbas, dmat, 2*nbas, z0, scra, 2*nbas)
  call zgemm ('c', 'n', 2*norb, 2*norb, 2*nbas, z1, dmat, &
       & 2*nbas, scra, 2*nbas, z0, amat, 2*norb)

  ! diagonalize A matrix
  ! set to negative to recover largest eigenvalues on top

  ! avec(1:2*norb,1:2*norb) = -amat(1:2*norb,1:2*norb)
  amat(1:2*norb,1:2*norb) = -amat(1:2*norb,1:2*norb)

  call zheevx ('v', 'a', 'u', 2*norb, amat, 2*norb, ddum, ddum, &
       & idum, idum, abstol, inum, aval, avec, 2*norb, wrk, &
       & lwrk, rwrk, iwork, ifail, info)

  ! call zheev ('v', 'u', 2*norb, avec, 2*norb, aval, wrk, lwrk, &
  !      & rwrk, info)

  if ( info /= 0 ) then
    iflg = 2
    return
  end if

  ! check that amat is positive definite
  ! build a^(-1/2)

  amat(1:2*norb,1:2*norb) = z0

  do k = 1, 2*norb
    if ( aval(k) < eps ) then
      iflg = 1
      return
    else
      amat(k,k) = d1/sqrt(-aval(k))
    end if
  end do

  ! build A^(-1/2) =  vecA . a^(-1/2) . vecA!

  call zgemm ('n', 'c', 2*norb, 2*norb, 2*norb, z1, amat, 2*norb, &
       & avec, 2*norb, z0, scra, 2*norb)
  call zgemm ('n', 'n', 2*norb, 2*norb, 2*norb, z1, avec, 2*norb, &
       & scra, 2*norb, z0, amat, 2*norb)

  ! form D = D . A^(-1/2)

  call zgemm ('n', 'n', 2*nbas, 2*norb, 2*norb, z1, dmat, &
       & 2*nbas, amat, 2*norb, z0, scra, 2*nbas)

  dmat(1:4*nbas*norb) = &
       & reshape (scra(1:2*nbas,1:2*norb), (/ 4*nbas*norb /))


  ! deallocate memory

  deallocate (smatX, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'smatX', 'gsproj2')

  deallocate (amat, avec, aval, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'amat', 'gsproj2')

  deallocate (scra, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scr?', 'gsproj2')

  deallocate (iwork, ifail, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'iwork', 'gsproj2')

  deallocate (wrk, rwrk, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'wrk', 'gsproj2')


  return
end subroutine gsproj2



subroutine gsproj1 (norb, nbas, xdim, dmat, smat, iflg)

! +----------------------------------------------------------------+
! |                                                                |
! | gsproj1  --  CAJH, 01.2013                                     |
! |                                                                |
! |                                                                |
! | Project a given matrix D of orbital coefficients into the      |
! | current std AO basis [for RHF, UHF type determinants].         |
! |                                                                |
! | That is, if D is the matrix of orbital coefficients provided   |
! | on input, this subroutine returns D' such that                 |
! |                                                                |
! |        =>   D!. S . D  =  A                                    |
! |        =>   D'! . S . D'  =  I,                                |
! |                                                                |
! | where D' = D . A^(-1/2).                                       |
! |                                                                |
! | Error codes:                                                   |
! |                                                                |
! |   iflg  =  0,  successful completion                           |
! |         =  1,  amat is not positive definite                   |
! |         =  2,  diagonalization of amat failed                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb - number of orbitals
  !   nbas - number of basis functions
  !   xdim - leading dimension of dmat
  !   dmat - matrix of orbital coefficients [ updated here ]
  !   smat - overlap matrix ( std AO basis )

  integer, intent(in) :: norb, nbas, xdim
  integer, intent(out) :: iflg

  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat
  complex(kind=dp), dimension(xdim), intent(inout) :: dmat


  ! other variables

  integer :: lwrk, nb, info, k
  integer :: istatus

  integer :: idum, inum
  real(kind=dp) :: ddum, abstol
  ! parameter ( abstol = 2.0e0_dp*0.149166814624004135e-153_dp )

  integer, dimension(:), allocatable :: iwork, ifail

  real(kind=dp), dimension(:), allocatable :: aval, rwrk

  complex(kind=dp), dimension(:), allocatable :: wrk

  complex(kind=dp), dimension(:,:), allocatable :: amat, avec
  complex(kind=dp), dimension(:,:), allocatable :: scra


  ! constants

  !   eps - estimate of machine precision

  real(kind=dp) :: eps
  parameter ( eps = 1.0-15_dp )


  ! external functions

  integer, external :: ilaenv
  real(kind=dp), external :: dlamch


  idum = 0
  ddum = d0

  abstol = dlamch('s')


  iflg = 0

  ! allocate amat, avec, aval

  allocate (amat(norb,norb), avec(norb,norb), &
          & aval(norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'amat', 'gsproj1')

  ! allocate space for scratch arrays

  allocate (scra(nbas,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scr?', 'gsproj1')

  ! determine lwrk
  ! allocate wrk, rwrk

  ! nb = 2
  nb = ilaenv (1, 'zhetrd', 'u', norb, -1, -1, -1)
  lwrk = max (1, (nb+1)*norb)

  allocate (iwork(5*norb), ifail(norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'iwork', 'gsproj1')

  ! allocate (wrk(lwrk), rwrk(3*norb-2), stat=istatus)
  allocate (wrk(lwrk), rwrk(7*norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wrk', 'gsproj1')


  ! compute D! . S . D = A
  ! store in amat

  call zgemm ('n', 'n', nbas, norb, nbas, z1, smat, &
       & nbas, dmat, nbas, z0, scra, nbas)
  call zgemm ('c', 'n', norb, norb, nbas, z1, dmat, &
       & nbas, scra, nbas, z0, amat, norb)

  ! diagonalize A matrix
  ! set to negative to recover largest eigenvalues on top

  ! avec(1:norb,1:norb) = -amat(1:norb,1:norb)
  amat(1:norb,1:norb) = -amat(1:norb,1:norb)

  call zheevx ('v', 'a', 'u', norb, amat, norb, ddum, ddum, &
       & idum, idum, abstol, inum, aval, avec, norb, wrk, &
       & lwrk, rwrk, iwork, ifail, info)

  ! call zheev ('v', 'u', norb, avec, norb, aval, wrk, lwrk, &
  !      & rwrk, info)

  if ( info /= 0 ) then
    iflg = 2
    return
  end if

  ! check that amat is positive definite
  ! build a^(-1/2)

  amat(1:norb,1:norb) = z0

  do k = 1, norb
    if ( aval(k) < eps ) then
      iflg = 1
      return
    else
      amat(k,k) = d1/sqrt(-aval(k))
    end if
  end do

  ! build A^(-1/2) =  vecA . a^(-1/2) . vecA!

  call zgemm ('n', 'c', norb, norb, norb, z1, amat, norb, &
       & avec, norb, z0, scra, norb)
  call zgemm ('n', 'n', norb, norb, norb, z1, avec, norb, &
       & scra, norb, z0, amat, norb)

  ! form D = D . A^(-1/2)

  call zgemm ('n', 'n', nbas, norb, norb, z1, dmat, &
       & nbas, amat, norb, z0, scra, nbas)

  dmat(1:nbas*norb) = &
       & reshape (scra(1:nbas,1:norb), (/ nbas*norb /))


  ! deallocate memory

  deallocate (amat, avec, aval, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'amat', 'gsproj1')

  deallocate (scra, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scr?', 'gsproj1')

  deallocate (iwork, ifail, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'iwork', 'gsproj1')

  deallocate (wrk, rwrk, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'wrk', 'gsproj1')


  return
end subroutine gsproj1


end module iguess


