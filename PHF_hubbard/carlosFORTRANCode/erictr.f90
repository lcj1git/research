

module erictr

  use constants
  use i2sint
  use util, only : error_alloc
  use purcart
  use linalg

  implicit none
  private
  save

  public :: setup_erictr, shutdown_erictr
  public :: erictr_ghf, erictr_uhf, erictr_rhf

! +----------------------------------------------------------------+
! |                                                                |
! | erictr                                                         |
! |                                                                |
! |                                                                |
! | A collection of subroutines that perform the contraction of    |
! | (transition) density matrices with electron-repulsion          |
! | integrals (ERIs).                                              |
! |                                                                |
! | NOTE:                                                          |
! |   In order to use this module,                                 |
! |     * setup_erictr *                                           |
! |   should be invoked before any calls to erictr_***.            |
! |   When the module is no longer needed, one may call            |
! |     * shutdown_erictr *                                        |
! |   to deallocate all scratch arrays.                            |
! |                                                                |
! +----------------------------------------------------------------+

! +----------------------------------------------------------------+
! |                                                                |
! | Some comments about ERIs                                       |
! |                                                                |
! | The subroutines in this module expect two-electron integrals   |
! | stored in Mulliken notation                                    |
! |                                                                |
! |   ( i j | k l )  =  < i k | j l >                              |
! |      Mulliken         Dirac                                    |
! |                                                                |
! | This is the most common form in which quantum-chemistry        |
! | packages build two-electron integrals. Integrals should be     |
! | stored as an array of (i2s) type, defined in i2sint:           |
! |                                                                |
! |   type(i2s)   integer*2 - ii, ij, ik, il                       |
! |                  real*8 - int                                  |
! |                                                                |
! | This allows efficient (sparse) storage of ERIs. No particular  |
! | ordering of two-electron integrals is assumed in this module.  |
! |                                                                |
! | Real two-electron integrals have the following symmetry        |
! | properties:                                                    |
! |                                                                |
! |   ( i j | k l )  =  ( k l | i j )   [ interchange 1-2 ]        |
! |                                                                |
! |   ( i j | k l )  =  ( j i | k l )   [ real integrals ]         |
! |                  =  ( i j | l k )                              |
! |                  =  ( j i | l k )                              |
! |   ( k l | i j )  =  ( k l | j i )                              |
! |                  =  ( l k | i j )                              |
! |                  =  ( l k | j i )                              |
! |                                                                |
! | Subroutines in this module understand this symmetry and use    |
! | it accordingly. Therefore, the array i2sv of ERIs provided     |
! | on input must contain a single copy of symmetry-related ints.  |
! | For instance,                                                  |
! |                                                                |
! |   if i2sv contains     ( 1 2 | 1 2 ),                          |
! |   it must not include  ( 1 2 | 2 1 ), ( 2 1 | 1 2 ),           |
! |                        ( 2 1 | 2 1 ).                          |
! |                                                                |
! | Some quantum chemistry packages build integrals including a    |
! | 'symmetry factor'. Subroutines in this module can handle that  |
! | by turning lsf = .true. on input. In particular, one can       |
! | classify ERIs according to their number of symmetry-related    |
! | partners:                                                      |
! |                                                                |
! |   ( i i | i i )  -  1                                          |
! |   ( i i | k k )  -  2                                          |
! |   ( i j | k k )  -  4                                          |
! |   ( i i | k l )  -  4                                          |
! |   ( i j | i j )  -  4                                          |
! |   ( i j | j i )  -  4                                          |
! |   ( i j | k l )  -  8                                          |
! |                                                                |
! | When lsf = .true., provided ERIs include that symmetry factor: |
! | an integral ( i j | k l ) should be divided by a factor of 8.  |
! |                                                                |
! | Lastly, this module expects integrals to be provided over      |
! | Cartesian basis functions. If input density matrices are in    |
! | a pure basis, the purcart module is used to perform pure to    |
! | Cartesian transformations.                                     |
! |                                                                |
! +----------------------------------------------------------------+


  ! scratch arrays

  integer :: nblk
  parameter ( nblk = 8192 )

  type (i2s), dimension(:), allocatable :: i2sblk

  complex(kind=dp), dimension(:,:), allocatable :: scrt

  complex(kind=dp), dimension(:,:,:), allocatable :: rho_ao_sp, gam_ao_sp
  complex(kind=dp), dimension(:,:,:), allocatable :: rho_ao2_sp, gam_ao2_sp
  complex(kind=dp), dimension(:,:,:), allocatable :: rhos_ao2_sp

  !$omp  threadprivate(scrt, rho_ao_sp, gam_ao_sp, rhos_ao2_sp, &
  !$omp&               rho_ao2_sp, gam_ao2_sp, i2sblk)


contains


subroutine setup_erictr (iwfnty, norb, nbas, nbct)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_erictr  --  CAJH, 12.2012                                |
! |                                                                |
! |                                                                |
! | Allocate memory for all scratch arrays used in erictr.         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iwfnty - type of wavefunction to use
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions

  integer, intent(in) :: iwfnty, norb, nbas, nbct


  ! other variables

  integer :: istatus


  ! scratch space for blocks of ERIs

  allocate (i2sblk(nblk), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'i2sblk', 'setup_erictr')


  ! scratch space for ort AO - std AO basis transformations

  allocate (scrt(nbas,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scrt', 'setup_erictr')


  ! allocate space for spin-blocks of rho_ao

  if ( iwfnty == 1 ) then
    allocate (rho_ao_sp(nbas,nbas,1), rho_ao2_sp(nbct,nbct,1), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (rho_ao_sp(nbas,nbas,2), rho_ao2_sp(nbct,nbct,2), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (rho_ao_sp(nbas,nbas,4), rho_ao2_sp(nbct,nbct,4), stat=istatus)
  end if
  if ( istatus /= 0 ) call error_alloc (1, 'rho_ao?_sp', 'setup_erictr')


  ! allocate space for spin-blocks of gam_ao

  if ( iwfnty == 1 ) then
    allocate (gam_ao_sp(nbas,nbas,1), gam_ao2_sp(nbct,nbct,1), stat=istatus)
  else if ( iwfnty == 2 ) then
    allocate (gam_ao_sp(nbas,nbas,2), gam_ao2_sp(nbct,nbct,2), stat=istatus)
  else if ( iwfnty == 3 ) then
    allocate (gam_ao_sp(nbas,nbas,4), gam_ao2_sp(nbct,nbct,4), stat=istatus)
  end if
  if ( istatus /= 0 ) call error_alloc (1, 'gam_ao?_sp', 'setup_erictr')

  ! other scratch blocks

  if ( iwfnty == 1 ) then
    allocate (rhos_ao2_sp(nbct,nbct,1), stat=istatus)
  else if ( iwfnty == 2 .or. iwfnty == 3 ) then
    allocate (rhos_ao2_sp(nbct,nbct,2), stat=istatus)
  end if
  if ( istatus /= 0 ) call error_alloc (1, 'rhos_ao2_sp', 'setup_erictr')


  return
end subroutine setup_erictr



subroutine shutdown_erictr

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_erictr  --  CAJH, 12.2012                             |
! |                                                                |
! |                                                                |
! | Deallocate memory for all scratch arrays used in erictr.       |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: istatus


  ! scratch space for blocks of ERIs

  deallocate (i2sblk, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'i2sblk', 'setup_erictr')

  ! scratch space for ort AO - std AO basis transformations

  deallocate (scrt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scrt', 'shutdown_erictr')

  ! deallocate space for spin-blocks of rho_ao

  deallocate (rho_ao_sp, rho_ao2_sp, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rho_ao?_sp', 'setup_erictr')

  ! deallocate space for spin-blocks of gam_ao

  deallocate (gam_ao_sp, gam_ao2_sp, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'gam_ao?_sp', 'setup_erictr')

  ! other scratch blocks

  deallocate (rhos_ao2_sp, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'rhos_ao2_sp', 'setup_erictr')


  return
end subroutine shutdown_erictr



subroutine erictr_ghf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | erictr_ghf  --  CAJH, 12.2012                                  |
! |                                                                |
! |                                                                |
! | Contract a GHF-type (transition) density matrix with anti-     |
! | symmetrized two-electron integrals. That is, form              |
! |                                                                |
! |   gam(i,k)  =  sum_lj  < ij || kl > * rho(l,j).                |
! |                                                                |
! |                                                                |
! | - The input density matrix should be in the ort AO basis.      |
! |                                                                |
! | - Electron repulsion integrals should be provided in the       |
! |   array i2sv with the following features:                      |
! |                                                                |
! |     o  stored non-antisymmetrized, in Mulliken notation        |
! |     o  integrals are stored in the standard AO basis           |
! |     o  only one copy of symmetry-related integrals is          |
! |        included in i2sv                                        |
! |                                                                |
! | - The output matrix, gam, is returned in the ort AO basis.     |
! |                                                                |
! |                                                                |
! | NOTE: The transformation matrix xmat [ = S^(-1/2) ] is used to |
! |   perform basis transformations inside this routine.           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb - number of orbitals
  !   nbas - number of basis functions
  !   nbct - number of Cartesian basis functions
  !   rho  - density matrix to contract
  !   gam  - contracted density matrix [ out ]
  !   xmat - transformation matrix [ = S^(-1/2) ]
  !   i2sv - vector of 2-electron integrals
  !   ni2s - number of 2-electron integrals stored in i2sv
  !   lsf  - whether stored integrals include symmetry factor

  integer, intent(in) :: norb, nbas, nbct, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(2*norb,2*norb), intent(inout) :: rho
  complex(kind=dp), dimension(2*norb,2*norb), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: j, k


  ! transform to std AO basis

  call aobstf (1, norb, nbas, xmat, rho(1,1), 2*norb, &
       & rho_ao_sp(1,1,1), nbas)
  call aobstf (1, norb, nbas, xmat, rho(norb+1,norb+1), 2*norb, &
       & rho_ao_sp(1,1,2), nbas)
  call aobstf (1, norb, nbas, xmat, rho(1,norb+1), 2*norb, &
       & rho_ao_sp(1,1,3), nbas)
  call aobstf (1, norb, nbas, xmat, rho(norb+1,1), 2*norb, &
       & rho_ao_sp(1,1,4), nbas)

  ! transform to Cartesian functions

  if ( nbct > nbas ) then
    do k = 1, 4
      call pur2cart (1, .true., nbas, nbct, rho_ao_sp(1,1,k), rho_ao2_sp(1,1,k))
    end do
  else
    rho_ao2_sp(1:nbct,1:nbct,1:4) = rho_ao_sp(1:nbas,1:nbas,1:4)
  end if

  ! symmetrize uu, dd blocks

  rhos_ao2_sp(1:nbct,1:nbct,1:2) = rho_ao2_sp(1:nbct,1:nbct,1:2)

  do k = 1, 2
    call zmat_symm (2, 'u', nbct, rhos_ao2_sp(1,1,k))
  end do

  ! clear gam_ao

  gam_ao2_sp = z0

  ! Coulomb contraction (only UT)

  call coul_ctr (nbct, 2, rhos_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  do j = 1, nbct
    do k = 1, j
      gam_ao2_sp(k,j,1) = gam_ao2_sp(k,j,1) + gam_ao2_sp(k,j,2)
      gam_ao2_sp(k,j,2) = gam_ao2_sp(k,j,1)
    end do
  end do

  ! expand gam to square

  do k = 1, 2
    call zmat_square (2, 'u', nbct, gam_ao2_sp(1,1,k))
  end do

  ! exchange contraction

  call exch_ctr (nbct, 4, rho_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  ! transform back to pure functions

  if ( nbct > nbas ) then
    do k = 1, 4
      call pur2cart (2, .false., nbas, nbct, gam_ao_sp(1,1,k), gam_ao2_sp(1,1,k))
    end do
  else
    gam_ao_sp(1:nbas,1:nbas,1:4) = gam_ao2_sp(1:nbct,1:nbct,1:4)
  end if

  ! transform back to ort AO basis

  call aobstf (2, norb, nbas, xmat, gam(1,1), 2*norb, &
       & gam_ao_sp(1,1,1), nbas)
  call aobstf (2, norb, nbas, xmat, gam(norb+1,norb+1), 2*norb, &
       & gam_ao_sp(1,1,2), nbas)
  call aobstf (2, norb, nbas, xmat, gam(1,norb+1), 2*norb, &
       & gam_ao_sp(1,1,3), nbas)
  call aobstf (2, norb, nbas, xmat, gam(norb+1,1), 2*norb, &
       & gam_ao_sp(1,1,4), nbas)


  return
end subroutine erictr_ghf



subroutine erictr_uhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | erictr_uhf  --  CAJH, 12.2012                                  |
! |                                                                |
! |                                                                |
! | ( UHF version of erictr_ghf. )                                 |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb - number of orbitals
  !   nbas - number of basis functions
  !   nbct - number of Cartesian basis functions
  !   rho  - density matrix to contract
  !   gam  - contracted density matrix [ out ]
  !   xmat - transformation matrix [ = S^(-1/2) ]
  !   i2sv - vector of 2-electron integrals
  !   ni2s - number of 2-electron integrals stored in i2sv
  !   lsf  - whether stored integrals include symmetry factor

  integer, intent(in) :: norb, nbas, nbct, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,2*norb), intent(inout) :: rho
  complex(kind=dp), dimension(norb,2*norb), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: j, k


  ! transform to std AO basis

  call aobstf (1, norb, nbas, xmat, rho(1,1), norb, &
       & rho_ao_sp(1,1,1), nbas)
  call aobstf (1, norb, nbas, xmat, rho(1,norb+1), norb, &
       & rho_ao_sp(1,1,2), nbas)

  ! transform to Cartesian functions

  if ( nbct > nbas ) then
    do k = 1, 2
      call pur2cart (1, .true., nbas, nbct, rho_ao_sp(1,1,k), rho_ao2_sp(1,1,k))
    end do
  else
    rho_ao2_sp(1:nbct,1:nbct,1:2) = rho_ao_sp(1:nbas,1:nbas,1:2)
  end if

  ! symmetrize uu, dd blocks

  rhos_ao2_sp = rho_ao2_sp

  do k = 1, 2
    call zmat_symm (2, 'u', nbct, rhos_ao2_sp(1,1,k))
  end do

  ! clear gam_ao

  gam_ao2_sp = z0

  ! Coulomb contraction (only UT)

  call coul_ctr (nbct, 2, rhos_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  do j = 1, nbct
    do k = 1, j
      gam_ao2_sp(k,j,1) = gam_ao2_sp(k,j,1) + gam_ao2_sp(k,j,2)
      gam_ao2_sp(k,j,2) = gam_ao2_sp(k,j,1)
    end do
  end do

  ! expand gam to square

  do k = 1, 2
    call zmat_square (2, 'u', nbct, gam_ao2_sp(1,1,k))
  end do

  ! exchange contraction

  call exch_ctr (nbct, 2, rho_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  ! transform back to pure functions

  if ( nbct > nbas ) then
    do k = 1, 2
      call pur2cart (2, .false., nbas, nbct, gam_ao_sp(1,1,k), gam_ao2_sp(1,1,k))
    end do
  else
    gam_ao_sp(1:nbas,1:nbas,1:2) = gam_ao2_sp(1:nbct,1:nbct,1:2)
  end if

  ! transform back to ort AO basis

  call aobstf (2, norb, nbas, xmat, gam(1,1), norb, &
       & gam_ao_sp(1,1,1), nbas)
  call aobstf (2, norb, nbas, xmat, gam(1,norb+1), norb, &
       & gam_ao_sp(1,1,2), nbas)


  return
end subroutine erictr_uhf



subroutine erictr_rhf (norb, nbas, nbct, rho, gam, xmat, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | erictr_rhf  --  CAJH, 12.2012                                  |
! |                                                                |
! |                                                                |
! | ( RHF version of erictr_ghf. )                                 |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   norb - number of orbitals
  !   nbas - number of basis functions
  !   nbct - number of Cartesian basis functions
  !   rho  - density matrix to contract
  !   gam  - contracted density matrix [ out ]
  !   xmat - transformation matrix [ = S^(-1/2) ]
  !   i2sv - vector of 2-electron integrals
  !   ni2s - number of 2-electron integrals stored in i2sv
  !   lsf  - whether stored integrals include symmetry factor

  integer, intent(in) :: norb, nbas, nbct, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(norb,norb), intent(inout) :: rho
  complex(kind=dp), dimension(norb,norb), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! transform to std AO basis

  call aobstf (1, norb, nbas, xmat, rho, norb, &
       & rho_ao_sp, nbas)

  ! transform to Cartesian functions

  if ( nbct > nbas ) then
    call pur2cart (1, .true., nbas, nbct, rho_ao_sp, rho_ao2_sp)
  else
    rho_ao2_sp(1:nbct,1:nbct,1) = rho_ao_sp(1:nbas,1:nbas,1)
  end if

  ! symmetrize uu, dd blocks

  rhos_ao2_sp = rho_ao2_sp

  call zmat_symm (2, 'u', nbct, rhos_ao2_sp)

  ! clear gam_ao

  gam_ao2_sp = z0

  ! Coulomb contraction (only UT)

  call coul_ctr (nbct, 1, rhos_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  gam_ao2_sp(1:nbct,1:nbct,1) = d2 * gam_ao2_sp(1:nbct,1:nbct,1)

  ! expand gam to square

  call zmat_square (2, 'u', nbct, gam_ao2_sp)

  ! exchange contraction

  call exch_ctr (nbct, 1, rho_ao2_sp, gam_ao2_sp, i2sv, ni2s, lsf)

  ! transform back to pure functions

  if ( nbct > nbas ) then
    call pur2cart (2, .false., nbas, nbct, gam_ao_sp, gam_ao2_sp)
  else
    gam_ao_sp(1:nbas,1:nbas,1) = gam_ao2_sp(1:nbct,1:nbct,1)
  end if

  ! transform back to ort AO basis

  call aobstf (2, norb, nbas, xmat, gam, norb, &
       & gam_ao_sp, nbas)


  return
end subroutine erictr_rhf



subroutine coulomb_ctr (nbct, rho, gam, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | coulomb_ctr  --  CAJH, 12.2012                                 |
! |                                                                |
! |                                                                |
! | Update gam with the Coulomb contraction of two-el integrals    |
! | with (transition) density matrix rho. That is,                 |
! |                                                                |
! |   gam(i,j)  =  gam(i,j)  +  sum_kl ( i j | k l ) * rho(l,k)    |
! |                                                                |
! |                                                                |
! | IMPORTANT NOTES :                                              |
! |                                                                |
! | - rho is assumed symmetric;                                    |
! |     only the upper triangle is referenced                      |
! |                                                                |
! | - contraction of rho with ints produces a symmetric gam;       |
! |     only the upper triangle of gam is updated                  |
! |                                                                |
! |                                                                |
! | WARNING : This subroutine assumes that only one copy of        |
! |   symmetry-related integrals is included in i2sv !             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nbct - number of Cartesian basis functions
  !   rho  - density matrix to contract
  !   gam  - contracted density matrix [ updated on output ]
  !   i2sv - vector of 2-electron integrals
  !   ni2s - number of 2-electron integrals stored in i2sv
  !   lsf  - whether stored integrals include symmetry factor

  integer, intent(in) :: nbct, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbct,nbct), intent(in) :: rho
  complex(kind=dp), dimension(nbct,nbct), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: k
  logical :: l1ij, l1kl, l2ij, l2kl

  integer*2 :: ji, jj, jk, jl


  ! constants

  real(kind=dp) :: f2, f12, f14
  parameter ( f2 = 2.0_dp, f12 = 0.50_dp, f14 = 0.25_dp )


  if ( lsf ) then

    do k = 1, ni2s
      ji = i2sv(k)%ii
      jj = i2sv(k)%ij
      jk = i2sv(k)%ik
      jl = i2sv(k)%il

      l1ij = ji == jj
      l1kl = jk == jl

      l2ij = ji > jj
      l2kl = jk > jl

      if ( l1ij .and. l1kl .and. ji == jk ) then   ! ( i i | i i )

        gam(ji,ji) = gam(ji,ji) + rho(ji,ji) * i2sv(k)%int

      else if ( l1ij .and. l1kl ) then             ! ( i i | k k )

        gam(ji,ji) = gam(ji,ji) + f12 * rho(jk,jk) * i2sv(k)%int
        gam(jk,jk) = gam(jk,jk) + f12 * rho(ji,ji) * i2sv(k)%int

      else if ( l1kl ) then                        ! ( i j | k k )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f14 * rho(jk,jk) * i2sv(k)%int
          gam(jk,jk) = gam(jk,jk) + f12 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f14 * rho(jk,jk) * i2sv(k)%int
          gam(jk,jk) = gam(jk,jk) + f12 * rho(ji,jj) * i2sv(k)%int
        end if

      else if ( l1ij ) then                        ! ( i i | k l )

        if ( l2kl ) then
          gam(jl,jk) = gam(jl,jk) + f14 * rho(ji,ji) * i2sv(k)%int
          gam(ji,ji) = gam(ji,ji) + f12 * rho(jl,jk) * i2sv(k)%int
        else
          gam(jk,jl) = gam(jk,jl) + f14 * rho(ji,ji) * i2sv(k)%int
          gam(ji,ji) = gam(ji,ji) + f12 * rho(jk,jl) * i2sv(k)%int
        end if

      else if ( ji == jk .and. jj == jl ) then     ! ( i j | i j )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f12 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f12 * rho(ji,jj) * i2sv(k)%int
        end if

      else if ( ji == jl .and. jj == jk ) then     ! ( i j | j i )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f12 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f12 * rho(ji,jj) * i2sv(k)%int
        end if

      else                                         ! ( i j | k l )

        if ( l2ij .and. l2kl ) then
          gam(jj,ji) = gam(jj,ji) + f14 * rho(jl,jk) * i2sv(k)%int
          gam(jl,jk) = gam(jl,jk) + f14 * rho(jj,ji) * i2sv(k)%int
        else if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f14 * rho(jk,jl) * i2sv(k)%int
          gam(jk,jl) = gam(jk,jl) + f14 * rho(jj,ji) * i2sv(k)%int
        else if ( l2kl ) then
          gam(ji,jj) = gam(ji,jj) + f14 * rho(jl,jk) * i2sv(k)%int
          gam(jl,jk) = gam(jl,jk) + f14 * rho(ji,jj) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f14 * rho(jk,jl) * i2sv(k)%int
          gam(jk,jl) = gam(jk,jl) + f14 * rho(ji,jj) * i2sv(k)%int
        end if

      end if

    end do

  else

    do k = 1, ni2s
      ji = i2sv(k)%ii
      jj = i2sv(k)%ij
      jk = i2sv(k)%ik
      jl = i2sv(k)%il

      l1ij = ji == jj
      l1kl = jk == jl

      l2ij = ji > jj
      l2kl = jk > jl

      if ( l1ij .and. l1kl .and. ji == jk ) then   ! ( i i | i i )

        gam(ji,ji) = gam(ji,ji) + rho(ji,ji) * i2sv(k)%int

      else if ( l1ij .and. l1kl ) then             ! ( i i | k k )

        gam(ji,ji) = gam(ji,ji) + rho(jk,jk) * i2sv(k)%int
        gam(jk,jk) = gam(jk,jk) + rho(ji,ji) * i2sv(k)%int

      else if ( l1kl ) then                        ! ( i j | k k )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + rho(jk,jk) * i2sv(k)%int
          gam(jk,jk) = gam(jk,jk) + f2 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + rho(jk,jk) * i2sv(k)%int
          gam(jk,jk) = gam(jk,jk) + f2 * rho(ji,jj) * i2sv(k)%int
        end if

      else if ( l1ij ) then                        ! ( i i | k l )

        if ( l2kl ) then
          gam(jl,jk) = gam(jl,jk) + rho(ji,ji) * i2sv(k)%int
          gam(ji,ji) = gam(ji,ji) + f2 * rho(jl,jk) * i2sv(k)%int
        else
          gam(jk,jl) = gam(jk,jl) + rho(ji,ji) * i2sv(k)%int
          gam(ji,ji) = gam(ji,ji) + f2 * rho(jk,jl) * i2sv(k)%int
        end if

      else if ( ji == jk .and. jj == jl ) then     ! ( i j | i j )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f2 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f2 * rho(ji,jj) * i2sv(k)%int
        end if

      else if ( ji == jl .and. jj == jk ) then     ! ( i j | j i )

        if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f2 * rho(jj,ji) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f2 * rho(ji,jj) * i2sv(k)%int
        end if

      else                                         ! ( i j | k l )

        if ( l2ij .and. l2kl ) then
          gam(jj,ji) = gam(jj,ji) + f2 * rho(jl,jk) * i2sv(k)%int
          gam(jl,jk) = gam(jl,jk) + f2 * rho(jj,ji) * i2sv(k)%int
        else if ( l2ij ) then
          gam(jj,ji) = gam(jj,ji) + f2 * rho(jk,jl) * i2sv(k)%int
          gam(jk,jl) = gam(jk,jl) + f2 * rho(jj,ji) * i2sv(k)%int
        else if ( l2kl ) then
          gam(ji,jj) = gam(ji,jj) + f2 * rho(jl,jk) * i2sv(k)%int
          gam(jl,jk) = gam(jl,jk) + f2 * rho(ji,jj) * i2sv(k)%int
        else
          gam(ji,jj) = gam(ji,jj) + f2 * rho(jk,jl) * i2sv(k)%int
          gam(jk,jl) = gam(jk,jl) + f2 * rho(ji,jj) * i2sv(k)%int
        end if

      end if

    end do
  end if


  return
end subroutine coulomb_ctr



subroutine exchange_ctr (nbct, rho, gam, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | exchange_ctr  --  CAJH, 12.2012                                |
! |                                                                |
! |                                                                |
! | Update gam with the exchange contraction of two-el integrals   |
! | with (transition) density matrix rho. That is,                 |
! |                                                                |
! |   gam(i,l)  =  gam(i,l)  -  sum_kl ( i j | k l ) * rho(j,k)    |
! |                                                                |
! |                                                                |
! | WARNING : This subroutine assumes that only one copy of        |
! |   symmetry-related integrals is included in i2sv !             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   nbct - number of Cartesian basis functions
  !   rho  - density matrix to contract
  !   gam  - contracted density matrix [ updated on output ]
  !   i2sv - vector of 2-electron integrals
  !   ni2s - number of 2-electron integrals stored in i2sv
  !   lsf  - whether stored integrals include symmetry factor

  integer, intent(in) :: nbct, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbct,nbct), intent(in) :: rho
  complex(kind=dp), dimension(nbct,nbct), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: k
  logical :: l1ij, l1kl

  integer*2 :: ji, jj, jk, jl


  ! constants

  real(kind=dp) :: f12, f14, f18
  parameter ( f12 = 0.50_dp, f14 = 0.25_dp, f18 = 0.125_dp )


  if ( lsf ) then

    do k = 1, ni2s
      ji = i2sv(k)%ii
      jj = i2sv(k)%ij
      jk = i2sv(k)%ik
      jl = i2sv(k)%il

      l1ij = ji == jj
      l1kl = jk == jl

      if ( l1ij .and. l1kl .and. ji == jk ) then   ! ( i i | i i )

        gam(ji,ji) = gam(ji,ji) - rho(ji,ji) * i2sv(k)%int

      else if ( l1ij .and. l1kl ) then             ! ( i i | k k )

        gam(ji,jk) = gam(ji,jk) - f12 * rho(ji,jk) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - f12 * rho(jk,ji) * i2sv(k)%int

      else if ( l1kl ) then                        ! ( i j | k k )

        gam(ji,jk) = gam(ji,jk) - f14 * rho(jj,jk) * i2sv(k)%int
        gam(jj,jk) = gam(jj,jk) - f14 * rho(ji,jk) * i2sv(k)%int
        gam(jk,jj) = gam(jk,jj) - f14 * rho(jk,ji) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - f14 * rho(jk,jj) * i2sv(k)%int

      else if ( l1ij ) then                        ! ( i i | k l )

        gam(ji,jl) = gam(ji,jl) - f14 * rho(ji,jk) * i2sv(k)%int
        gam(ji,jk) = gam(ji,jk) - f14 * rho(ji,jl) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - f14 * rho(jl,ji) * i2sv(k)%int
        gam(jl,ji) = gam(jl,ji) - f14 * rho(jk,ji) * i2sv(k)%int

      else if ( ji == jk .and. jj == jl ) then     ! ( i j | i j )

        gam(ji,jj) = gam(ji,jj) - f14 * rho(jj,ji) * i2sv(k)%int
        gam(jj,jj) = gam(jj,jj) - f14 * rho(ji,ji) * i2sv(k)%int
        gam(ji,ji) = gam(ji,ji) - f14 * rho(jj,jj) * i2sv(k)%int
        gam(jj,ji) = gam(jj,ji) - f14 * rho(ji,jj) * i2sv(k)%int

      else if ( ji == jl .and. jj == jk ) then     ! ( i j | j i )

        gam(jj,jj) = gam(jj,jj) - f14 * rho(ji,ji) * i2sv(k)%int
        gam(ji,jj) = gam(ji,jj) - f14 * rho(jj,ji) * i2sv(k)%int
        gam(jj,ji) = gam(jj,ji) - f14 * rho(ji,jj) * i2sv(k)%int
        gam(ji,ji) = gam(ji,ji) - f14 * rho(jj,jj) * i2sv(k)%int

      else                                         ! ( i j | k l )

        gam(ji,jl) = gam(ji,jl) - f18 * rho(jj,jk) * i2sv(k)%int
        gam(jj,jl) = gam(jj,jl) - f18 * rho(ji,jk) * i2sv(k)%int
        gam(ji,jk) = gam(ji,jk) - f18 * rho(jj,jl) * i2sv(k)%int
        gam(jj,jk) = gam(jj,jk) - f18 * rho(ji,jl) * i2sv(k)%int

        gam(jk,jj) = gam(jk,jj) - f18 * rho(jl,ji) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - f18 * rho(jl,jj) * i2sv(k)%int
        gam(jl,jj) = gam(jl,jj) - f18 * rho(jk,ji) * i2sv(k)%int
        gam(jl,ji) = gam(jl,ji) - f18 * rho(jk,jj) * i2sv(k)%int

      end if

    end do

  else

    ! loop over two-electron integrals

    do k = 1, ni2s
      ji = i2sv(k)%ii
      jj = i2sv(k)%ij
      jk = i2sv(k)%ik
      jl = i2sv(k)%il

      l1ij = ji == jj
      l1kl = jk == jl

      if ( l1ij .and. l1kl .and. ji == jk ) then   ! ( i i | i i )

        gam(ji,ji) = gam(ji,ji) - rho(ji,ji) * i2sv(k)%int

      else if ( l1ij .and. l1kl ) then             ! ( i i | k k )

        gam(ji,jk) = gam(ji,jk) - rho(ji,jk) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - rho(jk,ji) * i2sv(k)%int

      else if ( l1kl ) then                        ! ( i j | k k )

        gam(ji,jk) = gam(ji,jk) - rho(jj,jk) * i2sv(k)%int
        gam(jj,jk) = gam(jj,jk) - rho(ji,jk) * i2sv(k)%int
        gam(jk,jj) = gam(jk,jj) - rho(jk,ji) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - rho(jk,jj) * i2sv(k)%int

      else if ( l1ij ) then                        ! ( i i | k l )

        gam(ji,jl) = gam(ji,jl) - rho(ji,jk) * i2sv(k)%int
        gam(ji,jk) = gam(ji,jk) - rho(ji,jl) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - rho(jl,ji) * i2sv(k)%int
        gam(jl,ji) = gam(jl,ji) - rho(jk,ji) * i2sv(k)%int

      else if ( ji == jk .and. jj == jl ) then     ! ( i j | i j )

        gam(ji,jj) = gam(ji,jj) - rho(jj,ji) * i2sv(k)%int
        gam(jj,jj) = gam(jj,jj) - rho(ji,ji) * i2sv(k)%int
        gam(ji,ji) = gam(ji,ji) - rho(jj,jj) * i2sv(k)%int
        gam(jj,ji) = gam(jj,ji) - rho(ji,jj) * i2sv(k)%int

      else if ( ji == jl .and. jj == jk ) then     ! ( i j | j i )

        gam(jj,jj) = gam(jj,jj) - rho(ji,ji) * i2sv(k)%int
        gam(ji,jj) = gam(ji,jj) - rho(jj,ji) * i2sv(k)%int
        gam(jj,ji) = gam(jj,ji) - rho(ji,jj) * i2sv(k)%int
        gam(ji,ji) = gam(ji,ji) - rho(jj,jj) * i2sv(k)%int

      else                                         ! ( i j | k l )

        gam(ji,jl) = gam(ji,jl) - rho(jj,jk) * i2sv(k)%int
        gam(jj,jl) = gam(jj,jl) - rho(ji,jk) * i2sv(k)%int
        gam(ji,jk) = gam(ji,jk) - rho(jj,jl) * i2sv(k)%int
        gam(jj,jk) = gam(jj,jk) - rho(ji,jl) * i2sv(k)%int

        gam(jk,jj) = gam(jk,jj) - rho(jl,ji) * i2sv(k)%int
        gam(jk,ji) = gam(jk,ji) - rho(jl,jj) * i2sv(k)%int
        gam(jl,jj) = gam(jl,jj) - rho(jk,ji) * i2sv(k)%int
        gam(jl,ji) = gam(jl,ji) - rho(jk,jj) * i2sv(k)%int

      end if

    end do
  end if


  return
end subroutine exchange_ctr



subroutine coul_ctr (nbct, nd, rho, gam, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | coul_ctr  --  CAJH, 04.2013                                    |
! |                                                                |
! |                                                                |
! | Driver routine in calls to coulomb_ctr.                        |
! |                                                                |
! +----------------------------------------------------------------+

  integer, intent(in) :: nbct, nd, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbct,nbct,nd), intent(in) :: rho
  complex(kind=dp), dimension(nbct,nbct,nd), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: k, nk, z


  do k = 1, ni2s, nblk
    nk = min(nblk, ni2s-k+1)
    i2sblk(1:nk) = i2sv(k:k+nk-1)

    do z = 1, nd
      call coulomb_ctr (nbct, rho(1,1,z), gam(1,1,z), i2sblk, nk, lsf)
    end do
  end do


end subroutine coul_ctr



subroutine exch_ctr (nbct, nd, rho, gam, i2sv, ni2s, lsf)

! +----------------------------------------------------------------+
! |                                                                |
! | exch_ctr  --  CAJH, 04.2013                                    |
! |                                                                |
! |                                                                |
! | Driver routine in calls to exchange_ctr.                       |
! |                                                                |
! +----------------------------------------------------------------+

  integer, intent(in) :: nbct, nd, ni2s
  logical, intent(in) :: lsf

  complex(kind=dp), dimension(nbct,nbct,nd), intent(in) :: rho
  complex(kind=dp), dimension(nbct,nbct,nd), intent(inout) :: gam

  type (i2s), dimension(ni2s), intent(in) :: i2sv


  ! other variables

  integer :: k, nk, z


  do k = 1, ni2s, nblk
    nk = min(nblk, ni2s-k+1)
    i2sblk(1:nk) = i2sv(k:k+nk-1)

    do z = 1, nd
      call exchange_ctr (nbct, rho(1,1,z), gam(1,1,z), i2sblk, nk, lsf)
    end do
  end do


end subroutine exch_ctr



subroutine aobstf (ikey, norb, nbas, xmat, mat, ldm, mat_ao, ldm_ao)

! +----------------------------------------------------------------+
! |                                                                |
! | aobstf  --  CAJH, 12.2012                                      |
! |                                                                |
! |                                                                |
! | Perform othonormal AO - standard AO basis transformations:     |
! |                                                                |
! |   ikey  =  1,  ort AO -> std AO   mat_ao = X . mat . X!        |
! |         =  2,  std AO -> ort AO   mat = X! . mat_ao . X        |
! |                                                                |
! | Here, X is the transformation matrix [ = S^(-1/2) ] among the  |
! | two basis.                                                     |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   ikey   - type of operation to perform (see above)
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   mat    - matrix in orthonormal AO basis [ in or out ]
  !   mat_ao - matrix in standard AO basis [ in or out ]
  !   ldm    - leading dimension of mat
  !   ldm_ao - leading dimension of mat_ao

  integer, intent(in) :: ikey, norb, nbas, ldm, ldm_ao

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(ldm,norb), intent(inout) :: mat
  complex(kind=dp), dimension(ldm_ao,nbas), intent(inout) :: mat_ao


  select case (ikey)

    case (1)

      ! +---------------------------+
      ! |  mat_ao  =  X . mat . X!  |
      ! +---------------------------+

      call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
           & mat, ldm, z0, scrt, nbas)
      call zgemm ('n', 'c', nbas, nbas, norb, z1, scrt, nbas, &
           & xmat, nbas, z0, mat_ao, ldm_ao)

    case (2)

      ! +---------------------------+
      ! |  mat  =  X! . mat_ao . X  |
      ! +---------------------------+

      call zgemm ('n', 'n', nbas, norb, nbas, z1, mat_ao, ldm_ao, &
           & xmat, nbas, z0, scrt, nbas)
      call zgemm ('c', 'n', norb, norb, nbas, z1, xmat, nbas, &
           & scrt, nbas, z0, mat, ldm)

    case default
      write (6, *) 'error: unrecognized ikey in aobstf.'

  end select


  return
end subroutine aobstf


end module erictr


