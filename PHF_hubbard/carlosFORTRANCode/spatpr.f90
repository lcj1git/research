

module spatpr

  use constants
  use util, only : error_alloc

  implicit none
  private

  public :: chk_irrep_cs,  spatpr_cs
  public :: chk_irrep_c2,  spatpr_c2
  public :: chk_irrep_c2h, spatpr_c2h
  public :: chk_irrep_c2v, spatpr_c2v
  public :: chk_irrep_c4v, spatpr_c4v
  public :: chk_irrep_c6v, spatpr_c6v
  public :: chk_irrep_c8v, spatpr_c8v
  public :: chk_irrep_c12v,spatpr_c12v
  public :: chk_irrep_c16v,spatpr_c16v
  public :: chk_irrep_d2,  spatpr_d2
  public :: chk_irrep_d2h, spatpr_d2h
  public :: chk_irrep_d4h, spatpr_d4h
  public :: chk_irrep_d6h, spatpr_d6h
  public :: chk_irrep_d8h, spatpr_d8h
  public :: chk_irrep_d12h,spatpr_d12h
  public :: chk_irrep_d16h,spatpr_d16h


! +----------------------------------------------------------------+
! |                                                                |
! | spatpr                                                         |
! |                                                                |
! |                                                                |
! | A collection of subroutines required to perform spatial        |
! | symmetry [point group] restoration in phfmol.                  |
! |                                                                |
! | The subroutines spatpr_* are in charge of constructing the     |
! | weight matrices corresponding to the irrep selected. Please    |
! | check such routines for the labels used for each irrep and     |
! | the ordering of the symmetry operations.                       |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | Rotation matrices should be constructed elsewhere assuming an  |
! | 'orthonormal' basis. Most rotation matrices will then be       |
! | simply permutation or permutation-like matrices.               |
! |                                                                |
! |   ( see load_rotmat for further details )                      |
! |                                                                |
! | The correct rotation matrices (in the std AO basis) are        |
! | constructed here in load_rotmat. They are then transformed to  |
! | the ort AO basis as required by the phfscf module.             |
! |                                                                |
! | Currently, rotation matrices should be provided over Cartesian |
! | basis functions.                                               |
! |                                                                |
! +----------------------------------------------------------------+


contains


subroutine chk_irrep_cs (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_cs  --  CAJH, 01.2013                                |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the Cs point group.              |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A'
      ndsptl = 1

    case (2)   ! A''
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_cs.'
      stop

  end select


  return
end subroutine chk_irrep_cs



subroutine spatpr_cs (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_cs  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ Cs group ].                               |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 2 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   sg                                       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 2

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_cs')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_cs')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A'

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1

    case (2)   ! A''

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_cs.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('cs  ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_cs



subroutine chk_irrep_c2 (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c2  --  CAJH, 08.2013                                |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C2 point group.              |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A
      ndsptl = 1

    case (2)   ! B
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c2.'
      stop

  end select


  return
end subroutine chk_irrep_c2



subroutine spatpr_c2 (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c2  --  CAJH, 08.2013                                   |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C2 group ].                               |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 2 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C2                                       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 2

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c2')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c2')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1

    case (2)   ! B

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c2.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c2  ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c2



subroutine chk_irrep_c2h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c2h  --  CAJH, 08.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C2h point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_g
      ndsptl = 1

    case (2)   ! B_g
      ndsptl = 1

    case (3)   ! A_u
      ndsptl = 1

    case (4)   ! B_u
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c2h.'
      stop

  end select


  return
end subroutine chk_irrep_c2h



subroutine spatpr_c2h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c2h  --  CAJH, 08.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C2h group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 4 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C2                                       |
! |              =   3,   i                                        |
! |              =   4,   sg                                       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 4

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c2h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c2h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  + z1

    case (2)   ! B_g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  - z1

    case (3)   ! A_u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  - z1

    case (4)   ! B_u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  + z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c2h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c2h ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c2h



subroutine chk_irrep_c2v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c2v  --  CAJH, 01.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C2v point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c2v.'
      stop

  end select


  return
end subroutine chk_irrep_c2v



subroutine spatpr_c2v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c2v  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C2v group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 4 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C2(z)                                    |
! |              =   3,   sg(xz)                                   |
! |              =   4,   sg(yz)                                   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 4

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c2v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c2v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  + z1

    case (2)   ! A_2

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  - z1

    case (3)   ! B_1

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  - z1

    case (4)   ! B_2

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  + z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c2v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c2v ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c2v



subroutine chk_irrep_c4v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c4v  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C4v point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case (5)   ! E
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c4v.'
      stop

  end select


  return
end subroutine chk_irrep_c4v



subroutine spatpr_c4v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c4v  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C4v group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 8 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C4(z)                                    |
! |              =   3,   C4(z)^2                                  |
! |              =   4,   C4(z)^3                                  |
! |              =   5,   sg_v [xz]                                |
! |              =   6,   sg_d (xy+) = sg(xz) C4(z)                |
! |              =   7,   sg_v [yz]  = sg(xz) C4(z)^2              |
! |              =   8,   sg_d (xy-) = sg(xz) C4(z)^3              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 8

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c4v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c4v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = + z1
      end do

    case (2)   ! A_2

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = - z1
      end do

    case (3)   ! B_1

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = + z1
        wgtmat(1,1,i+5)  = - z1
      end do

    case (4)   ! B_2

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = - z1
        wgtmat(1,1,i+5)  = + z1
      end do

    case (5)  ! E

      k = pi/d2

      do i = 1, 4
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+4)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+4)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c4v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c4v ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c4v



subroutine chk_irrep_c6v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c6v  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C6v point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case (5:6) ! E
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c6v.'
      stop

  end select


  return
end subroutine chk_irrep_c6v



subroutine spatpr_c6v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c6v  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C6v group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 12 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C6(z)                                    |
! |              =   ...                                           |
! |              =   6,   C6(z)^5                                  |
! |              =   7,   sg(xz)                                   |
! |              =   8,   sg(xz) C6(z)                             |
! |              =   ...                                           |
! |              =  12,   sg(xz) C6(z)^5                           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 12

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c6v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c6v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = + z1
      end do

    case (2)   ! A_2

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = - z1
      end do

    case (3)   ! B_1

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = + z1
        wgtmat(1,1,i+7)  = - z1
      end do

    case (4)   ! B_2

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = - z1
        wgtmat(1,1,i+7)  = + z1
      end do

    case (5:6) ! E

      k = real(sptirr-4,dp)*pi/3.0_dp

      do i = 1, 6
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+6)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+6)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c6v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c6v ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c6v



subroutine chk_irrep_c8v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c8v  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C8v point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case (5:7) ! E
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c8v.'
      stop

  end select


  return
end subroutine chk_irrep_c8v



subroutine spatpr_c8v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c8v  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C8v group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 16 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C8(z)                                    |
! |              =   ...                                           |
! |              =   8,   C8(z)^7                                  |
! |              =   9,   sg(xz)                                   |
! |              =  10,   sg(xz) C8(z)                             |
! |              =   ...                                           |
! |              =  16,   sg(xz) C8(z)^7                           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 16

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c8v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c8v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = + z1
      end do

    case (2)   ! A_2

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = - z1
      end do

    case (3)   ! B_1

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+9)  = - z1
      end do

    case (4)   ! B_2

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+9)  = + z1
      end do

    case (5:7) ! E

      k = real(sptirr-4,dp)*pi/4.0_dp

      do i = 1, 8
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+8)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+8)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c8v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c8v ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c8v



subroutine chk_irrep_c12v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c12v  --  CAJH, 05.2013                              |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C12v point group.            |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case (5:9)  ! E
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c12v.'
      stop

  end select


  return
end subroutine chk_irrep_c12v



subroutine spatpr_c12v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c12v  --  CAJH, 05.2013                                 |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C12v group ].                             |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 24 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C12(z)                                   |
! |              =   ...                                           |
! |              =  12,   C12(z)^11                                |
! |              =  13,   sg(xz)                                   |
! |              =  14,   sg(xz) C12(z)                            |
! |              =   ...                                           |
! |              =  24,   sg(xz) C12(z)^11                         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 24

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c12v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c12v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = + z1
      end do

    case (2)   ! A_2

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = - z1
      end do

    case (3)   ! B_1

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1
      end do

    case (4)   ! B_2

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1
      end do

    case (5:9)  ! E

      k = real(sptirr-4,dp)*pi/6.0_dp

      do i = 1, 12
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c12v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c12v', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c12v



subroutine chk_irrep_c16v (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_c16v  --  CAJH, 05.2013                              |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the C16v point group.            |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1
      ndsptl = 1

    case (2)   ! A_2
      ndsptl = 1

    case (3)   ! B_1
      ndsptl = 1

    case (4)   ! B_2
      ndsptl = 1

    case (5:11) ! E
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_c16v.'
      stop

  end select


  return
end subroutine chk_irrep_c16v



subroutine spatpr_c16v (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_c16v  --  CAJH, 05.2013                                 |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ C16v group ].                             |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 32 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C16(z)                                   |
! |              =   ...                                           |
! |              =  16,   C16(z)^15                                |
! |              =  17,   sg(xz)                                   |
! |              =  18,   sg(xz) C16(z)                            |
! |              =   ...                                           |
! |              =  32,   sg(xz) C16(z)^15                         |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 32

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_c16v')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_c16v')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = + z1
      end do

    case (2)   ! A_2

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = - z1
      end do

    case (3)   ! B_1

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+17) = - z1
      end do

    case (4)   ! B_2

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+17) = + z1
      end do

    case (5:11) ! E

      k = real(sptirr-4,dp)*pi/8.0_dp

      do i = 1, 16
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+16) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+16) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_c16v.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('c16v', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_c16v



subroutine chk_irrep_d2 (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d2  --  CAJH, 01.2013                                |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D2 point group.              |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A
      ndsptl = 1

    case (2)   ! B_1
      ndsptl = 1

    case (3)   ! B_2
      ndsptl = 1

    case (4)   ! B_3
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d2.'
      stop

  end select


  return
end subroutine chk_irrep_d2



subroutine spatpr_d2 (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d2  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D2 group ].                               |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 4 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C2(z)                                    |
! |              =   3,   C2(y)                                    |
! |              =   4,   C2(x)                                    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 4

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d2')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d2')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  + z1

    case (2)   ! B_1

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  - z1

    case (3)   ! B_2

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  - z1

    case (4)   ! B_3

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  + z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d2.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d2  ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d2



subroutine chk_irrep_d2h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d2h  --  CAJH, 01.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D2h point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_g
      ndsptl = 1

    case (2)   ! B_1g
      ndsptl = 1

    case (3)   ! B_2g
      ndsptl = 1

    case (4)   ! B_3g
      ndsptl = 1

    case (5)   ! A_u
      ndsptl = 1

    case (6)   ! B_1u
      ndsptl = 1

    case (7)   ! B_2u
      ndsptl = 1

    case (8)   ! B_3u
      ndsptl = 1

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d2h.'
      stop

  end select


  return
end subroutine chk_irrep_d2h



subroutine spatpr_d2h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d2h  --  CAJH, 01.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D2h group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 8 ]                           |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C2(z)                                    |
! |              =   3,   C2(y)                                    |
! |              =   4,   C2(x)                                    |
! |              =   5,   i                                        |
! |              =   6,   sg(xy)                                   |
! |              =   7,   sg(xz)                                   |
! |              =   8,   sg(yz)                                   |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus


  nop = 8

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d2h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d2h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  + z1
      wgtmat(1,1,5) =  + z1
      wgtmat(1,1,6) =  + z1
      wgtmat(1,1,7) =  + z1
      wgtmat(1,1,8) =  + z1

    case (2)   ! B_1g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  - z1
      wgtmat(1,1,5) =  + z1
      wgtmat(1,1,6) =  + z1
      wgtmat(1,1,7) =  - z1
      wgtmat(1,1,8) =  - z1

    case (3)   ! B_2g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  - z1
      wgtmat(1,1,5) =  + z1
      wgtmat(1,1,6) =  - z1
      wgtmat(1,1,7) =  + z1
      wgtmat(1,1,8) =  - z1

    case (4)   ! B_3g

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  + z1
      wgtmat(1,1,5) =  + z1
      wgtmat(1,1,6) =  - z1
      wgtmat(1,1,7) =  - z1
      wgtmat(1,1,8) =  + z1

    case (5)   ! A_u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  + z1
      wgtmat(1,1,5) =  - z1
      wgtmat(1,1,6) =  - z1
      wgtmat(1,1,7) =  - z1
      wgtmat(1,1,8) =  - z1

    case (6)   ! B_1u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  + z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  - z1
      wgtmat(1,1,5) =  - z1
      wgtmat(1,1,6) =  - z1
      wgtmat(1,1,7) =  + z1
      wgtmat(1,1,8) =  + z1

    case (7)   ! B_2u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  + z1
      wgtmat(1,1,4) =  - z1
      wgtmat(1,1,5) =  - z1
      wgtmat(1,1,6) =  + z1
      wgtmat(1,1,7) =  - z1
      wgtmat(1,1,8) =  + z1

    case (8)   ! B_3u

      wgtmat(1,1,1) =  + z1
      wgtmat(1,1,2) =  - z1
      wgtmat(1,1,3) =  - z1
      wgtmat(1,1,4) =  + z1
      wgtmat(1,1,5) =  - z1
      wgtmat(1,1,6) =  + z1
      wgtmat(1,1,7) =  + z1
      wgtmat(1,1,8) =  - z1

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d2h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d2h ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d2h



subroutine chk_irrep_d4h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d4h  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D4h point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1g
      ndsptl = 1

    case (2)   ! A_2g
      ndsptl = 1

    case (3)   ! B_1g
      ndsptl = 1

    case (4)   ! B_2g
      ndsptl = 1

    case (5)   ! E_g
      ndsptl = 2

    case (6)   ! A_1u
      ndsptl = 1

    case (7)   ! A_2u
      ndsptl = 1

    case (8)   ! B_1u
      ndsptl = 1

    case (9)   ! B_2u
      ndsptl = 1

    case (10)  ! E_u
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d4h.'
      stop

  end select


  return
end subroutine chk_irrep_d4h



subroutine spatpr_d4h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d4h  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D4h group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 16 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C4(z)                                    |
! |              =   3,   C4(z)^2                                  |
! |              =   4,   C4(z)^3                                  |
! |              =   5,   sg_v [xz]                                |
! |              =   6,   sg_d (xy+) = sg(xz) C4(z)                |
! |              =   7,   sg_v [yz]  = sg(xz) C4(z)^2              |
! |              =   8,   sg_d (xy-) = sg(xz) C4(z)^3              |
! |              =   9,   sg_h [xy]                                |
! |              =  10,   S4(z)      = sg(xy) C4(z)                |
! |              =  11,   S4(z)^2    = sg(xy) C4(z)^2              |
! |              =  12,   S4(z)^3    = sg(xy) C4(z)^3              |
! |              =  13,   C2(x)      = sg(xy) sg(xz)               |
! |              =  14,   C2(xy+)    = sg(xy) sg(xz) C4(z)         |
! |              =  15,   C2(y)      = sg(xy) sg(xz) C4(z)^2       |
! |              =  16,   C2(xy-)    = sg(xy) sg(xz) C4(z)^3       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 16

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d4h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d4h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1g

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = + z1
        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+12) = + z1
      end do

    case (2)   ! A_2g

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = - z1
        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+12) = - z1
      end do

    case (3)   ! B_1g

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = + z1
        wgtmat(1,1,i+5)  = - z1

        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+9)  = - z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1
      end do

    case (4)   ! B_2g

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = - z1
        wgtmat(1,1,i+5)  = + z1

        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+9)  = - z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1
      end do

    case (5)   ! E_g

      k = pi/d2

      do i = 1, 4
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+4)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+4)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+8)  = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+8)  = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+12) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+12) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
      end do

    case (6)   ! A_1u

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = - z1
        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+12) = + z1
      end do

    case (7)   ! A_2u

      do i = 1, 4
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+4)  = + z1
        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+12) = - z1
      end do

    case (8)   ! B_1u

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = - z1
        wgtmat(1,1,i+5)  = + z1

        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+9)  = + z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1
      end do

    case (9)   ! B_2u

      do i = 1, 4, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+4)  = + z1
        wgtmat(1,1,i+5)  = - z1

        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+9)  = + z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1
      end do

    case (10)  ! E_u

      k = pi/d2

      do i = 1, 4
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+4)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+4)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+8)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+8)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d4h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d4h ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d4h



subroutine chk_irrep_d6h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d6h  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D6h point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1g
      ndsptl = 1

    case (2)   ! A_2g
      ndsptl = 1

    case (3)   ! B_1g
      ndsptl = 1

    case (4)   ! B_2g
      ndsptl = 1

    case (5:6)   ! E_g
      ndsptl = 2

    case (7)   ! A_1u
      ndsptl = 1

    case (8)   ! A_2u
      ndsptl = 1

    case (9)   ! B_1u
      ndsptl = 1

    case (10)  ! B_2u
      ndsptl = 1

    case (11:12) ! E_u
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d6h.'
      stop

  end select


  return
end subroutine chk_irrep_d6h



subroutine spatpr_d6h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d6h  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D6h group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 24 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C6(z)                                    |
! |              =   ...                                           |
! |              =   6,   C6(z)^5                                  |
! |              =   7,   sg(xz)                                   |
! |              =   8,   sg(xz) C6(z)                             |
! |              =   ...                                           |
! |              =  12,   sg(xz) C6(z)^5                           |
! |              =  13,   sg(xy)                                   |
! |              =  14,   sg(xy) C6(z)                             |
! |              =   ...                                           |
! |              =  18,   sg(xy) C6(z)^5                           |
! |              =  19,   sg(xy) sg(xz)                            |
! |              =  20,   sg(xy) sg(xz) C6(z)                      |
! |              =   ...                                           |
! |              =  24,   sg(xy) sg(xz) C6(z)^5                    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 24

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d6h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d6h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1g

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = + z1
        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+18) = + z1
      end do

    case (2)   ! A_2g

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = - z1
        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+18) = - z1
      end do

    case (3)   ! B_1g

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = + z1
        wgtmat(1,1,i+7)  = - z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1

        wgtmat(1,1,i+18) = + z1
        wgtmat(1,1,i+19) = - z1
      end do

    case (4)   ! B_2g

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = - z1
        wgtmat(1,1,i+7)  = + z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1

        wgtmat(1,1,i+18) = - z1
        wgtmat(1,1,i+19) = + z1
      end do

    case (5:6)   ! E_g

      k = real(sptirr-4,dp)*pi/3.0_dp

      do i = 1, 6
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+6)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+6)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-4,2) == 1 ) then
          wgtmat(1,1:2,i+12) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+12) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+18) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+18) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+18) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+18) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        end if
      end do

    case (7)   ! A_1u

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = - z1
        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+18) = + z1
      end do

    case (8)   ! A_2u

      do i = 1, 6
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+6)  = + z1
        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+18) = - z1
      end do

    case (9)   ! B_1u

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = - z1
        wgtmat(1,1,i+7)  = + z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1

        wgtmat(1,1,i+18) = + z1
        wgtmat(1,1,i+19) = - z1
      end do

    case (10)  ! B_2u

      do i = 1, 6, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+6)  = + z1
        wgtmat(1,1,i+7)  = - z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1

        wgtmat(1,1,i+18) = - z1
        wgtmat(1,1,i+19) = + z1
      end do

    case (11:12) ! E_u

      k = real(sptirr-10,dp)*pi/3.0_dp

      do i = 1, 6
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+6)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+6)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-10,2) == 1 ) then
          wgtmat(1,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+18) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+18) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+12) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+12) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+18) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+18) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        end if
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d6h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d6h ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d6h



subroutine chk_irrep_d8h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d8h  --  CAJH, 05.2013                               |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D8h point group.             |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1g
      ndsptl = 1

    case (2)   ! A_2g
      ndsptl = 1

    case (3)   ! B_1g
      ndsptl = 1

    case (4)   ! B_2g
      ndsptl = 1

    case (5:7)   ! E_g
      ndsptl = 2

    case (8)   ! A_1u
      ndsptl = 1

    case (9)   ! A_2u
      ndsptl = 1

    case (10)  ! B_1u
      ndsptl = 1

    case (11)  ! B_2u
      ndsptl = 1

    case (12:14) ! E_u
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d8h.'
      stop

  end select


  return
end subroutine chk_irrep_d8h



subroutine spatpr_d8h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d8h  --  CAJH, 05.2013                                  |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D8h group ].                              |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 32 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C8(z)                                    |
! |              =   ...                                           |
! |              =   8,   C8(z)^7                                  |
! |              =   9,   sg(xz)                                   |
! |              =  10,   sg(xz) C8(z)                             |
! |              =   ...                                           |
! |              =  16,   sg(xz) C8(z)^7                           |
! |              =  17,   sg(xy)                                   |
! |              =  18,   sg(xy) C8(z)                             |
! |              =   ...                                           |
! |              =  24,   sg(xy) C8(z)^7                           |
! |              =  25,   sg(xy) sg(xz)                            |
! |              =  26,   sg(xy) sg(xz) C8(z)                      |
! |              =   ...                                           |
! |              =  32,   sg(xy) sg(xz) C8(z)^7                    |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 32

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d8h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d8h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1g

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+24) = + z1
      end do

    case (2)   ! A_2g

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+24) = - z1
      end do

    case (3)   ! B_1g

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+9)  = - z1

        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+17) = - z1

        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+25) = - z1
      end do

    case (4)   ! B_2g

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+9)  = + z1

        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+17) = - z1

        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+25) = + z1
      end do

    case (5:7)   ! E_g

      k = real(sptirr-4,dp)*pi/4.0_dp

      do i = 1, 8
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+8)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+8)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-4,2) == 1 ) then
          wgtmat(1,1:2,i+16) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+16) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+24) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+16) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+16) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+24) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        end if
      end do

    case (8)   ! A_1u

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+24) = + z1
      end do

    case (9)   ! A_2u

      do i = 1, 8
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+24) = - z1
      end do

    case (10)  ! B_1u

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = - z1
        wgtmat(1,1,i+9)  = + z1

        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+17) = + z1

        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+25) = - z1
      end do

    case (11)  ! B_2u

      do i = 1, 8, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+8)  = + z1
        wgtmat(1,1,i+9)  = - z1

        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+17) = + z1

        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+25) = + z1
      end do

    case (12:14) ! E_u

      k = real(sptirr-11,dp)*pi/4.0_dp

      do i = 1, 8
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+8)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+8)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-11,2) == 1 ) then
          wgtmat(1,1:2,i+16) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+16) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+24) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+16) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+16) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+24) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        end if
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d8h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d8h ', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d8h



subroutine chk_irrep_d12h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d12h  --  CAJH, 05.2013                              |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D12h point group.            |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1g
      ndsptl = 1

    case (2)   ! A_2g
      ndsptl = 1

    case (3)   ! B_1g
      ndsptl = 1

    case (4)   ! B_2g
      ndsptl = 1

    case (5:9)   ! E_g
      ndsptl = 2

    case (10)  ! A_1u
      ndsptl = 1

    case (11)  ! A_2u
      ndsptl = 1

    case (12)  ! B_1u
      ndsptl = 1

    case (13)  ! B_2u
      ndsptl = 1

    case (14:18) ! E_u
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d12h.'
      stop

  end select


  return
end subroutine chk_irrep_d12h



subroutine spatpr_d12h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d12h  --  CAJH, 05.2013                                 |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D12h group ].                             |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 48 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C12(z)                                   |
! |              =   ...                                           |
! |              =  12,   C12(z)^11                                |
! |              =  13,   sg(xz)                                   |
! |              =  14,   sg(xz) C12(z)                            |
! |              =   ...                                           |
! |              =  24,   sg(xz) C12(z)^11                         |
! |              =  25,   sg(xy)                                   |
! |              =  26,   sg(xy) C12(z)                            |
! |              =   ...                                           |
! |              =  36,   sg(xy) C12(z)^11                         |
! |              =  37,   sg(xy) sg(xz)                            |
! |              =  38,   sg(xy) sg(xz) C12(z)                     |
! |              =   ...                                           |
! |              =  48,   sg(xy) sg(xz) C12(z)^11                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 48

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d12h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d12h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1g

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+36) = + z1
      end do

    case (2)   ! A_2g

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+36) = - z1
      end do

    case (3)   ! B_1g

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1

        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+25) = - z1

        wgtmat(1,1,i+36) = + z1
        wgtmat(1,1,i+37) = - z1
      end do

    case (4)   ! B_2g

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1

        wgtmat(1,1,i+24) = + z1
        wgtmat(1,1,i+25) = - z1

        wgtmat(1,1,i+36) = - z1
        wgtmat(1,1,i+37) = + z1
      end do

    case (5:9)   ! E_g

      k = real(sptirr-4,dp)*pi/6.0_dp

      do i = 1, 12
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-4,2) == 1 ) then
          wgtmat(1,1:2,i+24) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+36) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+36) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+24) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+36) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+36) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        end if
      end do

    case (10)  ! A_1u

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+36) = + z1
      end do

    case (11)  ! A_2u

      do i = 1, 12
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+36) = - z1
      end do

    case (12)  ! B_1u

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = - z1
        wgtmat(1,1,i+13) = + z1

        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+25) = + z1

        wgtmat(1,1,i+36) = + z1
        wgtmat(1,1,i+37) = - z1
      end do

    case (13)  ! B_2u

      do i = 1, 12, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+12) = + z1
        wgtmat(1,1,i+13) = - z1

        wgtmat(1,1,i+24) = - z1
        wgtmat(1,1,i+25) = + z1

        wgtmat(1,1,i+36) = - z1
        wgtmat(1,1,i+37) = + z1
      end do

    case (14:18) ! E_u

      k = real(sptirr-13,dp)*pi/6.0_dp

      do i = 1, 12
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+12) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+12) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-13,2) == 1 ) then
          wgtmat(1,1:2,i+24) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+36) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+36) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+24) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+24) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+36) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+36) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        end if
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d12h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d12h', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d12h



subroutine chk_irrep_d16h (sptirr, ndsptl)

! +----------------------------------------------------------------+
! |                                                                |
! | chk_irrep_d16h  --  CAJH, 05.2013                              |
! |                                                                |
! |                                                                |
! | Verify whether the input variable sptirr corresponds to an     |
! | irreducible-representation of the D16h point group.            |
! |                                                                |
! | If that is the case, the dimension of the irrep is returned    |
! | in ndsptl. Otherwise, the calculation is aborted.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep

  integer, intent(in) :: sptirr
  integer, intent(out) :: ndsptl


  ndsptl = 0

  select case (sptirr)

    case (1)   ! A_1g
      ndsptl = 1

    case (2)   ! A_2g
      ndsptl = 1

    case (3)   ! B_1g
      ndsptl = 1

    case (4)   ! B_2g
      ndsptl = 1

    case (5:11)  ! E_g
      ndsptl = 2

    case (12)  ! A_1u
      ndsptl = 1

    case (13)  ! A_2u
      ndsptl = 1

    case (14)  ! B_1u
      ndsptl = 1

    case (15)  ! B_2u
      ndsptl = 1

    case (16:22) ! E_u
      ndsptl = 2

    case default
      write (6, *) 'error: Incorrect sptirr in chk_irrep_d16h.'
      stop

  end select


  return
end subroutine chk_irrep_d16h



subroutine spatpr_d16h (sptirr, ndsptl, nop, rotmat, wgtmat, norb, &
     & nbas, nbct, xmat, smat, frot)

! +----------------------------------------------------------------+
! |                                                                |
! | spatpr_d16h  --  CAJH, 05.2013                                 |
! |                                                                |
! |                                                                |
! | Fill the rotation and weight matrices necessary to perform     |
! | spatial projection [ D16h group ].                             |
! |                                                                |
! | The order of the operations used in building the weight        |
! | matrices is as follows:  [ nop = 64 ]                          |
! |                                                                |
! |   operation  =   1,   E                                        |
! |              =   2,   C16(z)                                   |
! |              =   ...                                           |
! |              =  16,   C16(z)^15                                |
! |              =  17,   sg(xz)                                   |
! |              =  18,   sg(xz) C16(z)                            |
! |              =   ...                                           |
! |              =  32,   sg(xz) C16(z)^15                         |
! |              =  33,   sg(xy)                                   |
! |              =  34,   sg(xy) C16(z)                            |
! |              =   ...                                           |
! |              =  48,   sg(xy) C16(z)^15                         |
! |              =  49,   sg(xy) sg(xz)                            |
! |              =  50,   sg(xy) sg(xz) C16(z)                     |
! |              =   ...                                           |
! |              =  64,   sg(xy) sg(xz) C16(z)^15                  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptirr - index of irrep
  !   ndsptl - dimension of irrep
  !   nop    - number of operations in group [ out ]
  !   rotmat - rotation matrices for spatial projection [ out ]
  !   wgtmat - weight matrices for spatial projection [ out ]
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: sptirr, ndsptl, norb, nbas, nbct
  integer, intent(out) :: nop

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(:,:,:), allocatable, intent(out) :: rotmat, wgtmat

  character(len=*), intent(in) :: frot


  ! other variables

  integer :: istatus, i

  real(kind=dp) :: k, pi
  parameter ( pi = 4.0_dp*atan(d1) )


  nop = 64

  allocate (rotmat(norb,norb,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'rotmat', 'spatpr_d16h')

  allocate (wgtmat(ndsptl,ndsptl,nop), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'wgtmat', 'spatpr_d16h')


  ! fill wgtmat

  select case (sptirr)

    case (1)   ! A_1g

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+32) = + z1
        wgtmat(1,1,i+48) = + z1
      end do

    case (2)   ! A_2g

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+32) = + z1
        wgtmat(1,1,i+48) = - z1
      end do

    case (3)   ! B_1g

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+17) = - z1

        wgtmat(1,1,i+32) = + z1
        wgtmat(1,1,i+33) = - z1

        wgtmat(1,1,i+48) = + z1
        wgtmat(1,1,i+49) = - z1
      end do

    case (4)   ! B_2g

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+17) = + z1

        wgtmat(1,1,i+32) = + z1
        wgtmat(1,1,i+33) = - z1

        wgtmat(1,1,i+48) = - z1
        wgtmat(1,1,i+49) = + z1
      end do

    case (5:11)  ! E_g

      k = real(sptirr-4,dp)*pi/8.0_dp

      do i = 1, 16
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+16) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+16) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-4,2) == 1 ) then
          wgtmat(1,1:2,i+32) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+32) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+48) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+48) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+32) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+32) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+48) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+48) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        end if
      end do

    case (12)  ! A_1u

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+32) = - z1
        wgtmat(1,1,i+48) = + z1
      end do

    case (13)  ! A_2u

      do i = 1, 16
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+32) = - z1
        wgtmat(1,1,i+48) = - z1
      end do

    case (14)  ! B_1u

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = - z1
        wgtmat(1,1,i+17) = + z1

        wgtmat(1,1,i+32) = - z1
        wgtmat(1,1,i+33) = + z1

        wgtmat(1,1,i+48) = + z1
        wgtmat(1,1,i+49) = - z1
      end do

    case (15)  ! B_2u

      do i = 1, 16, 2
        wgtmat(1,1,i+0)  = + z1
        wgtmat(1,1,i+1)  = - z1

        wgtmat(1,1,i+16) = + z1
        wgtmat(1,1,i+17) = - z1

        wgtmat(1,1,i+32) = - z1
        wgtmat(1,1,i+33) = + z1

        wgtmat(1,1,i+48) = - z1
        wgtmat(1,1,i+49) = + z1
      end do

    case (16:22) ! E_u

      k = real(sptirr-15,dp)*pi/8.0_dp

      do i = 1, 16
        wgtmat(1,1:2,i+0)  = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+0)  = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

        wgtmat(1,1:2,i+16) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
        wgtmat(2,1:2,i+16) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)

        if ( mod(sptirr-15,2) == 1 ) then
          wgtmat(1,1:2,i+32) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+32) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+48) = (/ +sin(k*real(i-1,dp)), +cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+48) = (/ +cos(k*real(i-1,dp)), -sin(k*real(i-1,dp)) /)
        else
          wgtmat(1,1:2,i+32) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+32) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)

          wgtmat(1,1:2,i+48) = (/ -sin(k*real(i-1,dp)), -cos(k*real(i-1,dp)) /)
          wgtmat(2,1:2,i+48) = (/ -cos(k*real(i-1,dp)), +sin(k*real(i-1,dp)) /)
        end if
      end do

    case default
      write (6, *) 'error: Incorrect sptirr in spatpr_d16h.'
      stop

  end select

  ! fill rotmat

  call load_rotmat ('d16h', nop, norb, nbas, nbct, rotmat, xmat, &
       & smat, frot)


  return
end subroutine spatpr_d16h



subroutine load_rotmat (sptgr, nop, norb, nbas, nbct, rotmat, xmat, &
     & smat, frot)

  use purcart

! +----------------------------------------------------------------+
! |                                                                |
! | load_rotmat  --  CAJH, 01.2013                                 |
! |                                                                |
! |                                                                |
! | Load rotation matrices from file frot.                         |
! |                                                                |
! | They are transformed to the ort AO basis using the xmat        |
! | transformation matrix and loaded onto the array rotmat.        |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The rotation matrices provided in file frot should be stored   |
! | as real (double precision) matrices, of dimension nbas x nbas. |
! | [ nbas is the dimension of the std AO basis. ] The rotation    |
! | matrices should be provided on input assuming an orthonormal   |
! | basis set. In most cases, these rotation matrices are          |
! | permutation or permutation-like matrices.                      |
! |                                                                |
! | NOTE: Rotation matrices should be provided over Cartesian      |
! |   functions. This routine performs the transformation to       |
! |   pure functions.                                              |
! |                                                                |
! | The true (correct) rotation matrix corresponds to              |
! |                                                                |
! |   R_ij  =  < i | R | j >,                                      |
! |         =  < i | j'>                                           |
! |                                                                |
! | where R is the symmetry operation. Since the basis is in fact  |
! | non-orthonormal, the matrix elements should be proportional to |
! | matrix elements of the overlap matrix. That overlap factor is  |
! | accounted for in this routine.                                 |
! |                                                                |
! | That is, this routine explicitly forms                         |
! |                                                                |
! |    R  =  S . Rp,                                               |
! |                                                                |
! | where Rp is the rotation matrix in the 'orthonormal' basis     |
! | read from frot, and R is the rotation matrix in the proper     |
! | std AO basis.                                                  |
! |                                                                |
! | Lastly, the rotation matrices are transformed in this routine  |
! | to the ort AO basis by computing                               |
! |                                                                |
! |   X! . R . X,                                                  |
! |                                                                |
! | where X [ = S^(-1/2) ] is the transformation matrix. These     |
! | transformed rotation matrices are stored in the array rotmat.  |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   sptgr  - point group name
  !   nop    - number of operations in group
  !   norb   - number of orbitals
  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   rotmat - rotation matrices for spatial projection [ updated ]
  !   xmat   - transformation matrix [ = S^(-1/2) ]
  !   smat   - overlap matrix ( std AO basis )
  !   frot   - file name with rotation matrices

  integer, intent(in) :: nop, norb, nbas, nbct

  complex(kind=dp), dimension(nbas,norb), intent(in) :: xmat
  complex(kind=dp), dimension(nbas,nbas), intent(in) :: smat

  complex(kind=dp), dimension(norb,norb,nop), intent(inout) :: rotmat

  character(len=*), intent(in) :: sptgr, frot


  ! other variables

  integer :: idat
  parameter ( idat = 13 )

  character(len=4) :: ngrp

  integer :: nsymop, nbct1, istatus
  integer :: j, k, j1, j2
  logical :: ftest

  real(kind=dp), dimension(:,:), allocatable :: drot

  complex(kind=dp), dimension(:,:), allocatable :: scrt, smat1, smat2
  complex(kind=dp), dimension(:,:), allocatable :: rotao, rotao2


  ! allocate memory for scratch arrays

  !   drot   - real 'orthonormal' rotation matrix
  !   rotao2 - rotation matrix in std AO basis (Cartesian)
  !   rotao  - rotation matrix in std AO basis (pure)
  !   scrt   - scratch space for basis transformations

  allocate (drot(nbct,nbct), rotao(nbas,nbas), &
          & rotao2(nbct,nbct), scrt(nbas,norb), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'scratch', 'load_rotmat')

  !   smat1  - overlap matrix over pure functions
  !   smat2  - overlap matrix over Cartesian functions

  allocate (smat1(nbas,nbas), smat2(nbct,nbct), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'smat?', 'load_rotmat')

    if ( nbct > nbas ) then
      smat1(1:nbas,1:nbas) = smat(1:nbas,1:nbas)
      call pur2cart (1, .false., nbas, nbct, smat1, smat2)
    else
      smat2(1:nbct,1:nbct) = smat(1:nbas,1:nbas)
    end if


  ! open frot file

  inquire (file = frot, exist = ftest)

  if ( .not. ftest ) then
    write (6, *) 'error: Rotation matrix file ', frot, ' does not exist.'
    stop
  end if

  open (unit = idat, file = frot, status = 'old', form = 'unformatted')


  ! read some parameters:
  !   ngrp   - name of group
  !   nsymop - number of symmetry operations
  !   nbct   - size of the basis (Cartesian)

  read (idat) ngrp

  if ( ngrp /= sptgr ) then
    write (6, *) 'error: Inconsistent sptgr in frot file.'
    stop
  end if

  read (idat) nsymop, nbct1

  if ( nsymop /= nop ) then
    write (6, *) 'error: Incorrect nop in frot file.'
    stop
  end if

  if ( nbct1 /= nbct ) then
    write (6, *) 'error: Inconsistent nbct in frot file.'
    stop
  end if


  ! loop over rotation matrices

  do k = 1, nop

    ! load 'orthonormal' matrix from file frot

    read (idat) ((drot(j1,j2), j1 = 1, nbct), j2 = 1, nbct)

    ! transform to std AO basis

    rotao2 = z0

    do j2 = 1, nbct
      do j1 = 1, nbct
        do j = 1, nbct
          rotao2(j1,j2) = rotao2(j1,j2) + smat2(j1,j) * drot(j,j2)
        end do
      end do
    end do

    ! transform to pure functions

    if ( nbct > nbas ) then
      call pur2cart (2, .false., nbas, nbct, rotao, rotao2)
    else
      rotao(1:nbas,1:nbas) = rotao2(1:nbct,1:nbct)
    end if

    ! transform to ort AO basis

    call zgemm ('n', 'n', nbas, norb, nbas, z1, rotao, nbas, &
         & xmat, nbas, z0, scrt, nbas)
    call zgemm ('c', 'n', norb, norb, nbas, z1, xmat, nbas, &
         & scrt, nbas, z0, rotmat(1,1,k), norb)

  end do


  ! close frot file

  close (unit = idat)


  ! deallocate memory for scratch arrays

  deallocate (smat1, smat2, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'smat?', 'load_rotmat')

  deallocate (drot, rotao, scrt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'scratch', 'load_rotmat')


  return
end subroutine load_rotmat


end module spatpr


