

module purcart

  use constants
  use util, only : error_alloc

  implicit none
  private
  save

  public :: setup_purcart, shutdown_purcart
  public :: pur2cart

! +----------------------------------------------------------------+
! |                                                                |
! | purcart                                                        |
! |                                                                |
! |                                                                |
! | A module used to oversee the transformation of matrices from   |
! | pure to Cartesian functions and viceversa.                     |
! |                                                                |
! | Once the transformation matrices have been loaded once in      |
! | the program (by a call to setup_purcart), the subroutine       |
! | pur2cart can be accessed anywhere else in the program to       |
! | perform the transformation.                                    |
! |                                                                |
! | Arrays can be deallocated (and memory cleared) by a call to    |
! | shutdown_purcart.                                              |
! |                                                                |
! +----------------------------------------------------------------+


  ! store pure to Cartesian transformation matrices
  ! scratch arrays

  complex(kind=dp), dimension(:), allocatable :: x2pure, x2cart
  complex(kind=dp), dimension(:), allocatable :: scrt

  !$omp  threadprivate(x2pure, x2cart, scrt)


contains


subroutine setup_purcart (nbas, nbct, x2pur1, x2crt1)

! +----------------------------------------------------------------+
! |                                                                |
! | setup_purcart  --  CAJH, 02.2013                               |
! |                                                                |
! |                                                                |
! | Load transformation matrices (pure to Cartesian and viceversa) |
! | into the module. Allocate memory for a scratch array.          |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   nbas   - number of basis functions
  !   nbct   - number of Cartesian basis functions
  !   x2pur1 - transformation matrix Cartesian -> pure
  !   x2crt1 - transformation matrix pure -> Cartesian

  integer, intent(in) :: nbas, nbct

  complex(kind=dp), dimension(nbct,nbas), intent(in) :: x2pur1
  complex(kind=dp), dimension(nbas,nbct), intent(in) :: x2crt1


  ! other variables

  integer :: istatus


  allocate (x2pure(nbct*nbas), x2cart(nbct*nbas), &
          & scrt(nbct*nbas), stat=istatus)
  if ( istatus /= 0 ) call error_alloc (1, 'all', 'setup_purcart')

  x2pure(1:nbct*nbas) = &
       reshape (x2pur1(1:nbct,1:nbas), (/ nbct*nbas /))
  x2cart(1:nbct*nbas) = &
       reshape (x2crt1(1:nbas,1:nbct), (/ nbct*nbas /))


  return
end subroutine setup_purcart



subroutine shutdown_purcart

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_purcart  --  CAJH, 02.2013                            |
! |                                                                |
! |                                                                |
! | Deallocate memory for arrays used in the purcart module.       |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: istatus

  deallocate (x2pure, x2cart, scrt, stat=istatus)
  if ( istatus /= 0 ) call error_alloc (2, 'all', 'shutdown_purcart')

  return
end subroutine shutdown_purcart



subroutine pur2cart (ikey, lden, nbas, nbct, mtpur, mtct)

! +----------------------------------------------------------------+
! |                                                                |
! | pur2cart  --  CAJH, 02.2013                                    |
! |                                                                |
! |                                                                |
! | Perform pure <-> Cartesian AO basis transformations:           |
! |                                                                |
! |   ikey  =  1,  pure -> Cartesian                               |
! |         =  2,  Cartesian -> pure                               |
! |                                                                |
! | The logical variable lden is used to determine the type of     |
! | transformation to perform. lden = .true. implies that the      |
! | matrix transforms as a density matrix:                         |
! |                                                                |
! |   lden = F,  ikey  =  1,  mtct  =  x2cart! . mtpur . x2cart    !
! |                    =  2,  mtpur =  x2pure! . mtct  . x2pure    |
! |                                                                |
! |   lden = T,  ikey  =  1,  mtct  =  x2pure . mtpur . x2pure!    |
! |                           mtpur =  x2cart . mtct  . x2cart!    |
! |                                                                |
! | Here, x2pure and x2cart are transformation matrices previously |
! | loaded into the module.                                        |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   ikey  - type of transformation to perform (see above)
  !   lden  - whether input matrix transforms as a density matrux
  !   nbas  - number of basis functions
  !   nbct  - number of Cartesian basis functions
  !   mtpur - input / output matrix in pure basis
  !   mtct  - input / output matrix in Cartesian basis

  integer, intent(in) :: ikey, nbas, nbct
  logical, intent(in) :: lden

  complex(kind=dp), dimension(nbas,nbas), intent(inout) :: mtpur
  complex(kind=dp), dimension(nbct,nbct), intent(inout) :: mtct


  if ( .not. lden ) then

    select case (ikey)

      case (1)

        ! +-------------------------------------+
        ! |  mtct  =  x2cart! . mtpur . x2cart  |
        ! +-------------------------------------+

        call zgemm ('n', 'n', nbas, nbct, nbas, z1, mtpur, nbas, &
             & x2cart, nbas, z0, scrt, nbas)
        call zgemm ('c', 'n', nbct, nbct, nbas, z1, x2cart, nbas, &
             & scrt, nbas, z0, mtct, nbct)

      case (2)

        ! +-------------------------------------+
        ! |  mtpur  =  x2pure! . mtct . x2pure  |
        ! +-------------------------------------+

        call zgemm ('n', 'n', nbct, nbas, nbct, z1, mtct, nbct, &
             & x2pure, nbct, z0, scrt, nbct)
        call zgemm ('c', 'n', nbas, nbas, nbct, z1, x2pure, nbct, &
             & scrt, nbct, z0, mtpur, nbas)

      case default
        write (6, *) 'error: unrecognized ikey in pur2cart.'

    end select

  else

    select case (ikey)

      case (1)

        ! +-------------------------------------+
        ! |  mtct  =  x2pure . mtpur . x2pure!  |
        ! +-------------------------------------+

        call zgemm ('n', 'c', nbas, nbct, nbas, z1, mtpur, nbas, &
             & x2pure, nbct, z0, scrt, nbas)
        call zgemm ('n', 'n', nbct, nbct, nbas, z1, x2pure, nbct, &
             & scrt, nbas, z0, mtct, nbct)

      case (2)

        ! +-------------------------------------+
        ! |  mtpur  =  x2cart . mtct . x2cart!  |
        ! +-------------------------------------+

        call zgemm ('n', 'c', nbct, nbas, nbct, z1, mtct, nbct, &
             & x2cart, nbas, z0, scrt, nbct)
        call zgemm ('n', 'n', nbas, nbas, nbct, z1, x2cart, nbas, &
             & scrt, nbct, z0, mtpur, nbas)

      case default
        write (6, *) 'error: unrecognized ikey in pur2cart.'

    end select

  end if


  return
end subroutine pur2cart


end module purcart


