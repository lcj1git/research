

program dirty

! +----------------------------------------------------------------+
! |                                                                |
! | dirty  --  CAJH, 01.2013                                       |
! |                                                                |
! |                                                                |
! | Dirty program to read the data file produced by phfmol         |
! | [ phfrun.dat ] and print the Hamiltonian eigenvalues with      |
! | more digits than what phfmol prints.                           |
! |                                                                |
! +----------------------------------------------------------------+

  implicit none

  integer :: dp
  parameter ( dp = selected_real_kind (15, 307) )

  integer :: mxstt
  parameter ( mxstt = 20 )

  real(kind=dp), dimension(mxstt) :: stval

  integer :: nbas, nbct, norb, nup, ndn, iwfnty, imethd
  integer :: ndet, nstat
  integer :: k

  character(len=*) :: datfile
  parameter ( datfile = 'phfrun.dat' )

  integer :: idat
  parameter ( idat = 13 )

  logical :: ftest


  ! open data file

  inquire (file = datfile, exist = ftest)

  if ( .not. ftest ) then
    write (6, *) 'error: File ', datfile, ' does not exist.'
    stop
  end if

  open (unit = idat, file = datfile, status = 'old', form = 'unformatted')


  ! read integer scalars

  read (idat) nbas, nbct, norb, nup, ndn, iwfnty, imethd

  ! spin projection

  if ( mod(imethd-1,2) == 1 ) then
    read (idat)
    read (idat)
    read (idat)
  end if

  ! spatial projection

  if ( mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3 ) then
    read (idat)
    read (idat)
  end if

  ! xmat

  read (idat)

  ! determinant information

  read (idat) ndet, nstat
  read (idat)

  ! overlap, Hamiltonian matrices
  ! eigenvectors and eigenvalues

  read (idat)
  read (idat)
  read (idat)
  read (idat)
  read (idat)
  read (idat) (stval(k), k = 1, nstat)

  ! write eigenvalues to stdout

  write (6, *)
  write (6, 10) 'hmstat eigenvalues (long) :'
  write (6, *)

  10 format (X, A)

  do k = 1, nstat
    write (6, 100) k, stval(k)
  end do

  100 format (5X, I4, 3X, F20.12)

  write (6, *)


  ! close data file

  close (unit = idat)


end program dirty


