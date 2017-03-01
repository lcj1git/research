

module util

  use constants, only : dp

  implicit none
  public

! +----------------------------------------------------------------+
! |                                                                |
! | util                                                           |
! |                                                                |
! |                                                                |
! | A collection of standard utilities used by phfmol and related  |
! | programs.                                                      |
! |                                                                |
! +----------------------------------------------------------------+

contains


subroutine print_zmat (iout, nmat, idim1, idim2, str, xmat)

! +----------------------------------------------------------------+
! |                                                                |
! | print_zmat  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | Print an array, xmat, of nmat rectangular matrices (dimension  |
! | idim1 x idim2). That is, xmat = xmat(idim1,idim2,nmat).        |
! |                                                                |
! | The array zmat is a double complex array and the matrix is     |
! | printed to an output file (defined by out) in a nice format.   |
! | The header string 'str' is used for all matrices.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iout  - output file unit
  !   nmat  - number of matrices in array xmat
  !   idim1 - leading dimension of matrices in array xmat
  !   idim2 - secondary dimension (columns) of matrices in array xmat
  !   str   - header string to use for matrices in array xmat
  !   xmat  - double complex array of matrices to be printed

  integer, intent(in) :: iout, nmat, idim1, idim2

  character(len=*), intent(in) :: str

  complex(kind=dp), dimension(idim1,idim2,nmat), intent(in) :: xmat


  ! some parameters

  !   pdim   - number of columns to print at a time **
  !   thresh - printing threshold

  !  ** should be consistent with format statements

  integer :: pdim
  parameter ( pdim = 6 )

  real(kind=dp) :: thresh
  parameter ( thresh = 1.0e-14_dp )


  ! other variables

  integer :: j, y, z, k1, k2, imat
  integer :: maxy, maxk2

  real(kind=dp), dimension(pdim) :: arrprt


  ! format statements

  201 format (X, 8X, 6(12X, I4))
  202 format (X, 4X, I4, 1P, 6(2X, E14.6))

  301 format (X, A, ' (real) :')
  302 format (X, A, ' (imag) :')
  311 format (X, A, ' (real, imat = ', I3, ') :')
  312 format (X, A, ' (imag, imat = ', I3, ') :')


  ! loop over matrices in array

  do imat = 1, nmat

    ! loop over real, imaginary

    do z = 1, 2

      ! print header string of matrix
    
      if ( nmat == 1 ) then
        if ( z == 1 )  write (iout, 301)  str
        if ( z == 2 )  write (iout, 302)  str
      else
        if ( z == 1 )  write (iout, 311)  str, imat
        if ( z == 2 )  write (iout, 312)  str, imat
      end if
    
      ! external loop over columns
    
      do k2 = 1, idim2, pdim
    
        maxk2 = min (k2+pdim-1, idim2)
        maxy  = min (maxk2-k2+1, pdim)
    
        write (iout, 201)  (j, j = k2, maxk2)
    
        ! internal loop over rows
    
        do k1 = 1, idim1
    
          ! prepare printing array
    
          do y = 1, maxy
            if ( z == 1 ) then
              if ( abs (real (xmat(k1,k2-1+y,imat))) >= thresh ) then
                arrprt(y) = real (xmat(k1,k2-1+y,imat))
              else
                arrprt(y) = 0.0e0_dp
              end if

            else
              if ( abs (aimag (xmat(k1,k2-1+y,imat))) >= thresh ) then
                arrprt(y) = aimag (xmat(k1,k2-1+y,imat))
              else
                arrprt(y) = 0.0e0_dp
              end if
            end if
          end do
    
          write (iout, 202)  k1, (arrprt(j), j = 1, maxy)
    
        end do
      end do
    
      write (iout, *)
    end do
  end do


  return
end subroutine print_zmat



subroutine print_dmat (iout, nmat, idim1, idim2, str, xmat)

! +----------------------------------------------------------------+
! |                                                                |
! | print_dmat  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | Print an array, xmat, of nmat rectangular matrices (dimension  |
! | idim1 x idim2). That is, xmat = xmat(idim1,idim2,nmat).        |
! |                                                                |
! | The array zmat is a double precision array and the matrix is   |
! | printed to an output file (defined by out) in a nice format.   |
! | The header string 'str' is used for all matrices.              |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   iout  - output file unit
  !   nmat  - number of matrices in array xmat
  !   idim1 - leading dimension of matrices in array xmat
  !   idim2 - secondary dimension (columns) of matrices in array xmat
  !   str   - header string to use for matrices in array xmat
  !   xmat  - double complex array of matrices to be printed

  integer, intent(in) :: iout, nmat, idim1, idim2

  character(len=*), intent(in) :: str

  real(kind=dp), dimension(idim1,idim2,nmat), intent(in) :: xmat


  ! some parameters

  !   pdim   - number of columns to print at a time **
  !   thresh - printing threshold

  !  ** should be consistent with format statements

  integer :: pdim
  parameter ( pdim = 6 )

  real(kind=dp) :: thresh
  parameter ( thresh = 1.0e-14_dp )


  ! other variables

  integer :: j, y, k1, k2, imat
  integer :: maxy, maxk2

  real(kind=dp), dimension(pdim) :: arrprt


  ! format statements

  201 format (X, 8X, 6(12X, I4))
  202 format (X, 4X, I4, 1P, 6(2X, E14.6))

  301 format (X, A, ' :')
  311 format (X, A, ' (imat = ', I3, ') :')


  ! loop over matrices in array

  do imat = 1, nmat

    ! print header string of matrix
  
    if ( nmat == 1 ) then
      write (iout, 301)  str
    else
      write (iout, 311)  str, imat
    end if
  
    ! external loop over columns
  
    do k2 = 1, idim2, pdim
  
      maxk2 = min (k2+pdim-1, idim2)
      maxy  = min (maxk2-k2+1, pdim)
  
      write (iout, 201)  (j, j = k2, maxk2)
  
      ! internal loop over rows
  
      do k1 = 1, idim1
  
        ! prepare printing array
  
        do y = 1, maxy
          if ( abs (xmat(k1,k2-1+y,imat)) >= thresh ) then
            arrprt(y) = xmat(k1,k2-1+y,imat)
          else
            arrprt(y) = 0.0e0_dp
          end if
        end do
  
        write (iout, 202)  k1, (arrprt(j), j = 1, maxy)
  
      end do
    end do
  
    write (iout, *)
  end do


  return
end subroutine print_dmat



subroutine error_alloc (itype, array, subr)

! +----------------------------------------------------------------+
! |                                                                |
! | error_alloc  --  CAJH, 11.2012                                 |
! |                                                                |
! |                                                                |
! | Kill a job if there are problems with allocation (itype = 1)   |
! | or deallocation (itype = 2) of memory.                         |
! |                                                                |
! | The character strings 'array' and 'subr' are used to provide   |
! | the name of the array for which memory allocation failed and   |
! | the name of the subroutine where failure occurred.             |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  !   itype - type of error message to print (see above)
  !   array - name of array for which error was produced
  !   subr  - name of subroutine in which error was produced

  integer, intent(in) :: itype

  character(len=*), intent(in) :: array, subr


  select case (itype)

    case (1)
      write (6, *) 'error: Could not allocate memory for array ', &
                 & array, ' in subroutine ', subr, '.'

    case (2)
      write (6, *) 'error: Could not de-allocate memory for array ', &
                 & array, ' in subroutine ', subr, '.'

    case default
      write (6, *) 'error: itype not recognized in error_alloc.'

  end select


  stop
end subroutine error_alloc


end module util


