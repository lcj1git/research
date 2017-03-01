

program linmol

! +----------------------------------------------------------------+
! |                                                                |
! | linmol  --  CAJH, 05.2013                                      |
! |                                                                |
! |                                                                |
! | Program to construct the rotation matrices needed to perform   |
! | spatial projection in phfmol for LINEAR MOLECULES.             |
! |                                                                |
! |   ( This program replaces rotmal_linmol. )                     |
! |                                                                |
! | This program assumes that the molecule is oriented along the   |
! | z-axis. That is, it assumes that the coordinates of atoms are  |
! |                                                                |
! |    A[1]    0.0   0.0    z(1)                                   |
! |    A[2]    0.0   0.0    z(2)                                   |
! |    ...                                                         |
! |    A[N-1]  0.0   0.0    z(n-1)                                 |
! |    A[N]    0.0   0.0    z(n)                                   |
! |                                                                |
! | The program supports two types of linear systems:              |
! |                                                                |
! |   linear molecules with sig[xy] symmetry  [ PG = Dinfh ]       |
! |   linear molecules w/o  sig[xy] symmetry  [ PG = Cinfv ]       |
! |                                                                |
! | The code lets the user provides the following information:     |
! |                                                                |
! |   - whether the molecule has sig[xy] symmetry                  |
! |   - number of atoms                                            |
! |   - number of shells per atom                                  |
! |   - angular momentum of each shell                             |
! |   - number of grid points to use **                            |
! |                                                                |
! | ** This determines the number of angles to use in the Cinf     |
! |    rotations allowed for linear molecules. The number of grid  |
! |    points determines the actual group used:                    |
! |                                                                |
! |      ngrd  =  1,  =>   group :  Cs,  D2                        |
! |            =  2,  =>            C2v, D2h                       |
! |            =  4,  =>            C4v, D4h                       |
! |            ...                                                 |
! |                                                                |
! | The code assumes that Cartesian basis functions are used.      |
! | ( Below we describe the ordering of the basis functions for    |
! | each type of shell. ) The code can support shells with         |
! | arbitrary angular momentum.                                    |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | The code outputs the rotation matrices in the format required  |
! | by phfmol in the 'data' (unformatted) file                     |
! |                                                                |
! |   linmol.rot                                                   |
! |                                                                |
! | Additionally, if debug is set to true, the code prints the     |
! | rotation matrices in ASCII format to the 'output' file         |
! |                                                                |
! |   linmol.out                                                   |
! |                                                                |
! +----------------------------------------------------------------+

! +----------------------------------------------------------------+
! |                                                                |
! | ordering of states in each shell                               |
! | ================================                               |
! |                                                                |
! | This is the order used in Gaussian. For P, D, F shells, a      |
! | particular ordering is used if lpdf is set to true.            |
! |                                                                |
! | L = 0  ==  S functions                                         |
! |                                                                |
! |  0                                                             |
! |                                                                |
! | L = 1  ==  P functions                                         |
! |                                                                |
! |  Z    if  lpdf = .true. , then   X                             |
! |  Y                               Y                             |
! |  X                               Z                             |
! |                                                                |
! | L = 2  ==  D functions                                         |
! |                                                                |
! |  ZZ   if  lpdf = .true. , then   XX                            |
! |  YZ                              YY                            |
! |  YY                              ZZ                            |
! |  XZ                              XY                            |
! |  XY                              XZ                            |
! |  XX                              YZ                            |
! |                                                                |
! | L = 3  ==  F functions                                         |
! |                                                                |
! |  ZZZ  if  lpdf = .true. , then   XXX                           |
! |  YZZ                             YYY                           |
! |  YYZ                             ZZZ                           |
! |  YYY                             XYY                           |
! |  XZZ                             XXY                           |
! |  XYZ                             XXZ                           |
! |  XYY                             XZZ                           |
! |  XXZ                             YZZ                           |
! |  XXY                             YYZ                           |
! |  XXX                             XYZ                           |
! |                                                                |
! +----------------------------------------------------------------+

  implicit none

  ! constants

  integer :: dp
  parameter ( dp = selected_real_kind (15, 307) )

  real(kind=dp) :: d0, d1
  parameter ( d0 = 0.0e0_dp, d1 = 1.0e0_dp )


  ! some parameters

    ! maxnat - maximum number of atoms
    ! maxats - maximum shells per atom
    ! maxang - maximum angular momentum supported
    ! maxgrd - maximum number of grid points allowed

  integer :: maxnat
  parameter ( maxnat = 10 )

  integer :: maxats
  parameter ( maxats = 40 )

  integer :: maxang, maxdim
  parameter ( maxang = 6 )
  parameter ( maxdim = (maxang+1)*(maxang+2)/2 )

  integer :: maxgrd
  parameter ( maxgrd = 32 )

    ! lpdf  - whether to use special ordering for P, D, F shells
    ! debug - whether to turn debugging on

  logical :: lpdf, debug
  parameter ( lpdf = .true., debug = .false. )


  ! output, data files

  character(len=*) :: outfile, datfile
  parameter ( outfile = 'linmol.out', datfile = 'linmol.rot' )

  integer :: iout, idat
  parameter ( iout = 12, idat = 13 )

  logical :: ftest


  ! other variables

  integer :: nat, nop, ndim, ngrd
  integer :: maxj, maxshl, maxshd
  integer :: j, k, ind, ind1, ind2, ind3, ind4
  integer :: iop, igrd, jshl
  logical :: lxy, lxz, lsgxy

  real(kind=dp) :: angle, pi
  parameter ( pi = 4.0_dp*atan(d1) )

  character(len=1) :: cgrd1
  character(len=2) :: cgrd2
  character(len=4) :: namgrp

  integer, dimension(maxnat) :: nshlat, indat
  integer, dimension(maxats,maxnat) :: shl

  integer, dimension(3,maxdim,0:maxang) :: shlxyz
  real(kind=dp), dimension(maxdim,0:maxang) :: shlnorm

  real(kind=dp), dimension(3,3) :: rotxyz
  real(kind=dp), dimension(:,:,:), allocatable :: rotmat, shlrot


  ! number of states per shell

  integer :: l, ishl

  ishl(l) = (l+1)*(l+2)/2


  ! build x-y-z indices of basis states
  ! build normalization factors for basis states

    shlxyz  = 0
    shlnorm = d0

  do j = 0, maxang
    call ordering (j, shlxyz(1,1,j))
    call normalize (j, shlxyz(1,1,j), shlnorm(1,j))
  end do


  ! plane of symmetry
  ! =================

  write (6, *)
  write (6, *) 'Is the molecule symmetric under sg(xy)?'
  read  (5, *) lsgxy


  ! number of atoms
  ! ===============

  write (6, *)
  write (6, *) 'Provide the number of atoms:'
  read  (5, *) nat

  if ( nat > maxnat ) then
    write (6, *)
    write (6, *) 'error: nat > maxnat.'
    stop
  else if ( nat < 1 ) then
    write (6, *)
    write (6, *) 'error: nat < 1.'
    stop
  end if

    if ( lsgxy ) then
      maxj = (nat+1)/2
    else
      maxj = nat
    end if


  ! number of shells per atom
  ! =========================

  do j = 1, maxj
    write (6, *)
    write (6, 11) j
    read  (5, *) nshlat(j)

    if ( lsgxy ) then
      nshlat(nat+1-j) = nshlat(j)
    end if
  end do

  11 format (X, 'Provide the number of shells for atom ', I2, ':')

  ! check that nshlat is within bounds

  do j = 1, maxj
    if ( nshlat(j) < 0 ) then
      write (6, *)
      write (6, 101) j
      stop
    else if ( nshlat(j) > maxats ) then
      write (6, *)
      write (6, 102) j
      stop
    end if
  end do

  101 format (X, 'error: nshlat (', I2, ') < 0.')
  102 format (X, 'error: nshlat (', I2, ') > maxats.')


  ! angular momentum of each shell
  ! ==============================

  do j = 1, maxj
    write (6, *)
    write (6, 12) j
    read  (5, *) (shl(k,j), k = 1, nshlat(j))

    if ( lsgxy ) then
      shl(1:nshlat(j),nat+1-j) = shl(1:nshlat(j),j)
    end if
  end do

  12 format (X, 'Provide the ang mom of each shell for atom ', I2, ':')

  ! check that angular momentum is within bounds

    maxshl = 0

  do j = 1, maxj
    do k = 1, nshlat(j)
      if ( shl(k,j) < 0 ) then
        write (6, *)
        write (6, 111) k, j
        stop
      else if ( shl(k,j) > maxang ) then
        write (6, *)
        write (6, 112) k, j
        stop
      end if

      maxshl = max (maxshl, shl(k,j))
    end do
  end do

    maxshd = ishl(maxshl)

  111 format (X, 'error: shl (', I3, ',', I2, ') < 0.')
  112 format (X, 'error: shl (', I3, ',', I2, ') > maxang.')


  ! determine number of states
  ! fill indat array (first state of each atom)

    ndim = 0

  do j = 1, nat
    indat(j) = ndim + 1
    do k = 1, nshlat(j)
      jshl = shl(k,j)
      ndim = ndim + ishl(jshl)
    end do
  end do


  ! number of grid points
  ! =====================

  write (6, *)
  write (6, *) 'Provide the number of grid points to use:'
  read  (5, *) ngrd

  if ( ngrd > maxgrd ) then
    write (6, *)
    write (6, *) 'error: ngrd > maxgrd.'
    stop
  else if ( ngrd < 1 ) then
    write (6, *)
    write (6, *) 'error: ngrd < 1.'
    stop
  end if


  ! determine number of operations

  if ( lsgxy ) then
    nop = 4*ngrd
  else
    nop = 2*ngrd
  end if


  ! allocate memory for rotation matrix

  allocate (rotmat(ndim,ndim,nop))
  rotmat = d0

  ! allocate memory for shell rotation matrices

  allocate (shlrot(maxshd,maxshd,0:maxshl))
  shlrot = d0


  ! loop over rotation operations
  ! =============================

  do iop = 1, nop

    igrd = mod(iop-1,ngrd) + 1
    lxz  = (iop-1)/ngrd == 1 .or. (iop-1)/ngrd == 3
    lxy  = (iop-1)/ngrd == 2 .or. (iop-1)/ngrd == 3


    ! construct rotation matrix

    angle = 2.0_dp * pi * real(igrd-1,dp) / real(ngrd,dp)

    rotxyz(1,:) = (/ +cos(angle), -sin(angle), d0 /)
    rotxyz(2,:) = (/ +sin(angle), +cos(angle), d0 /)
    rotxyz(3,:) = (/ d0, d0, d1 /)

      if ( lxz ) then
        rotxyz(:,2) = -rotxyz(:,2)
      end if

      if ( lxy ) then
        rotxyz(:,3) = -rotxyz(:,3)
      end if


    ! construct shell rotation matrices

    do j = 0, maxshl
      call shell_rotation (j, shlxyz(1,1,j), rotxyz, shlrot(1,1,j), maxshd)
      call fix_norm (j, shlnorm(1,j), shlrot(1,1,j), maxshd)
    end do


    ! loop over atoms and shells to construct operation

    do j = 1, nat
      ind = 0

      do k = 1, nshlat(j)
        jshl = shl(k,j)

        ! determine indices

        ind1 = indat(j) + ind
        ind2 = indat(j) + ind + ishl(jshl)-1

        if ( .not. lxy ) then
          ind3 = ind1
          ind4 = ind2
        else
          ind3 = indat(nat+1-j) + ind
          ind4 = indat(nat+1-j) + ind + ishl(jshl)-1
        end if

        rotmat(ind1:ind2,ind3:ind4,iop) = &
             shlrot(1:ishl(jshl),1:ishl(jshl),jshl)

        ind = ind + ishl(jshl)
      end do
     end do

  end do


  ! save rotation matrices
  ! ======================

  ! determine the name of the group

  if ( ngrd == 1 ) then
    if ( lsgxy ) then
      namgrp = 'd2  '
    else
      namgrp = 'cs  '
    end if

  else if ( ngrd < 10 ) then
      write (cgrd1, '(I1)') ngrd
    if ( lsgxy ) then
      namgrp = 'd' // cgrd1 // 'h '
    else
      namgrp = 'c' // cgrd1 // 'v '
    end if

  else
      write (cgrd2, '(I2)') ngrd
    if ( lsgxy ) then
      namgrp = 'd' // cgrd2 // 'h'
    else
      namgrp = 'c' // cgrd2 // 'v'
    end if
  end if


  ! save rotation matrices to dat file

  inquire (file = datfile, exist = ftest)

  if ( ftest ) then
    open (unit = idat, file = datfile, &
        & status = 'old', form = 'unformatted')
    close (unit = idat, status = 'delete')
  end if

  open (unit = idat, file = datfile, &
      & status = 'new', form = 'unformatted')

  write (idat) namgrp
  write (idat) nop, ndim

  if ( namgrp == 'd2h ' ) then

    ! order for D2h: E, C2(z), C2(y), C2(x), i, sg(xy), sg(xz), sg(yz)

    write (idat) ((rotmat(j,k,1), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,2), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,8), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,7), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,6), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,5), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,3), j = 1, ndim), k = 1, ndim)
    write (idat) ((rotmat(j,k,4), j = 1, ndim), k = 1, ndim)

  else
    do iop = 1, nop
      write (idat) ((rotmat(j,k,iop), j = 1, ndim), k = 1, ndim)
    end do
  end if

  close (unit = idat)


  ! print matrices to output file

  if ( debug ) then
    inquire (file = outfile, exist = ftest)

    if ( ftest ) then
      open (unit = iout, file = outfile, &
          & status = 'old', form = 'formatted')
      close (unit = iout, status = 'delete')
    end if

    open (unit = iout, file = outfile, &
        & status = 'new', form = 'formatted')

    write (iout, *)

    if ( namgrp == 'd2h ' ) then
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   1)', rotmat(1,1,1))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   2)', rotmat(1,1,2))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   3)', rotmat(1,1,8))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   4)', rotmat(1,1,7))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   5)', rotmat(1,1,6))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   6)', rotmat(1,1,5))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   7)', rotmat(1,1,3))
      call print_dmat (iout, 1, ndim, ndim, 'rotmat (imat =   8)', rotmat(1,1,4))

    else
      call print_dmat (iout, nop, ndim, ndim, 'rotmat', rotmat)
    end if

    close (unit = iout)
  end if


  ! deallocate memory

  deallocate (rotmat, shlrot)


  stop

contains


subroutine shell_rotation (l, shlxyz, rotxyz, shlrot, sdim)

! +----------------------------------------------------------------+
! |                                                                |
! | shell_rotation  --  CAJH, 05.2013                              |
! |                                                                |
! |                                                                |
! | Given an x-y-z rotation matrix and the x-y-z indices of the    |
! | basis states in a shell, construct the rotation matrix among   |
! | the basis states in the shell.                                 |
! |                                                                |
! +----------------------------------------------------------------+


  ! input / output variables

  !   l      - angular momentum of shell
  !   shlxyz - x-y-z indices of basis states in shell
  !   rotxyz - x-y-z rotation matrix
  !   shlrot - rotation matrix among basis states in shell
  !   sdim   - dimension of shlrot array

  integer, intent(in) :: l, sdim
  integer, dimension(3,(l+1)*(l+2)/2), intent(in) :: shlxyz
  real(kind=dp), dimension(3,3), intent(in) :: rotxyz
  real(kind=dp), dimension(sdim,sdim), intent(inout) :: shlrot

  ! other variables

  integer :: j1, j2, jdim
  integer :: nx, kx, lx, ny, ky, ly, nz, kz, lz
  integer :: cfix, cfiy, cfiz
  real(kind=dp) :: cfdx, cfdy, cfdz, cft

  integer, dimension(3) :: shlx, shly, shlz, shlt


  if ( l == 0 ) then
    shlrot(1,1) = d1
    return
  end if

  jdim = (l+1)*(l+2)/2

    shlrot = d0

  do j1 = 1, jdim
    nx = shlxyz(1,j1)
    ny = shlxyz(2,j1)
    nz = shlxyz(3,j1)

    do kx = 0, nx
    do lx = 0, kx

      cfix = trinomial (nx, nx-kx, kx-lx, lx)
      cfdx = rotxyz(1,1)**(nx-kx) * &
           & rotxyz(1,2)**(kx-lx) * &
           & rotxyz(1,3)**(lx)

      shlx = (/ nx-kx, kx-lx, lx /)

      do ky = 0, ny
      do ly = 0, ky

        cfiy = trinomial (ny, ny-ky, ky-ly, ly)
        cfdy = rotxyz(2,1)**(ny-ky) * &
             & rotxyz(2,2)**(ky-ly) * &
             & rotxyz(2,3)**(ly)

        shly = (/ ny-ky, ky-ly, ly /)

        do kz = 0, nz
        do lz = 0, kz

          cfiz = trinomial (nz, nz-kz, kz-lz, lz)
          cfdz = rotxyz(3,1)**(nz-kz) * &
               & rotxyz(3,2)**(kz-lz) * &
               & rotxyz(3,3)**(lz)

          shlz = (/ nz-kz, kz-lz, lz /)

          cft = real(cfix*cfiy*cfiz,dp) * cfdx*cfdy*cfdz
          shlt(:) = shlx(:) + shly(:) + shlz(:)

          do j2 = 1, jdim
            if ( shlt(1) == shlxyz(1,j2) .and. &
                 shlt(2) == shlxyz(2,j2) .and. &
                 shlt(3) == shlxyz(3,j2) ) then

              shlrot(j2,j1) = shlrot(j2,j1) + cft
            end if
          end do

        end do
        end do
      end do
      end do
    end do
    end do

  end do


  return
end subroutine shell_rotation



subroutine fix_norm (l, shlnorm, shlrot, sdim)

! +----------------------------------------------------------------+
! |                                                                |
! | fix_norm  --  CAJH, 05.2013                                    |
! |                                                                |
! |                                                                |
! | Fix rotation matrix among basis states in a given shell by     |
! | accounting for the normalization factors of primitives from    |
! | the angular momentum exponents.                                |
! |                                                                |
! +----------------------------------------------------------------+


  ! input / output variables

  !   l       - angular momentum of shell
  !   shlnorm - normalization factor for basis states in shell
  !   shlrot  - rotation matrix among basis states in shell
  !   sdim    - dimension of shlrot array

  integer, intent(in) :: l, sdim
  real(kind=dp), dimension((l+1)*(l+2)/2), intent(in) :: shlnorm
  real(kind=dp), dimension(sdim,sdim), intent(inout) :: shlrot

  ! other variables

  integer :: j, k, jdim


  if ( l == 0 ) then
    return
  end if

  jdim = (l+1)*(l+2)/2

  do j = 1, jdim
    do k = 1, jdim
      shlrot(j,k) = shlrot(j,k) * shlnorm(k) / shlnorm(j)
    end do
  end do


  return
end subroutine fix_norm



subroutine ordering (l, shlxyz)

! +----------------------------------------------------------------+
! |                                                                |
! | ordering  --  CAJH, 05.2013                                    |
! |                                                                |
! |                                                                |
! | Given the angular momentum of a shell, return the x-y-z        |
! | indices of basis states in the shell according to the order    |
! | described in the top of the program.                           |
! |                                                                |
! +----------------------------------------------------------------+

  ! input / output variables

  !   l      - angular momentum of shell
  !   shlxyz - x-y-z indices of basis states in shell

  integer, intent(in) :: l
  integer, dimension(3,(l+1)*(l+2)/2), intent(inout) :: shlxyz

  ! other variables

  integer :: i, jx, jy


  if ( l == 0 ) then
    shlxyz(1:3,1) = (/ 0, 0, 0 /)
    return
  end if

  if ( lpdf .and. l <= 3 ) then

    select case (l)
      case (1)
        shlxyz(1:3,1)  = (/ 1, 0, 0 /)
        shlxyz(1:3,2)  = (/ 0, 1, 0 /)
        shlxyz(1:3,3)  = (/ 0, 0, 1 /)

      case (2)
        shlxyz(1:3,1)  = (/ 2, 0, 0 /)
        shlxyz(1:3,2)  = (/ 0, 2, 0 /)
        shlxyz(1:3,3)  = (/ 0, 0, 2 /)
        shlxyz(1:3,4)  = (/ 1, 1, 0 /)
        shlxyz(1:3,5)  = (/ 1, 0, 1 /)
        shlxyz(1:3,6)  = (/ 0, 1, 1 /)

      case (3)
        shlxyz(1:3,1)  = (/ 3, 0, 0 /)
        shlxyz(1:3,2)  = (/ 0, 3, 0 /)
        shlxyz(1:3,3)  = (/ 0, 0, 3 /)
        shlxyz(1:3,4)  = (/ 1, 2, 0 /)
        shlxyz(1:3,5)  = (/ 2, 1, 0 /)
        shlxyz(1:3,6)  = (/ 2, 0, 1 /)
        shlxyz(1:3,7)  = (/ 1, 0, 2 /)
        shlxyz(1:3,8)  = (/ 0, 1, 2 /)
        shlxyz(1:3,9)  = (/ 0, 2, 1 /)
        shlxyz(1:3,10) = (/ 1, 1, 1 /)
    end select

  else
      i = 1
    do jx = 0, l
      do jy = 0, l-jx
        shlxyz(1:3,i) = (/ jx, jy, l-jx-jy /)
        i = i+1
      end do
    end do
  end if


  return
end subroutine ordering



subroutine normalize (l, shlxyz, shlnorm)

! +----------------------------------------------------------------+
! |                                                                |
! | normalize  --  CAJH, 05.2013                                   |
! |                                                                |
! |                                                                |
! | Compute normalization factors (due to angular momentum         |
! | exponents) for basis states in a given shell. Here, the        |
! | normalization factors are given by                             |
! |                                                                |
! |   1 / sqrt [ (2*nx-1)!! (2*ny-1)!! (2*nz-1)!! ],               |
! |                                                                |
! | where !! denotes a double factorial and nx,ny,nz are the       |
! | angular momentum exponents of a basis state.                   |
! |                                                                |
! +----------------------------------------------------------------+


  ! input / output variables

  !   l       - angular momentum of shell
  !   shlxyz  - x-y-z indices of basis states in shell
  !   shlnorm - normalization factor for basis states in shell

  integer, intent(in) :: l
  integer, dimension(3,(l+1)*(l+2)/2), intent(in) :: shlxyz
  real(kind=dp), dimension((l+1)*(l+2)/2), intent(inout) :: shlnorm

  ! other variables

  integer :: j, jdim
  integer :: jx, jy, jz


  if ( l == 0 ) then
    shlnorm(1) = d1
    return
  end if

  jdim = (l+1)*(l+2)/2

  do j = 1, jdim
    jx = shlxyz(1,j)
    jy = shlxyz(2,j)
    jz = shlxyz(3,j)

    shlnorm(j) = d1 / sqrt (real (double_fact (2*jx-1) * &
                                & double_fact (2*jy-1) * &
                                & double_fact (2*jz-1), dp) )
  end do


  return
end subroutine normalize



function trinomial (n, k1, k2, k3)

! +----------------------------------------------------------------+
! |                                                                |
! | trinomial  --  CAJH, 05.2013                                   |
! |                                                                |
! |                                                                |
! | Compute the trinomial :                                        |
! |                                                                |
! |               n!                                               |
! |   z  =  --------------- ,                                      |
! |          k1!  k2!  k3!                                         |
! |                                                                |
! | with (k1,k2,k3) >= 0 and n = k1 + k2 + k3.                     |
! |                                                                |
! +----------------------------------------------------------------+


  ! input / output variables

  integer, intent(in) :: n, k1, k2, k3
  integer :: trinomial

  ! other variables

  integer :: j1, j2, j3
  integer :: jt, z, i


  ! error checking

  if ( k1 < 0 .or. k2 < 0 .or. k3 < 0 ) then
    write (6, *)
    write (6, *) 'error: Incorrect input arguments in trinomial.'
    stop
  end if

  if ( n /= k1 + k2 + k3 ) then
    write (6, *)
    write (6, *) 'error: Incorrect input arguments in trinomial.'
    stop
  end if

  ! order k1, k2, k3

  j1 = k1;  j2 = k2;  j3 = k3

  if ( j2 < j3 ) then
    jt = j2
    j2 = j3
    j3 = jt
  end if

  if ( j1 < j2 ) then
    jt = j1
    j1 = j2

    if ( jt < j3 ) then
      j2 = j3
      j3 = jt
    else
      j2 = jt
    end if
  end if

    z = 1

  ! compute z = n! / j1!

  do i = n, j1+1, -1
    z = z*i
  end do

  ! compute z = z / j2!

  do i = j2, 2, -1
    z = z/i
  end do

  ! compute z = z / j3!

  do i = j3, 2, -1
    z = z/i
  end do

  trinomial = z


  return
end function trinomial



function double_fact (n)

! +----------------------------------------------------------------+
! |                                                                |
! | double_fact  --  CAJH, 05.2013                                 |
! |                                                                |
! |                                                                |
! | Evaluate the double factorial of an odd number n, n >= -1.     |
! |                                                                |
! +----------------------------------------------------------------+


  ! input / output variables

  integer, intent(in) :: n
  integer :: double_fact

  ! other variables

  integer :: i, z


  ! error checking

  if ( mod(n,2) == 0 .or. n < -1 ) then
    write (6, *)
    write (6, *) 'error: Incorrect input argument in double_fact.'
    stop
  end if

    z = 1

  do i = n, 3, -2
    z = z*i
  end do

  double_fact = z


  return
end function double_fact



subroutine print_dmat (iout, nmat, idim1, idim2, str, xmat)

! +----------------------------------------------------------------+
! |                                                                |
! | print_dmat  --  CAJH, 11.2012                                  |
! |                                                                |
! |                                                                |
! | Print an array, xmat, of nmat rectangular matrices (dimension  |
! | idim1 x idim2). That is, xmat = xmat(idim1,idim2,nmat).        |
! |                                                                |
! | The array xmat is a double precision array and the matrix is   |
! | printed to an output file (defined by iout) in a nice format.  |
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


end program linmol


