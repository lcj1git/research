

program gendet

  use constants
  use grid
  use spin
  use spatpr
  use purcart

! +----------------------------------------------------------------+
! |                                                                |
! | gendet  --  CAJH, 02.2013                                      |
! |                                                                |
! |                                                                |
! | A program to read the data file created by the program         |
! | phflat and save the determinantal expansion of the (symmetry   |
! | projected) HF state to an output file.                         |
! |                                                                |
! | NOTE: This code currently supports only a 1-det state!         |
! |                                                                |
! +----------------------------------------------------------------+
! |                                                                |
! | outfile should be read as                                      |
! |                                                                |
! |   read(unit)  ndet    [4-byte integer]                         |
! |   read(unit)  ci_vec  [complex*16, dim(ndet) vector]           |
! |   do k = 1, ndet                                               |
! |     read(unit)  det k                                          |
! |   end do                                                       |
! |                                                                |
! | All determinants are saved in GHF-sytle, as 2*nbs*(nup+ndn)    |
! | arrays of complex*16 entries.                                  |
! |                                                                |
! +----------------------------------------------------------------+

  implicit none

  integer :: nbas, nbct
  integer :: nup, ndn
  integer :: iwfnty, imethd
  integer :: ngrda, ngrdb, ngrdg, imult
  integer :: igrda, igrdb, igrdg, sptirr

  character(len=4) :: sptgr
  character(len=40) :: frot


  ! output files

  character(len=40) :: datfile
  character(len=40) :: outfile

  integer :: idat
  parameter ( idat = 11 )


  ! other variables

  integer :: norb, xdim, edim
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: npar, len

  integer :: ndspin, ndsptl, ndcplx
  integer :: ndet, nstat, mxdtst, istat
  integer, dimension(:), allocatable :: idetst

  real(kind=dp), dimension(:), allocatable :: stval
  complex(kind=dp), dimension(:,:), allocatable :: ovdet, hmdet
  complex(kind=dp), dimension(:,:), allocatable :: ovstat, hmstat, stvec
  complex(kind=dp), dimension(:,:), allocatable :: mdvec
  complex(kind=dp), dimension(:,:), allocatable :: darr

  complex(kind=dp), dimension(:,:), allocatable :: smat, xmat

  ! spatial symmetry projection

  integer :: nop

  complex(kind=dp), dimension(:,:,:), allocatable :: rotmat, wgtmat

  ! spin integration grid

  real(kind=dp), dimension(:), allocatable :: grda, grdb, grdg
  real(kind=dp), dimension(:), allocatable :: wgta, wgtb, wgtg


  ! read input parameters

  npar = command_argument_count ()

  if ( npar < 2 ) then
    write (6, *) 'error: Incorrect number of command line arguments.'
    stop
  end if

    len = 40
  call get_command_argument (1, datfile, len)
    len = 40
  call get_command_argument (2, outfile, len)
  if ( npar > 2 ) then
      len = 40
    call get_command_argument (3, frot, len)
  end if

  ! read data file

  call read_data

  istat = 1  ! state to use

  ! load pure to Cartesian transformation matrices

  ! if ( nbct > nbas )  call load_purcart


  ! build projected state

  call build_grid

  call build_proj_state (idat, outfile)

  ! if ( nbct > nbas )  call shutdown_purcart


  ! shutdown

  call shutdown_grid
  call shutdown_gendet


contains


subroutine read_data

! +----------------------------------------------------------------+
! |                                                                |
! | read_data  --  CAJH, 01.2013                                   |
! |                                                                |
! |                                                                |
! | Read all important variables from previous data file.          |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: j, k, ib, nflg
  logical :: lchk, ftest
  integer :: nosq, ndtt


  ! open data file

  inquire (file = datfile, exist = ftest)

  if ( .not. ftest ) then
    write (6, *) 'error: Data file ', datfile, ' does not exist.'
    stop
  end if

  open (unit = idat, file = datfile, status = 'old', form = 'unformatted')


  ! read several integer scalars

  read (idat) nbas, nbct, norb, nup, ndn, iwfnty, imethd

    ! some logical variables

    lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
    lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
    lspin = mod(imethd-1,2) == 1

    lsuhf = .false.
    lsghf = .false.

    if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
    if ( lspin .and. iwfnty == 3 )  lsghf = .true.

    nosq = norb*norb
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


  ! read spin projection parameters

  if ( mod(imethd-1,2) == 1 ) then
    read (idat) imult
    read (idat) ngrda, ngrdb, ngrdg
    read (idat) igrda, igrdb, igrdg
  else
    ngrda = 1
    ngrdb = 1
    ngrdg = 1
  end if


  ! read spatial projection parameters

  if ( mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3 ) then
    read (idat) sptgr
    read (idat) sptirr
  end if


  ! determine ndcplx

  if ( lcplx ) then
    ndcplx = 2
  else
    ndcplx = 1
  end if

  ! determine ndspin, ndsptl

  if ( lsghf ) then
    ndspin = imult
  else
    ndspin = 1
  end if

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
  else
    ndsptl = 1
  end if


  ! build S as identity matrix
  allocate (smat(nbas,nbas))
  smat = d0

  do ib = 1, nbas
    smat(ib,ib) = d1
  end do    

  ! read xmat

  allocate (xmat(nbas,norb))
  read (idat) ((xmat(j,k), j = 1, nbas), k = 1, norb)

  ! test whether X! . S . X = I

  call check_xmat (xmat, lchk)

  if ( .not. lchk ) then
    write (6, *) 'error: Test X!.S.X = I failed in read_data.'
    stop
  end if


  ! read previously stored number of determinants

  read (idat) ndet, nstat


  ! allocate idetst and fill

  allocate (idetst(nstat))
  read (idat) (idetst(k), k = 1, nstat)


  ! determine several important variables

  ndtt   = ndet*ndspin*ndsptl*ndcplx
  mxdtst = maxval (idetst(1:nstat))


  ! allocate important arrays

  allocate (ovdet(ndtt,ndtt), hmdet(ndtt,ndtt))
  allocate (ovstat(nstat,nstat), hmstat(nstat,nstat))
  allocate (stvec(nstat,nstat), stval(nstat))
  allocate (mdvec(ndspin*ndsptl*ndcplx*mxdtst,nstat))

  ! fill in arrays

  read (idat) ((ovdet(j,k), j = 1, ndtt), k = 1, ndtt)
  read (idat) ((hmdet(j,k), j = 1, ndtt), k = 1, ndtt)
  read (idat) ((ovstat(j,k), j = 1, nstat), k = 1, nstat)
  read (idat) ((hmstat(j,k), j = 1, nstat), k = 1, nstat)
  read (idat) ((stvec(j,k), j = 1, nstat), k = 1, nstat)
  read (idat) (stval(k), k = 1, nstat)

  ! CI expansion of states among its dets

  read (idat) ((mdvec(j,k), j = 1, ndspin*ndsptl*ndcplx*mxdtst), &
             & k = 1, nstat)


  ! read array of determinants

  allocate (darr(xdim,ndet))

  do k = 1, ndet
    read (idat) (darr(j,k), j = 1, xdim)
  end do


  ! read nflag and decide whether to proceed

  read (idat) nflg

  ! if ( nflg /= 0 ) then
  !   write (6, *) 'error: Aborting calculation.'
  !   write (6, *) 'It appears like the previous calculation did not converge.'
  !   stop
  ! end if


  ! close data file

  close (unit = idat)


  return
end subroutine read_data



subroutine build_grid

! +----------------------------------------------------------------+
! |                                                                |
! | build_grid  --  CAJH, 02.2013                                  |
! |                                                                |
! |                                                                |
! | Construct the spin projection grid as well as the rotation     |
! | matrices necessary for spatial symmetry projection.            |
! |                                                                |
! +----------------------------------------------------------------+

  ! other variables

  integer :: j

  ! constants

  real(kind=dp) :: pi
  parameter ( pi = 4.0e0_dp * atan(d1) )


  if ( lsptl ) then
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
        write (6, *) 'error: Point group not supported.'
    end select
  else
    nop = 1

    allocate (wgtmat(1,1,1))
    allocate (rotmat(norb,norb,1))

    wgtmat(1,1,1) = z0
    rotmat(1:norb,1:norb,1) = z0
  end if

  ! prepare spin projection grid

  allocate (grda(ngrda), grdb(ngrdb), grdg(ngrdg))
  allocate (wgta(ngrda), wgtb(ngrdb), wgtg(ngrdg))

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

  return
end subroutine build_grid



subroutine shutdown_grid

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_grid  --  CAJH, 02.2013                               |
! |                                                                |
! |                                                                |
! | Deallocate memory for arrays built in build_grid.              |
! |                                                                |
! +----------------------------------------------------------------+


  ! shutdown spin projection grid

  deallocate (grda, grdb, grdg)
  deallocate (wgta, wgtb, wgtg)


  ! shutdown spatial projection grid

  deallocate (wgtmat, rotmat)


  return
end subroutine shutdown_grid



subroutine build_proj_state (iout, outfile)

! +----------------------------------------------------------------+
! |                                                                |
! | build_proj_state  --  CAJH, 02.2013                            |
! |                                                                |
! |                                                                |
! | Subroutine to construct the determinant expansion of the       |
! | (symmetry projected) HF state and save it to the output file   |
! | outfile.                                                       |
! |                                                                |
! +----------------------------------------------------------------+

  ! input variables

  integer, intent(in) :: iout
  character(len=*), intent(in) :: outfile


  ! other variables

  integer :: idum
  integer :: nosq, ntot, iwfnt1
  integer :: ngrds, nevl, ndtt
  logical :: lcplx, lsptl, lspin, lsuhf, lsghf

  integer :: jval, mval, nval
  integer :: j, j1, j1a, j1b
  integer :: icnt, k1, k2
  integer :: is1, ik1, indz
  integer :: ind1, ind2, ind3, ind4
  integer :: aind, bind, gind
  logical :: lint, ftest

  real(kind=dp) :: alpha, beta, gamma, wgtY
  complex(kind=dp) :: wgtX

  integer, dimension(:,:), allocatable :: indx

  complex(kind=dp), dimension(:), allocatable :: cicf

  complex(kind=dp), dimension(:), allocatable :: dst, dstX, dstY


  lcplx = mod(imethd-1,8) >= 4 .and. mod(imethd-1,8) <= 7
  lsptl = mod(imethd-1,4) >= 2 .and. mod(imethd-1,4) <= 3
  lspin = mod(imethd-1,2) == 1

  lsuhf = .false.
  lsghf = .false.

  if ( lspin .and. iwfnty == 2 )  lsuhf = .true.
  if ( lspin .and. iwfnty == 3 )  lsghf = .true.


  ! set some important variables

  nosq = norb*norb
  ntot = nup + ndn

  iwfnt1 = iwfnty
  if ( lsuhf ) iwfnt1 = 3

  ngrds = ngrda*ngrdb*ngrdg
  nevl  = ndcplx*ndet*nop*ngrds
  ndtt  = nevl*ndspin*ndsptl


  ! decide spin to use in projection

  if ( lspin .and. mod(imult,2) == 0 ) then
    jval = imult-1
    lint = .false.
  else if ( lspin ) then
    jval = (imult-1)/2
    lint = .true.
  end if


  ! prepare indexing array indx:
  !   first index contains the state index of the det
  !   second index contains the index of det in state

  allocate (indx(ndet,2))

  icnt = 1

  do k1 = 1, nstat
    do k2 = 1, idetst(k1)
      indx(icnt,1) = k1
      indx(icnt,2) = k2

      icnt = icnt + 1
    end do
  end do


  ! allocate matrix of CI coefficients

  allocate (cicf(ndtt))


  ! allocate dst

  if ( iwfnt1 == 1 ) then
    allocate (dst(nosq))
  else if ( iwfnt1 == 2 ) then
    allocate (dst(2*nosq))
  else if ( iwfnt1 == 3 ) then
    allocate (dst(4*nosq))
  end if

  ! some scratch versions of dst

  if ( lsptl ) then
    if ( iwfnty == 1 ) then
      allocate (dstY(nosq))
    else if ( iwfnty == 2 ) then
      allocate (dstY(2*nosq))
    else if ( iwfnty == 3 ) then
      allocate (dstY(4*nosq))
    end if
  end if

  if ( lspin ) then
    allocate (dstX(4*nosq))
  end if


  ! open file, save ndet

  inquire (file = outfile, exist = ftest)
    
  if ( ftest ) then
    write (6, *) 'error: File ', outfile, ' exists in current directory.'
    stop
  end if
    
  open (unit = iout, file = outfile, status = 'new', form = 'unformatted')

  write (iout) ndtt
  !jag20
  write (*,*) ndtt


  ! loop over determinants
  ! do two passes:
  !   - record CI coefficients in first pass
  !   - store determinants in second pass

    aind = 0
    bind = 0
    gind = 0

  do j = 1, nevl

    ind1 = mod(j-1,ngrds) + 1                 ! spin grid index
    ind2 = mod((j-1)/ngrds,nop) + 1           ! spatial symm index
    ind3 = mod((j-1)/(ngrds*nop),ndcplx) + 1  ! complex left index
    ind4 = (j-1)/(ngrds*nop*ndcplx) + 1       ! det index

    if ( lspin ) then
      gind = mod(ind1-1,ngrdg) + 1
      bind = mod((ind1-1)/ngrdg,ngrdb) + 1
      aind = (ind1-1)/(ngrdg*ngrdb) + 1

      alpha = grda(aind)
      beta  = grdb(bind)
      gamma = grdg(gind)
    end if

    ! loop over ndspin, ndsptl

    is1 = indx(ind4,1)

    do j1 = 1, ndspin*ndsptl
      j1a = mod(j1-1,ndspin) + 1
      j1b = (j1-1)/ndspin + 1

      ik1  = (indx(ind4,2)-1)*ndspin*ndsptl*ndcplx &
         & + (ind3-1)*ndspin*ndsptl + j1

      indz = (j-1)*ndspin*ndsptl*ndcplx &
         & + (ind3-1)*ndspin*ndsptl + j1

      if ( lsghf .and. .not. lint ) then
        mval = jval
        nval = jval-2*(j1a-1)
      else if ( lsghf ) then
        mval = jval
        nval = jval-(j1a-1)
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
        wgtX = wgtX * wgtmat(1,j1b,ind2)
      end if

      wgtX = wgtX * stvec(is1,istat) * mdvec(ik1,is1)

      cicf(indz) = wgtX
    end do
  end do

  write (iout) cicf

  do j = 1, nevl

    ind1 = mod(j-1,ngrds) + 1                 ! spin grid index
    ind2 = mod((j-1)/ngrds,nop) + 1           ! spatial symm index
    ind3 = mod((j-1)/(ngrds*nop),ndcplx) + 1  ! complex index
    ind4 = (j-1)/(ngrds*nop*ndcplx) + 1       ! det index

    if ( lspin ) then
      gind = mod(ind1-1,ngrdg) + 1
      bind = mod((ind1-1)/ngrdg,ngrdb) + 1
      aind = (ind1-1)/(ngrdg*ngrdb) + 1

      alpha = grda(aind)
      beta  = grdb(bind)
      gamma = grdg(gind)
    end if


    ! load determinants (we actually load D*)
    ! multiply by spatial rotation matrices

    if ( iwfnty == 1 ) then
      if ( .not. lsptl ) then
        if ( ind3 == 1 ) then
          dst(1:nosq) = conjg (darr(1:nosq,ind4))
        else
          dst(1:nosq) = darr(1:nosq,ind4)
        end if

      else
        if ( ind3 == 1 ) then
          dstY(1:nosq) = conjg (darr(1:nosq,ind4))
        else
          dstY(1:nosq) = darr(1:nosq,ind4)
        end if

        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY, norb, z0, dst, norb)
      end if

    else if ( iwfnty == 2 ) then
      if ( .not. lsptl .and. .not. lspin ) then
        if ( ind3 == 1 ) then
          dst(1:2*nosq) = conjg (darr(1:2*nosq,ind4))
        else
          dst(1:2*nosq) = darr(1:2*nosq,ind4)
        end if

      else if ( .not. lsptl .and. lspin ) then
        if ( ind3 == 1 ) then
          dstX(1:2*nosq) = conjg (darr(1:2*nosq,ind4))
        else
          dstX(1:2*nosq) = darr(1:2*nosq,ind4)
        end if

      else
        if ( ind3 == 1 ) then
          dstY(1:2*nosq) = conjg (darr(1:2*nosq,ind4))
        else
          dstY(1:2*nosq) = darr(1:2*nosq,ind4)
        end if

        if ( lspin ) then
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(1), norb, z0, dstX(1), norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(nosq+1), norb, z0, dstX(nosq+1), norb)
        else
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(1), norb, z0, dst(1), norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(nosq+1), norb, z0, dst(nosq+1), norb)
        end if
      end if

    else if ( iwfnty == 3 ) then
      if ( .not. lsptl .and. .not. lspin ) then
        if ( ind3 == 1 ) then
          dst(1:4*nosq) = conjg (darr(1:4*nosq,ind4))
        else
          dst(1:4*nosq) = darr(1:4*nosq,ind4)
        end if

      else if ( .not. lsptl .and. lspin ) then
        if ( ind3 == 1 ) then
          dstX(1:4*nosq) = conjg (darr(1:4*nosq,ind4))
        else
          dstX(1:4*nosq) = darr(1:4*nosq,ind4)
        end if

      else
        if ( ind3 == 1 ) then
          dstY(1:4*nosq) = conjg (darr(1:4*nosq,ind4))
        else
          dstY(1:4*nosq) = darr(1:4*nosq,ind4)
        end if

        if ( lspin ) then
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(1), 2*norb, z0, dstX(1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(norb+1), 2*norb, z0, dstX(norb+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(2*nosq+1), 2*norb, z0, dstX(2*nosq+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(2*nosq+norb+1), 2*norb, z0, dstX(2*nosq+norb+1), 2*norb)

        else
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(1), 2*norb, z0, dst(1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(norb+1), 2*norb, z0, dst(norb+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(2*nosq+1), 2*norb, z0, dst(2*nosq+1), 2*norb)
        call zgemm ('n', 'n', norb, norb, norb, z1, rotmat(1,1,ind2), norb, &
             & dstY(2*nosq+norb+1), 2*norb, z0, dst(2*nosq+norb+1), 2*norb)
        end if
      end if
    end if


    ! form full spin-orbital matrices for SUHF

    if ( lsuhf ) then
      call uhf_to_ghf (norb, nup, ndn, dstX)
    end if


    ! multiply D by spin rotation matrix

    if ( lspin ) then
      call spinrot (norb, alpha, beta, gamma, dstX, dst, 1)
    end if

    do j1 = 1, ndspin*ndsptl
      call save_det (iout, iwfnt1, dst)
    end do
  end do


  ! close outfile

  close (unit = iout)


  ! clear memory

  deallocate (indx)
  deallocate (cicf)
  deallocate (dst)

  if ( lsptl ) then
    deallocate (dstY)
  end if

  if ( lspin ) then
    deallocate (dstX)
  end if


  return
end subroutine build_proj_state



subroutine save_det (iux, iwfnX, dmat)

! +----------------------------------------------------------------+
! |                                                                |
! | save_det  --  CAJH, 01.2013                                    |
! |                                                                |
! |                                                                |
! | Save determinant dmat into unit iux.                           |
! |                                                                |
! | Note that the determinant is saved in the std AO basis for     |
! | proper interaction with guess_read or other programs. This is  |
! | accomplished by letting                                        |
! |                                                                |
! |    D  <-  X . D,                                               |
! |                                                                |
! | where X is the transformation matrix [ = S^(-1/2) ].           |
! |                                                                |
! +----------------------------------------------------------------+

  integer, intent(in) :: iux, iwfnX

  complex(kind=dp), dimension(:), intent(in) :: dmat

  integer :: xdimx, nosq

  complex(kind=dp), dimension(:), allocatable :: dscr

  complex(kind=dp), dimension(:,:), allocatable :: scr_s1, scr_s2
  complex(kind=dp), dimension(:,:), allocatable :: dmat_s1, dmat_s2


  nosq = norb*norb

  xdimx = 0
  if ( iwfnX == 1 ) then
    xdimx = nbas*norb
  else if ( iwfnX == 2 ) then
    xdimx = 2*nbas*norb
  else if ( iwfnX == 3 ) then
    xdimx = 4*nbas*norb
  end if


  ! allocate space for D in std AO basis

  allocate (dscr(xdimx))

  ! scratch arrays for GHF wfns

  if ( iwfnX == 3 ) then
    allocate (dmat_s1(2*nbas,2*norb), dmat_s2(2*norb,2*norb), &
            & scr_s1(nbas,norb), scr_s2(norb,norb))
  else
    allocate (dmat_s1(1,1), dmat_s2(1,1), &
            & scr_s1(1,1), scr_s2(1,1))
  end if


  ! transform to std AO basis
  ! perform  D  <-  X . D

  if ( iwfnX == 1 ) then
    call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
         & dmat(1), norb, z0, dscr, nbas)

  else if ( iwfnX == 2 ) then
    call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
         & dmat(1), norb, z0, dscr(1), nbas)
    call zgemm ('n', 'n', nbas, norb, norb, z1, xmat, nbas, &
         & dmat(nosq+1), norb, z0, dscr(nbas*norb+1), nbas)

  else if ( iwfnX == 3 ) then
    dmat_s2(1:2*norb,1:2*norb) = &
         & reshape (dmat(1:4*nosq), (/ 2*norb, 2*norb /))

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

  write (iux) dscr


  deallocate (dscr)
  deallocate (dmat_s1, dmat_s2, scr_s1, scr_s2)


  return
end subroutine save_det



subroutine shutdown_gendet

! +----------------------------------------------------------------+
! |                                                                |
! | shutdown_gendet  --  CAJH, 02.2013                             |
! |                                                                |
! |                                                                |
! | Deallocate all scratch arrays.                                 |
! |                                                                |
! +----------------------------------------------------------------+

  deallocate (idetst)
  deallocate (ovdet, hmdet)
  deallocate (ovstat, hmstat)
  deallocate (stvec, stval)
  deallocate (mdvec)
  deallocate (darr)
  deallocate (smat, xmat)


  return
end subroutine shutdown_gendet



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


  allocate (scr1(nbas,norb), scr2(norb,norb))


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

  deallocate (scr1, scr2)


  return
end subroutine check_xmat



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


end program gendet


