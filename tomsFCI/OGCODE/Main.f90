
      Program FCI
      Use Precision
      Use HamiltonianData
      Use IO
      Use DoIntegrals
      Use DIIS
      Use UHFWrap
      Use IntTrans
      Use IndexCI
      Use Hamiltonian
      Use CI
      Implicit None
! Dimensioning variables
      Integer :: NOccA, NOccB, NAO, NDet, NDetA, NDetB, NRoots
! AO Integrals
      Real (Kind=pr) :: ENuc
      Real (Kind=pr), Allocatable :: Olap(:,:), OneH(:,:)
      Real (Kind=pr), Allocatable :: ERI(:,:,:,:)
! MO Integrals
      Real (Kind=pr), Allocatable :: H1a(:,:)
      Real (Kind=pr), Allocatable :: H1b(:,:)
      Real (Kind=pr), Allocatable :: H2aa(:,:,:,:)
      Real (Kind=pr), Allocatable :: H2ab(:,:,:,:)
      Real (Kind=pr), Allocatable :: H2bb(:,:,:,:)
! SCF output
      Real (Kind=pr), Allocatable :: EvecA(:,:)
      Real (Kind=pr), Allocatable :: EvecB(:,:)
      Real (Kind=pr) :: ESCF
! CI output
      Real (Kind=pr), Allocatable :: CIVec(:,:), ECI(:)
! Error checking variables
      Integer :: IAlloc
! Timer stuf
      Integer :: RunTime(8)
      Character (len=10) :: Time
      Character (len=8)  :: Date
      Character (len=5)  :: Zone

!==============================================================!
!  This is an FCI code.                                        !
!                                                              !
!  While I store integrals in the full GHF spinorbital style,  !
!  this is actually a UHF-based FCI code only at present.      !
!--------------------------------------------------------------!
!  Read the basic information for the calculation.             !
!==============================================================!

      Call ReadInput(NOccA,NOccB,NAO,NRoots,ENuc)


!====================!
!  Allocate memory.  !
!====================!

      IAlloc = 0
      Allocate(Olap(NAO,NAO),          &
               OneH(NAO,NAO),          &
               ERI(NAO,NAO,NAO,NAO),   &
               EvecA(NAO,NAO),         &
               EvecB(NAO,NAO),         &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in main"


!=========================================================!
!  Do the one- and two-electron integrals, then the SCF.  !
!=========================================================!

      Call Drv1E(Olap,OneH,NAO)
      Call Drv2E(ERI,NAO)
      Call DoUHF(Olap,OneH,ERI,EvecA,EvecB,NAO,NOccA,NOccB,ENuc,ESCF)


!=====================================================================!
!  Allocate and build the MO integrals, deallocate the AO integrals.  !
!=====================================================================!

! Allocate the MO integrals
      Allocate(H1a(NAO,NAO),              &
               H1b(NAO,NAO),              &
               H2aa(NAO,NAO,NAO,NAO),     &
               H2ab(NAO,NAO,NAO,NAO),     &
               H2bb(NAO,NAO,NAO,NAO),     &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in main"

! Build them
      H1a = OneH;   Call IntTran2(H1a,EvecA,NAO)
      H1b = OneH;   Call IntTran2(H1b,EvecB,NAO)
      H2aa = ERI;   Call IntTran4(H2aa,EvecA,EvecA,NAO,.true.)
      H2bb = ERI;   Call IntTran4(H2bb,EvecB,EvecB,NAO,.true.)
      H2ab = ERI;   Call IntTran4(H2ab,EvecA,EvecB,NAO,.false.)

! Deallocate the AO integrals
      DeAllocate(Olap, OneH, ERI, EvecA, EvecB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in main"


!=======================================!
!  Construct the list of determinants.  !
!=======================================!

      Call Date_and_Time(Date,Time,Zone,RunTime)
      Call SetUpIndexCI(NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Call Timer(RunTime,"SetUpIndexCI")


!===========================!
!  Allocate the CI vector.  !
!===========================!

      Allocate(CIVec(NDet,NRoots), ECI(NRoots), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in main"


!=======================================!
!  Allocate and build the Hamiltonian.  !
!=======================================!

      Call Date_and_Time(Date,Time,Zone,RunTime)
      Call SetUpHam(NDetA,NDetB,NDet,NAO)
      Call BuildHam(H1a,H1b,H2aa,H2ab,H2bb,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      HDiag = HDiag + ENuc
      Call Timer(RunTime,"SetUpHam and BuildHam")


!==================!
!  Solve the FCI.  !
!==================!

      Call Date_and_Time(Date,Time,Zone,RunTime)
      Call SolveCI(H2aa,H2bb,H2ab,CIVec,ESCF,ECI,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NRoots)
      Call Timer(RunTime,"SolveCI")


!==============================================!
!  Lastly, deallocate memory and exit safely.  !
!==============================================!

      If(HamiltonianType == 0) Call ShutDownHamiltonianData
      Call ShutDownIndexCI
      Call ShutDownHam
      DeAllocate(H1a, H1b, H2aa, H2ab, H2bb, CIVec, ECI, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in main"

      Stop
      End Program FCI

