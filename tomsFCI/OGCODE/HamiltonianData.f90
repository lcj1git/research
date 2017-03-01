
   Module HamiltonianData
   Use Precision
   Use Constants
   Implicit None
! This flag just sets which kind of Hamiltonian we have.
! HamiltonianType = (0,1,2) => Molecular/Pairing/Hubbard
   Integer :: HamiltonianType
! This stuff is needed for the pairing Hamiltonian
   Real (Kind=pr) :: PairingG
! This stuff is needed for the Hubbard Hamiltonian
   Integer :: NSitesX, NSitesY
   Logical :: DoPBCX, DoPBCY
   Real (Kind=pr) :: HubbardT, HubbardU
! Thus stuff is needed for the electronc Hamiltonian
   Integer :: NPrim, NAtom
   Character (len=2),  Dimension(:),   Allocatable :: AtomName
   Real (Kind=pr),     Dimension(:,:), Allocatable :: AtomCenter
   Real (Kind=pr),     Dimension(:),   Allocatable :: AtomCharge
   Character (len=12), Dimension(:),   Allocatable :: AtomBasis
   Real (Kind=pr),     Dimension(:),   Allocatable :: AOExp
   Real (Kind=pr),     Dimension(:),   Allocatable :: AOCoeff
   Real (Kind=pr),     Dimension(:),   Allocatable :: AONorm
   Real (Kind=pr),     Dimension(:,:), Allocatable :: AOCenter
   Integer,            Dimension(:,:), Allocatable :: LMNAO
   Integer,            Dimension(:),   Allocatable :: IAOAtom
   Integer,            Dimension(:),   Allocatable :: FirstPrimInCont
   Integer,            Dimension(:),   Allocatable :: LastPrimInCont
   Public
   Contains

      Subroutine SetUpGeometry(NAtom)
      Implicit None
      Integer, Intent(In) :: NAtom
      Integer :: IAlloc(4)
      Allocate(AtomName(NAtom),     Stat=IAlloc(1))
      Allocate(AtomCenter(NAtom,3), Stat=IAlloc(2))
      Allocate(AtomCharge(NAtom),   Stat=IAlloc(3))
      Allocate(AtomBasis(NAtom),    Stat=IAlloc(4))
      If(Any(IAlloc /= 0)) Stop 'Could not allocate atomic data'
      Return
      End Subroutine SetUpGeometry






      Subroutine SetUpBasis(NAO)
      Implicit None
      Integer, Intent(In) :: NAO
      Integer :: IAlloc(8)
      Allocate(AOExp(NPrim),           Stat=IAlloc(1))
      Allocate(AOCoeff(NPrim),         Stat=IAlloc(2))
      Allocate(AONorm(NAO),            Stat=IAlloc(3))
      Allocate(AOCenter(NAO,3),        Stat=IAlloc(4))
      Allocate(LMNAO(NAO,3),           Stat=IAlloc(5))
      Allocate(IAOAtom(NAO),           Stat=IAlloc(6))
      Allocate(FirstPrimInCont(NAO),   Stat=IAlloc(7))
      Allocate(LastPrimInCont(NAO),    Stat=IAlloc(8))
      If(Any(IAlloc /= 0)) Stop 'Could not allocate basis data'
      Return
      End Subroutine SetUpBasis






      Subroutine ReOrderD(NAO)
      Implicit None
      Integer, Intent(In) :: NAO
      Integer :: I, J
      Integer :: LVecOld(3), LVecNew(3)
      Integer :: LMe(3,6), LGau(3,6)
      Integer, Parameter, Dimension(18) ::         &
        LVecMe  = (/2,0,0,    1,1,0,   1,0,1,      &
                    0,2,0,    0,1,1,   0,0,2/),    &
        LVecGau = (/2,0,0,    0,2,0,   0,0,2,      &
                    1,1,0,    1,0,1,   0,1,1/)

!==================================================!
!  This subroutine reorders the d orbitals.        !
!  My order:         dxx, dxy, dxz, dyy, dyz, dzz  !
!  Gaussian's order: dxx, dyy, dzz, dxy, dxz, dyz  !
!==================================================!

! Start by filling in LMe and LGau.
      LMe  = ReShape(LVecMe,(/3,6/))
      LGau = ReShape(LVecGau,(/3,6/))

! Fix the functions.
! We have to fix the angular momentum vector, obviously
! The normalization of the primitives changes, which goes into contraction coefficients
! The normalization of the contracted function could also change
      Do I = 1,NAO
        If(Sum(LMNAO(I,:)) == 2  .and. LMNAO(I,1) == 2) Then   ! It's a dxx.  Start from here
          Do J = 1,6
            LVecOld = LMe(:,J)                     ! The angular momentum vector from my code
            LVecNew = LGau(:,J)                    ! The angular momentum vector from Gaussian
            LMNAO(I+J-1,:) = LVecNew               ! Fix the angular momtnum vector
            Call FixCoeff(LVecOld,LVecNew,I+J-1)   ! Fix the contraction coefficients
            Call FixNorm(I+J-1)                    ! Renormalize the wave function
          End Do
        End If
      End Do

      Return
      End Subroutine ReOrderD






      Subroutine ReOrderF(NAO)
      Implicit None
      Integer, Intent(In) :: NAO
      Integer :: I, J
      Integer :: LVecOld(3), LVecNew(3)
      Integer :: LMe(3,10), LGau(3,10)
      Integer, Parameter, Dimension(30) ::                        &
        LVecMe  = (/3,0,0,    2,1,0,   2,0,1,   1,2,0,   1,1,1,   &
                    1,0,2,    0,3,0,   0,2,1,   0,1,2,   0,0,3/), &
        LVecGau = (/3,0,0,    0,3,0,   0,0,3,   1,2,0,   2,1,0,   &
                    2,0,1,    1,0,2,   0,1,2,   0,2,1,   1,1,1/)

!================================================!
!  This subroutine reorders the f orbitals.      !
!  My order:         xxx, xxy, xxz, xyy, xyz,    !
!                    xzz, yyy, yyz, yzz, zzz     !
!  Gaussian's order: xxx, yyy, zzz, xyy, xxy,    !
!                    xxz, xzz, yzz, yyz, xyz     !
!================================================!

! Start by filling in LMe and LGau.
      LMe  = ReShape(LVecMe,(/3,10/))
      LGau = ReShape(LVecGau,(/3,10/))

! Fix the functions.
! We have to fix the angular momentum vector, obviously
! The normalization of the primitives changes, which goes into contraction coefficients
! The normalization of the contracted function could also change
      Do I = 1,NAO
        If(Sum(LMNAO(I,:)) == 3  .and. LMNAO(I,1) == 3) Then   ! It's a fxxx.  Start from here.
          Do J = 1,10
            LVecOld = LMe(:,J)                     ! The angular momentum vector from my code
            LVecNew = LGau(:,J)                    ! The angular momentum vector from Gaussian
            LMNAO(I+J-1,:) = LVecNew               ! Fix the angular momtnum vector
            Call FixCoeff(LVecOld,LVecNew,I+J-1)   ! Fix the contraction coefficients
            Call FixNorm(I+J-1)                    ! Renormalize the wave function
          End Do
        End If
      End Do
 
      Return
      End Subroutine ReOrderF






      Subroutine FixCoeff(LVecOld,LVecNew,IAO)
      Implicit None
      Integer :: LVecOld(3), LVecNew(3), IAO
      Integer :: J
      Integer :: LOld, MOld, NOld
      Integer :: LNew, MNew, NNew
      Real (Kind=pr) :: A, DNrmNew, DNrmOld, F1

      LOld = LVecOld(1)
      MOld = LVecOld(2)
      NOld = LVecOld(3)
      LNew = LVecNew(1)
      MNew = LVecNew(2)
      NNew = LVecNew(3)
      Do J = FirstPrimInCont(IAO), LastPrimInCont(IAO)
        A = AOExp(J)
        DNrmNew = F1(2*LNew,2*A)*F1(2*MNew,2*A)*F1(2*NNew,2*A)
        DNrmOld = F1(2*LOld,2*A)*F1(2*MOld,2*A)*F1(2*NOld,2*A)
        AOCoeff(J) = AOCoeff(J)*Sqrt(DNrmOld/DNrmNew)
      End Do

      Return
      End Subroutine FixCoeff






      Subroutine FixNorm(IAO)
      Implicit None
      Integer :: IAO
      Integer :: L, M, N, I, J, IStart, IEnd
      Real (Kind=pr) :: A, C, DNorm, F1

      DNorm = Zero
      IStart = FirstPrimInCont(IAO)
      IEnd   = LastPrimInCont(IAO)
      L      = LMNAO(IAO,1)
      M      = LMNAO(IAO,2)
      N      = LMNAO(IAO,3)
      Do I = IStart,IEnd
      Do J = IStart,IEnd
        A = (AOExp(I)   + AOExp(J))/Two
        C = AOCoeff(I) * AOCoeff(J)
        DNorm = DNorm + C*F1(2*L,2*A)*F1(2*M,2*A)*F1(2*N,2*A)
      End Do
      End Do
      AONorm(IAO) = One/Sqrt(DNorm)

      Return
      End Subroutine FixNorm






      Subroutine ShutDownHamiltonianData
      Implicit None
      Integer :: IAlloc(12)
      DeAllocate(AtomName,          Stat=IAlloc(1))
      DeAllocate(AtomCenter,        Stat=IAlloc(2))
      DeAllocate(AtomCharge,        Stat=IAlloc(3))
      DeAllocate(AtomBasis,         Stat=IAlloc(4))
      DeAllocate(AOExp,             Stat=IAlloc(5))
      DeAllocate(AOCoeff,           Stat=IAlloc(6))
      DeAllocate(AONorm,            Stat=IAlloc(7))
      DeAllocate(AOCenter,          Stat=IAlloc(8))
      DeAllocate(LMNAO,             Stat=IAlloc(9))
      DeAllocate(IAOAtom,           Stat=IAlloc(10))
      DeAllocate(FirstPrimInCont,   Stat=IAlloc(11))
      DeAllocate(LastPrimInCont,    Stat=IAlloc(12))
      If(Any(IAlloc /= 0)) Stop 'Could not deallocate basis data'
      Return
      End Subroutine ShutDownHamiltonianData

   End Module HamiltonianData

