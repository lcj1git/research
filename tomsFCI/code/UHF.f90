
   Module UHF
   Use Precision
   Use Constants
   Use DIIS
   Private
   Public  :: DrvUHF, UHFMix, SimpleUHF
   Contains

      Subroutine DrvUHF(Olap,OneH,ERI,EvecsA,EvecsB,ENuc,ESCF,NOccA,NOccB,NAO,UseDIIS,GuessType)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NOccB, NAO, GuessType
      Logical,        Intent(In)    :: UseDIIS
      Real (Kind=pr), Intent(In)    :: Olap(NAO,NAO)
      Real (Kind=pr), Intent(In)    :: OneH(NAO,NAO)
      Real (Kind=pr), Intent(In)    :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(InOut) :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Intent(In)    :: ENuc
      Real (Kind=pr), Intent(Out)   :: ESCF
      Real (Kind=pr), Parameter     :: TolMax1 = 1.0E-08_pr
      Real (Kind=pr), Parameter     :: TolMax2 = 1.0E-08_pr
      Real (Kind=pr), Allocatable   :: Trans(:,:)
      Real (Kind=pr), Allocatable   :: Scr1(:,:), Scr2(:,:)
      Real (Kind=pr), Allocatable   :: FockA(:,:), DenMatA(:,:), EvalsA(:)
      Real (Kind=pr), Allocatable   :: FockB(:,:), DenMatB(:,:), EvalsB(:)
      Real (Kind=pr), Allocatable   :: FockADIIS(:,:), FockAVec(:,:,:)
      Real (Kind=pr), Allocatable   :: FockBDIIS(:,:), FockBVec(:,:,:)
      Real (Kind=pr), Parameter     :: LevelShift = 0.0_pr
      Integer           :: IAlloc(13)
      Real (Kind=pr)    :: RMS
      Integer           :: NIter

!===================================================!
!  This routine does open-shell SCF calculations.   !
!  At the end of NewDen, DenHF is the new density.  !
!===================================================!

      Open(7,File='Output',Position='Append')
      Allocate(Scr1(NAO,NAO),             Stat=IAlloc(1))
      Allocate(Scr2(NAO,NAO),             Stat=IAlloc(2))
      Allocate(FockA(NAO,NAO),            Stat=IAlloc(3))
      Allocate(DenMatA(NAO,NAO),          Stat=IAlloc(4))
      Allocate(EvalsA(NAO),               Stat=IAlloc(5))
      Allocate(FockB(NAO,NAO),            Stat=IAlloc(6))
      Allocate(DenMatB(NAO,NAO),          Stat=IAlloc(7))
      Allocate(EvalsB(NAO),               Stat=IAlloc(8))
      Allocate(FockADIIS(NAO,NAO),        Stat=IAlloc(9))
      Allocate(FockBDIIS(NAO,NAO),        Stat=IAlloc(10))
      Allocate(FockAVec(NAO,NAO,NDIIS),   Stat=IAlloc(11))
      Allocate(FockBVec(NAO,NAO,NDIIS),   Stat=IAlloc(12))
      Allocate(Trans(NAO,NAO),            Stat=IAlloc(13))
      If(Any(IAlloc /= 0)) Stop 'Could not Allocate in UHF'


!=========================================================!
!  Get the transformation matrix and prepare for output.  !
!=========================================================!

      Write(7,*)
      Write(7,1000)
      Write(7,1010)
      Write(7,1020)
      Write(7,1030)
      Call GetTrans(Trans,Olap,NAO)
 

!===============!
!  Do the SCF.  !
!===============!

! Initial Guess - read in density matrix!
      NIter = 0
      Select Case(GuessType)
        Case(0) ! Core guess
          DenMatA = Zero
          DenMatB = Zero
        Case(1) ! Read orbitals and build density
          DenMatA = Zero
          DenMatB = Zero
          Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,RMS,NIter,NAO,Scr1,Scr2)
        Case(2) ! Fancy Hubbard guess
          Call GuessDMat(NAO,DenMatA,.true.)
          Call GuessDMat(NAO,DenMatB,.false.)
          RMS = Max(MaxVal(Abs(DenMatA)),MaxVal(Abs(DenMatB)))
        Case Default
          Stop "Unknown Guess Type!"
      End Select
      Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
      Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
      Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
      Write(7,1040) ESCF,NIter,RMS


! Initialize DIIS
      FockADIIS = FockA
      FockBDIIS = FockB
      FockAVec  = Zero
      FockBVec  = Zero
      ErrVecs   = Zero


! Iterate to convergence
      Do
        Call FockDiag(FockADIIS,EvecsA,EvalsA,Trans,NAO,Scr1,Scr2)
        Call FockDiag(FockBDIIS,EvecsB,EvalsB,Trans,NAO,Scr1,Scr2)
        Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,   &
                    RMS,NIter,NAO,Scr1,Scr2)
        Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
        Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
        Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
        Write(7,1040) ESCF,NIter,RMS
        If(RMS < TolMax1) Exit
        Call DrvDIIS(Olap,Trans,FockA,FockB,FockADIIS,FockBDIIS,FockAVec,FockBVec,  &
                     DenMatA,DenMatB,Scr1,Scr2,RMS,NAO,NIter)
        If(.not. UseDIIS) Then
          FockADIIS = FockA
          FockBDIIS = FockB
        End If
      End Do


! TMH, experimental - turn off DIIS once we hit "convergence"
      Do
        Call FockDiag(FockA,EvecsA,EvalsA,Trans,NAO,Scr1,Scr2)
        Call FockDiag(FockB,EvecsB,EvalsB,Trans,NAO,Scr1,Scr2)
        Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,   &
                    RMS,NIter,NAO,Scr1,Scr2)
        Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
        Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
        Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
        Write(7,1040) ESCF,NIter,RMS
        If(RMS < TolMax2) Exit
      End Do

      Write(7,1020)
      Write(7,1050) NIter
      Write(7,1060) ESCF
      Write(7,1000)
      Close(7)


!===============================!
!  Deallocate and exit safely.  !
!===============================!

      DeAllocate(Scr1,        Stat=IAlloc(1))
      DeAllocate(Scr2,        Stat=IAlloc(2))
      DeAllocate(FockA,       Stat=IAlloc(3))
      DeAllocate(DenMatA,     Stat=IAlloc(4))
      DeAllocate(EvalsA,      Stat=IAlloc(5))
      DeAllocate(FockB,       Stat=IAlloc(6))
      DeAllocate(DenMatB,     Stat=IAlloc(7))
      DeAllocate(EvalsB,      Stat=IAlloc(8))
      DeAllocate(FockADIIS,   Stat=IAlloc(9))
      DeAllocate(FockBDIIS,   Stat=IAlloc(10))
      DeAllocate(FockAVec,    Stat=IAlloc(11))
      DeAllocate(FockBVec,    Stat=IAlloc(12))
      DeAllocate(Trans,       Stat=IAlloc(13))
      If(Any(IAlloc /= 0)) Stop 'Could not DeAllocate in UHF'

1000  Format(14x,'**************************************************')
1010  Format(14X,'*               UHF summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'*   UHF Energy     Iteration    Biggest change   *')
1040  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* SCF has converged in ',I4,' iterations',11x,'*')
1060  Format(14X,'* Final UHF Energy is ',2x,F18.12,' a.u.',2x,'*')

      Return
      End Subroutine DrvUHF






      Subroutine UHFMix(EvecsA,EvecsB,NAO,NOccA,NOccB)
      Implicit None
      Integer,        Intent(In)    :: NoccA, NOccB, NAO
      Real (Kind=pr), Intent(InOut) :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Allocatable :: V1(:), V2(:)
      Real (Kind=pr), Parameter :: Theta = Pi/12
      Real (Kind=pr), Parameter :: C = Cos(Theta), S = Sin(Theta)

!==================================================!
!  Mixes alpha HOMO and LUMO, beta HOMO and LUMO.  !
!==================================================!

      Allocate(V1(NAO), V2(NAO))
! Do alpha
      V1 = EvecsA(:,NOccA)
      V2 = EvecsA(:,NOccA+1)
      EvecsA(:,NOccA)   = C*V1 + S*V2
      EvecsA(:,NOccA+1) = C*V2 - S*V1
! Do beta
      If(NOccB == 0) GoTo 10
      V1 = EvecsB(:,NOccB)
      V2 = EvecsB(:,NOccB+1)
      EvecsB(:,NOccB)   = C*V1 - S*V2
      EvecsB(:,NOccB+1) = C*V2 + S*V1
10    DeAllocate(V1, V2)

      Return
      End Subroutine UHFMix






      Subroutine SimpleUHF(OneH,ERI,EvecsA,EvecsB,ENuc,ESCF,NOccA,NOccB,NAO)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB,NAO
      Real (Kind=pr), Intent(In)  :: OneH(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: ENuc
      Real (Kind=pr), Intent(Out) :: ESCF
      Real (Kind=pr), Allocatable :: PA(:,:), PB(:,:), FA(:,:), FB(:,:)
      Integer :: IAlloc

!======================================================================!
!  Presuming we already know the orbitals, just calculate the energy.  !
!======================================================================!

! Reuse existing routines
      Allocate(PA(NAO,NAO), PB(NAO,NAO), FA(NAO,NAO), FB(NAO,NAO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in SimpleUHF"


! Build the density matrices PA and PB using FA and FB for scratch
      Call NewDen(PA,PB,EvecsA,EvecsB,NOccA,NOccB,ESCF,IAlloc,NAO,FA,FB)

! Build the Fock matrices FA and FB
      Call MkFock(OneH,ERI,PA,PB,FA,NAO)
      Call MkFock(OneH,ERI,PB,PA,FB,NAO)

! Get the SCF Energy
      Call SCFEnergy(OneH,FA,FB,PA,PB,ENuc,ESCF,NAO)

! Write to disk
      Open(7,File="Output",Position="Append")
        Write(7,*)
        Write(7,1000)
        Write(7,1010)
        Write(7,1020)
        Write(7,1060) ESCF
        Write(7,1000)
      Close(7)

! Deallocate and exit safely
      DeAllocate(FA, FB, PA, PB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in SimpleUHF"
1000  Format(14x,'**************************************************')
1010  Format(14X,'*               UHF summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1060  Format(14X,'* Final UHF Energy is ',2x,F18.12,' a.u.',2x,'*')

      Return
      End Subroutine SimpleUHF






      Subroutine GetTrans(Trans,Olap,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: Olap(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Trans(NAO,NAO)
      Real (Kind=pr), Allocatable :: UMat(:,:), SVals(:), Scr(:,:)
      Real (Kind=pr), Parameter   :: TolMax = 1.0E-9_pr
      Integer :: I

!=================================================!
!  This generates the canonically orthogonalized  !
!  transformation matrix, X = U*S^(-1/2).         !
!=================================================!

! Get U, s
      Allocate(UMat(NAO,NAO), SVals(NAO), Scr(NAO,NAO))
      Call DiagR(Olap,SVals,UMat,NAO)

! Form s^(-1/2)
      Scr = Zero
      Do I=1,NAO
       If(SVals(I) <= TolMax) Stop 'Linear dependence in basis set'
       Scr(I,I) = One/Sqrt(SVals(I))
      End Do

! Build X, which is block-diagonal
      Trans = MatMul(UMat,Scr)
      DeAllocate(UMat,SVals,Scr)

      Return
      End Subroutine GetTrans






      Subroutine NewDen(DenA,DenB,EvecsA,EvecsB,NOccA,NOccB,     &
                        RMS,NIter,NAO,DenHF2,DiffM)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NoccB, NAO
      Integer,        Intent(InOut) :: NIter
      Real (Kind=pr), Intent(In)    :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Intent(InOut) :: DenA(NAO,NAO),   DenB(NAO,NAO)
      Real (Kind=pr), Intent(Out)   :: RMS
      Real (Kind=pr) :: DenHF2(NAO,NAO), DiffM(NAO,NAO)
      Integer, Parameter :: SCFMaxCyc = 15000

!=======================================================================!
!  This generates the density matrix from the SCF eigenvectors, forms   !
!  the convergence criterion (I'm using largest absolute change in the  !
!  density matrix elements) and updates the density matrix.             !
!=======================================================================!

! Do the alpha MOs
      DiffM  = Transpose(EvecsA)
      DenHF2 = MatMul(EvecsA(:,1:NOccA),DiffM(1:NOccA,:))
      DiffM  = DenHF2-DenA
      DenA   = DenHF2
      RMS    = MaxVal(Abs(DiffM))

! Do the beta MOs
      If(NOccB == 0) GoTo 10
      DiffM  = Transpose(EvecsB)
      DenHF2 = MatMul(EvecsB(:,1:NOccB),DiffM(1:NOccB,:))
      DiffM  = DenHF2-DenB
      DenB   = DenHF2
      RMS    = Max(MaxVal(Abs(DiffM)),RMS)

! Check convergence and out!
10    NIter = NIter + 1
      If(NIter >= SCFMaxCyc) Stop 'SCF did not converge'

      Return
      End Subroutine NewDen






      Subroutine MKFock(OneH,ERI,Den1,Den2,Fock,NAO)
      Implicit None
      Integer,           Intent(In)  :: NAO
      Real (Kind=pr),    Intent(In)  :: OneH(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr),    Intent(In)  :: Den1(NAO,NAO), Den2(NAO,NAO)
      Real (Kind=pr),    Intent(Out) :: Fock(NAO,NAO)
      Integer :: Mu, Nu, Lam, Kap

!=================================================!
!  This forms the Fock matrix, obviously enough.  !
!=================================================!

      Fock = OneH
      Do Mu = 1,NAO
      Do Nu = Mu,NAO
        Do Lam = 1,NAO
        Do Kap = 1,NAO
          Fock(Mu,Nu) = Fock(Mu,Nu)                                            &
                      + ERI(Mu,Nu,Lam,Kap)*(Den1(Lam,Kap) + Den2(Lam,Kap))     &
                      - ERI(Mu,Kap,Lam,Nu)*Den1(Lam,Kap)
        End Do   
        End Do   
        Fock(Nu,Mu) = Fock(Mu,Nu)
      End Do
      End Do

      Return
      End Subroutine MKFock






      Subroutine SCFEnergy(OneH,FockA,FockB,DenA,DenB,ENuc,ESCF,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: OneH(NAO*NAO)
      Real (Kind=pr), Intent(In)  :: FockA(NAO*NAO), DenA(NAO*NAO)
      Real (Kind=pr), Intent(In)  :: FockB(NAO*NAO), DenB(NAO*NAO)
      Real (Kind=pr), Intent(In)  :: ENuc
      Real (Kind=pr), Intent(Out) :: ESCF

!======================================!
!  This gets the total energy from     !
!  E = 1/2*P(i,j) * [F(i,j) + H(i,j)]  ! 
!======================================!

      ESCF = Dot_Product(DenA,FockA+OneH)/2   &
           + Dot_Product(DenB,FockB+OneH)/2   &
           + ENuc

      Return
      End Subroutine SCFEnergy






      Subroutine FockDiag(Fock,Evec,Evals,Trans,NAO,FockMO,Scr)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: Trans(NAO,NAO), Fock(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Evec(NAO,NAO), Evals(NAO)
      Real (Kind=pr) :: FockMO(NAO,NAO), Scr(NAO,NAO)

!==================================================================!
!  This generates the MO basis Fock matrix, FockMO                 !
!  We diagonalize it in EIG (leaving eigenvalues on the diagonal)  !
!  In this EIG call, the eigenvectors are put in Scr               !
!  The real eigenvectors are then made by another transformation   ! 
!==================================================================!

      FockMO = MatMul(Fock,Trans)
      Scr    = Transpose(Trans)
      FockMO = MatMul(Scr,FockMO)
      Call DiagR(FockMO,Evals,Scr,NAO)
      Evec = MatMul(Trans,Scr)

      Return
      End Subroutine FockDiag






      Subroutine DrvDIIS(Olap,Trans,FockA,FockB,FockADIIS,FockBDIIS,FockAVec,FockBVec,  &
                         DenMatA,DenMatB,Scr1,Scr2,RMS,NAO,NIter)
      Implicit None
      Integer,        Intent(In)  :: NAO, NIter
      Real (Kind=pr), Intent(In)  :: RMS
      Real (Kind=pr), Intent(In)  :: Olap(NAO,NAO),  Trans(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FockA(NAO,NAO), DenMatA(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FockB(NAO,NAO), DenMatB(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: FockADIIS(NAO,NAO), FockBDIIS(NAO,NAO)
      Real (Kind=pr) :: FockAVec(NAO,NAO,NDIIS), FockBVec(NAO,NAO,NDIIS)
      Real (Kind=pr) :: Scr1(NAO,NAO), Scr2(NAO,NAO)
! Local variables
      Integer :: I
      Real (Kind=pr) :: Coeffs(NDIIS)

!====================================!
!  This subroutine drives the DIIS.  !
!====================================!

! Step 1: Cycle the list of error vectors and Fock matrices
      Do I = 1,NDIIS-1
        ErrVecs(:,I)    = ErrVecs(:,I+1)
        FockAVec(:,:,I) = FockAVec(:,:,I+1)
        FockBVec(:,:,I) = FockBVec(:,:,I+1)
      End Do

! Step 2: Append to the list of Fock matrices and error vectors
      FockAVec(:,:,NDIIS) = FockA
      FockBVec(:,:,NDIIS) = FockB
      Call GetErrVec(ErrVecs(:,NDIIS),Olap,Trans,FockA,FockB,DenMatA,DenMatB,Scr1,Scr2,NAO)

! Step 3: Construct FockDIIS, the Fock matrix to be diagonalized next
      If(NIter >= StartDIIS) Then
        Call DoDIIS(Coeffs)
        FockADIIS = Zero
        FockBDIIS = Zero
        Do I = 1,NDIIS
          FockADIIS = FockADIIS + FockAVec(:,:,I)*Coeffs(I)
          FockBDIIS = FockBDIIS + FockBVec(:,:,I)*Coeffs(I)
        End Do
      Else
        FockADIIS = FockA           ! By default, we diagonalize the new Fock matrix
        FockBDIIS = FockB           ! By default, we diagonalize the new Fock matrix
      End If

      Return
      End Subroutine DrvDIIS






      Subroutine GetErrVec(Vec,S,X,FA,FB,PA,PB,Scr1,Scr2,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: S(NAO,NAO), X(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FA(NAO,NAO), FB(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: PA(NAO,NAO), PB(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Vec(2*NAO*NAO)
      Real (Kind=pr) :: Scr1(NAO,NAO), Scr2(NAO,NAO)

!=======================================================!
!  The error vector is the real part of the commutator  !
!  done in an orthogonalized basis                      !
!=======================================================!

! The commutators are FPS - SPF = C.  Do alpha first.
      Scr1 = MatMul(FA,PA)
      Scr2 = MatMul(Scr1,S)             ! C = F.P.S
      Scr1 = MatMul(PA,FA)
      Scr2 = Scr2 - MatMul(S,Scr1)      ! C = F.P.S - S.P.F
      Scr1 = MatMul(Scr2,X)             ! C.X
      Scr2 = MatMul(Transpose(X),Scr1)  ! X+.C.X
      Vec(1:NAO*NAO) = ReShape(Scr2, (/ NAO*NAO /))

! Do beta.
      Scr1 = MatMul(FB,PB)
      Scr2 = MatMul(Scr1,S)             ! C = F.P.S
      Scr1 = MatMul(PB,FB)
      Scr2 = Scr2 - MatMul(S,Scr1)      ! C = F.P.S - S.P.F
      Scr1 = MatMul(Scr2,X)             ! C.X
      Scr2 = MatMul(Transpose(X),Scr1)  ! X+.C.X
      Vec(NAO*NAO+1:2*NAO*NAO) = ReShape(Scr2, (/NAO*NAO /))

      Return
      End Subroutine GetErrVec






      Subroutine GuessDMat(N,DenMat,Alpha)
      Use Precision
      Use HamiltonianData
      Implicit None
      Integer,        Intent(In)  :: N
      Logical,        Intent(In)  :: Alpha
      Real (Kind=pr), Intent(Out) :: DenMat(N,N)
      Integer :: Indx, IX, IY, I
      Real (Kind=pr) :: X, V, Z

!=================================================================!
!  Forms an initial guess density matrix which is just diagonal.  !
!  Define y = U/(U+T) = Z/(Z+1).  Define two occupation numbers   !
!    x   = (1+y)/2                                                !
!    1-x = (1-y)/2                                                !
!  For large Z, x and 1-x go to 1 and 0.                          !
!  For small Z, x and 1-x go to 1/2.                              !
!                                                                 !
!  We define the occupation number in site 1,1 to be x.           !
!  We have                                                        !
!    n(IX,IY+1) = n(IX+1,IY) = -n(IX,IY).                         !
!  Beta occupation numbers are just 1 - Alpha occupations.        !
!=================================================================!

! Define an index function
      Indx(IX,IY) = IX + (IY-1)*NSitesX


! Exit if we get here and are not doing Hubbard
      If(HamiltonianType /= 2) Stop "Cannot use Hubbard guess for this Hamiltonian!"


! Define the occupation numbers, switching for the beta density matrix.
      Z = HubbardU/HubbardT
      X = (Z + F12)/(Z +1)
      If(.not. Alpha) X = 1-X


! Build the density matrix.  V is the occupation number at site 1,1.
! We switch across a row by flipping the sign of V.
! We switch across a column by noting 
      DenMat = Zero
      Do IX = 1,NSitesX
        Select Case(Mod(IX,2))
          Case(1);  V = X     ! In rows 1, 3, 5, we start with occupation x
          Case(0);  V = 1-X   ! In rows 2, 4, 6, we start with occupation 1-x
        End Select
        Do IY = 1,NSitesY
          I = Indx(IX,IY)
          DenMat(I,I) = V
          V = 1-V             ! As we go across a row, we switch x and 1-x
        End Do
      End Do

      Return
      End Subroutine GuessDMat

   End Module UHF

