
   Module DoIntegrals
   Use Precision
   Use ObaraSaika
   Use HamiltonianData
   Implicit None
   Private
   Public :: Drv1E, Drv2E
   Contains

      Subroutine Drv1E(Olap,OneH,NAO)
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: Olap(NAO,NAO), OneH(NAO,NAO)

!===========================================================!
!  This just calls the appropriate one-electron integrals.  !
!===========================================================!

      Select Case(HamiltonianType)
        Case(0)
          Call Drv1EMolecular(Olap,OneH,NAO,NAtom)
        Case(1)
          Call Drv1EPairing(Olap,OneH,NAO)
        Case(2)
          Call Drv1EHubbard(Olap,OneH,NSitesX,NSitesY,NAO,DoPBCX,DoPBCY)
          OneH = OneH*HubbardT
        Case Default
          Stop "Unrecognized HamiltonianType in Drv1E"
      End Select

      Return
      End Subroutine Drv1E






      Subroutine Drv2E(ERI,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: ERI(NAO,NAO,NAO,NAO)

!===========================================================!
!  This just calls the appropriate two-electron integrals.  !
!===========================================================!

      Select Case(HamiltonianType)
        Case(0)
          Call Drv2EMolecular(ERI,NAO)
        Case(1)
          Call Drv2EPairing(ERI,NAO)
          ERI = -ERI*PairingG
        Case(2)
          Call Drv2EHubbard(ERI,NAO)
          ERI = ERI*HubbardU
        Case Default
          Stop "Unrecognized HamiltonianType in Drv2E"
      End Select

      Return
      End Subroutine Drv2E






      Subroutine Drv1EMolecular(Olap,OneH,NAO,NAtom)
      Implicit None
      Integer,        Intent(In)  :: NAO, NAtom
      Real (Kind=pr), Intent(Out) :: Olap(NAO,NAO), OneH(NAO,NAO)
! Basis function indices
      Integer        :: Mu, Nu, IMu, INu, IAtom
! Centers
      Real (Kind=pr) ::  A(3),  B(3),  C(3)
! Angular momentum vectors
      Integer        :: LA(3), LB(3)
! Normalization constants
      Real (Kind=pr) :: XMu, XNu
! Exponents
      Real (Kind=pr) :: ZA, ZB
! Contraction coefficients and nuclear charge
      Real (Kind=pr) :: CA, CB, Q
! Integral values
      Real (Kind=pr) :: ValS, ValV, ValT, SMuNu, HMuNu

!=========================================!
!  One-electron integrals for molecules.  !
!=========================================!

      Olap = Zero
      OneH = Zero
      Do Mu = 1,NAO
      Do Nu = Mu,NAO
! Initialize the contracted integrals for this basis function pair
        SMuNu = Zero
        HMuNu = Zero
! Set the contracted data for Obara-Saika
        A   = AOCenter(Mu,:)
        B   = AOCenter(Nu,:)
        LA  = LMNAO(Mu,:)
        LB  = LMNAO(Nu,:)
        XMu = AONorm(Mu)
        XNu = AONorm(Nu)
! Loop over primitive functions
        Do IMu = FirstPrimInCont(Mu), LastPrimInCont(Mu)
        Do INu = FirstPrimInCont(Nu), LastPrimInCont(Nu)
! Set the primitive data for Obara-Saika
          ZA = AOExp(IMu)
          ZB = AOExp(INu)
          CA = AOCoeff(IMu)
          CB = AOCoeff(INu)
! Calculate primitive integrals
          ValS = Olap2c(ZA,CA,A,LA,ZB,CB,B,LB)
          ValT = TInt(ZA,CA,A,LA,ZB,CB,B,LB)
          ValV = Zero
          Do IAtom = 1,NAtom
            C = AtomCenter(IAtom,:)
            Q = AtomCharge(IAtom)
            ValV = ValV + VInt(ZA,CA,A,LA,ZB,CB,B,LB,Q,C)
          End Do
! Increment the contracted integrals
          SMuNu = SMuNu + ValS*XMu*XNu
          HMuNu = HMuNu + (ValT - ValV)*XMu*XNu
        End Do
        End Do
        Olap(Mu,Nu) = SMuNu
        Olap(Nu,Mu) = SMuNu
        OneH(Mu,Nu) = HMuNu
        OneH(Nu,Mu) = HMuNu
      End Do
      End Do

      Return
      End Subroutine Drv1EMolecular






      Subroutine Drv1EPairing(Olap,OneH,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: Olap(NAO,NAO), OneH(NAO,NAO)
      Integer :: Mu

!=======================================================!
!  One-electron integrals for the pairing Hamiltonian.  !
!=======================================================!

      Olap = Zero
      OneH = Zero
      Do Mu = 1,NAO
        OneH(Mu,Mu) = Mu
        Olap(Mu,Mu) = One
      End Do

      Return
      End Subroutine Drv1EPairing






      Subroutine Drv1EHubbard(Olap,OneH,NX,NY,NAO,DoPBCX,DoPBCY)
      Implicit None
      Integer,        Intent(In)  :: NX, NY, NAO
      Logical,        Intent(In)  :: DoPBCX, DoPBCY
      Real (Kind=pr), Intent(Out) :: Olap(NAO,NAO), OneH(NAO,NAO)
      Integer :: IX, IY, I, JX, JY, J, Indx

!=======================================================!
!  One-electron integrals for the Hubbard Hamiltonian.  !
!=======================================================!

      Indx(IX,IY) = IX + (IY-1)*NX

      Olap = Zero
      OneH = Zero
      Do I = 1,NAO
        Olap(I,I) = One
      End Do

      Do IX = 1,NX
      Do IY = 1,NY
        I = Indx(IX,IY)

! Set T with the site to the right (with IX + 1)
! At IX = NX, we interact with 1 if the system is periodic in X and NX /= 1
        If(IX /= NX) Then
          JX = IX+1
        Else If(DoPBCX  .and.  NX /= 1) Then
          JX = 1
        Else
          GoTo 10
        End If
        J = Indx(JX,IY)
        OneH(I,J) = -1
        OneH(J,I) = -1
        If(DoPBCX .and. NX == 2 .and. JX == 1) Then
          OneH(I,J) = -2
          OneH(J,I) = -2
        End If


! Set T with the site to the bottom (with IY + 1)
! At IY = NY, we interact with 1 if the system is peridoic in Y and NY /= 1
! At IY = NY, we interact with 1 instead if NY /= 1.
10      If(IY /= NY) Then
          JY = IY+1
        Else If(DoPBCY  .and.  NY /= 1) Then
          JY = 1
        Else
          Cycle
        End If
        J = Indx(IX,JY)
        OneH(I,J) = -1
        OneH(J,I) = -1
        If(DoPBCY .and. NY == 2 .and. JY == 1) Then
          OneH(I,J) = -2
          OneH(J,I) = -2
        End If
      End Do
      End Do

      Return
      End Subroutine Drv1EHubbard






      Subroutine Drv2EMolecular(ERI,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: ERI(NAO,NAO,NAO,NAO)
! Contracted function indices
      Integer        :: Mu, Nu, Lam, Kap
      Integer        :: Iota
! Primitive function indices
      Integer        :: IMu, INu, ILam, IKap
! Centers
      Real (Kind=pr) ::  A(3),  B(3),  C(3),  D(3)
! Angular momentum vectors
      Integer        :: LA(3), LB(3), LC(3), LD(3)
! Normalization constants
      Real (Kind=pr) :: XA, XB, XC, XD
! Exponents
      Real (Kind=pr) :: ZA, ZB, ZC, ZD
! Contraction coefficients
      Real (Kind=pr) :: CA, CB, CC, CD
! Integral values
      Real (Kind=pr) :: Val, VMuNuLamKap
   
!=========================================!
!  Two-electron integrals for molecules.  !
!=========================================!

      Iota(Lam,Mu,Nu) = Lam + Lam/Mu * (Nu-Lam)

      ERI = Zero
      Do Mu = 1,NAO
      Do Nu = 1,Mu
      Do Lam = 1,Mu
      Do Kap = 1,Iota(Lam,Mu,Nu)
! Initialize the contracted integral for this basis function quartet
        VMuNuLamKap = Zero
! Set the contracted data for Obara-Saika
        A  = AOCenter(Mu,:)
        B  = AOCenter(Nu,:)
        C  = AOCenter(Lam,:)
        D  = AOCenter(Kap,:)
        LA = LMNAO(Mu,:)
        LB = LMNAO(Nu,:)
        LC = LMNAO(Lam,:)
        LD = LMNAO(Kap,:)
        XA = AONorm(Mu)
        XB = AONorm(Nu)
        XC = AONorm(Lam)
        XD = AONorm(Kap)
! Loop over primitive functions
        Do IMu  = FirstPrimInCont(Mu),  LastPrimInCont(Mu)
        Do INu  = FirstPrimInCont(Nu),  LastPrimInCont(Nu)
        Do ILam = FirstPrimInCont(Lam), LastPrimInCont(Lam)
        Do IKap = FirstPrimInCont(Kap), LastPrimInCont(Kap)
! Set the primitive data for Obara-Saika
          ZA = AOExp(IMu)
          ZB = AOExp(INu)
          ZC = AOExp(ILam)
          ZD = AOExp(IKap)
          CA = AOCoeff(IMu)
          CB = AOCoeff(INu)
          CC = AOCoeff(ILam)
          CD = AOCoeff(IKap)
! Calculate the primitive integral and increment the contracted integral
          Val = ERInt(ZA,CA,A,LA,ZB,CB,B,LB,ZC,CC,C,LC,ZD,CD,D,LD)
          VMuNuLamKap = VMuNuLamKap + Val*XA*XB*XC*XD
        End Do
        End Do
        End Do
        End Do
        ERI(Mu,Nu,Lam,Kap) = VMuNuLamKap
        ERI(Nu,Mu,Lam,Kap) = VMuNuLamKap
        ERI(Mu,Nu,Kap,Lam) = VMuNuLamKap
        ERI(Nu,Mu,Kap,Lam) = VMuNuLamKap
        ERI(Lam,Kap,Mu,Nu) = VMuNuLamKap
        ERI(Lam,Kap,Nu,Mu) = VMuNuLamKap
        ERI(Kap,Lam,Mu,Nu) = VMuNuLamKap
        ERI(Kap,Lam,Nu,Mu) = VMuNuLamKap
      End Do
      End Do
      End Do
      End Do

      Return
      End Subroutine Drv2EMolecular






      Subroutine Drv2EPairing(ERI,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: ERI(NAO,NAO,NAO,NAO)
      Integer :: Mu, Nu
   
!=======================================================!
!  Two-electron integrals for the pairing Hamiltonian.  !
!=======================================================!

      ERI = Zero
      Do Mu = 1,NAO
      Do Nu = 1,NAO
        ERI(Mu,Nu,Mu,Nu) = One
      End Do
      End Do

      Return
      End Subroutine Drv2EPairing






      Subroutine Drv2EHubbard(ERI,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(Out) :: ERI(NAO,NAO,NAO,NAO)
      Integer :: I

!=======================================================!
!  Two-electron integrals for the Hubbard Hamiltonian.  !
!=======================================================!

      ERI = Zero
      Do I = 1,NAO
        ERI(I,I,I,I) = One
      End Do

      Return
      End Subroutine Drv2EHubbard

   End Module DoIntegrals

