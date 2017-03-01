
   Module Hamiltonian
   Use Precision
   Use Constants
   Use IndexCI
   Implicit None
   Real (Kind=pr), Allocatable :: HDiag(:)
   Real (Kind=pr), Allocatable :: FockSSa(:,:,:)
   Real (Kind=pr), Allocatable :: FockSSb(:,:,:)
   Real (Kind=pr), Allocatable :: FockOSa(:,:,:)
   Real (Kind=pr), Allocatable :: FockOSb(:,:,:)
   Real (Kind=pr), Allocatable :: HDoubA(:,:)
   Real (Kind=pr), Allocatable :: HDoubB(:,:)
   Real (Kind=pr), Allocatable :: FockA(:,:)
   Real (Kind=pr), Allocatable :: FockB(:,:)
   Contains

!=====================================================================!
!  This stores most everything we need to build the FCI Hamiltonian.  !
!                                                                     !
!  Diagonal elements are stored in HDiag.                             !
!                                                                     !
!  Single excitations are complicated.  Suppose I have a determinant  !
!     |IDet> = |IDetA IDetB>                                          !
!  and I want to do Pa! Qa |IDet>.  First of all, I know that         !
!     Pa! Qa |IDetA> = |JDetA>                                        !
!  complete with a possible sign, and I can read all of that from     !
!  IndPQonA(P,Q,IDetA).  Let's keep track explicitly of the sign in   !
!  these notes; we'll call it #.  So                                  !
!     Pa! Qa |IDet> = # |JDetA IDetB> = # |JDet>                      !
!  where I can extract |JDet> from IndAB2Full.                        !
                                                                      !
!  Now, this means that                                               !
!     <Pa! Qa IDet|H|IDet> = # <JDet|H|IDet>                          !
!  and accordingly                                                    !
!     <JDet|H|IDet> = # <Pa|F(IDetA,IDetB)|Qa>                        !
!  Break that Fock operator into same spin and opposite spin pieces:  !
!     Fss(IDetA,P,Q) = <P|ha|Q> + Sum_j ERIaa(P,J,Q,J)                !
!     Fos(IDetB,P,Q) = Sum_j ERIab(P,J,Q,J)                           !
!  where we note that the sum over J in the same-spin part goes over  !
!  J occupied in IDetA and the sum over J in the opposite-spin part   !
!  goes over J occupied in IDetB.  My final result is                 !
!     <JDet|H|IDet> = # Fss(IDetA,P,Q) + # Fos(IDetB,P,Q)             !
!  So for every alpha determinant IDetA I store its same-spin and     !
!  it's opposite-spin components in FockSSa and FockOSa, and the      !
!  same of course holds for beta.  These components can be combined   !
!  to build the single excitation blocks of the Hamiltonian.          !
!                                                                     !
!  Making progress!  The same-spin double excitations are easy.       !
!     Pa! Qa! Sa Ra |IDetA IDetB> = # |JDetA IDetB>                   !
!  which means that                                                   !
!     <JDetA IDetB|H|IDetA IDetB> = # Vaa(P,Q,R,S)                    !
!  The opposite-spin double excitations are a touch worse:            !
!     Pa! Sa |IDetA> = #a |JDetA>                                     !
!     Qb! Rb |IDetB> = #b |JDetB>                                     !
!     <JDetA JDetB|H|IDetA IDetB> = #a #b Vab(P,Q,S,R)                !
!  We'll store the same-spin doubles and compute the opposite-spin    !
!  doubles on the fly.                                                !
!=====================================================================!

      Subroutine SetUpHam(NDetA,NDetB,NDet,NAO)
      Implicit None
      Integer, Intent(In) :: NDetA, NDetB, NDet, NAO
      Integer :: IAlloc

!============================================!
!  Fills in space for Hamiltonian matrices.  !
!============================================!

      Allocate(HDiag(NDet),                &
               FockSSa(NDetA,NAO,NAO),     &
               FockSSb(NDetB,NAO,NAO),     &
               FockOSa(NDetB,NAO,NAO),     &
               FockOSb(NDetA,NAO,NAO),     &
               FockA(NAO,NAO),             &
               FockB(NAO,NAO),             &
               HDoubA(NDetA,NDoubA),       &
               HDoubB(NDetB,NDoubB),       &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in SetUpHam"

      Return
      End Subroutine SetUpHam






      Subroutine BuildHam(H1a,H1b,H2aa,H2ab,H2bb,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: H1a(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H1b(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2aa(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2ab(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2bb(NAO,NAO,NAO,NAO)
      Integer,        Allocatable :: OccsA(:), OccsB(:)
      Real (Kind=pr) :: HElement, HElement1, HElement2
      Integer :: IDet, JDet, KDet, IDetA, IDetB, JDetA, JDetB, KDetA, KDetB
      Integer :: IDoubA, IDoubB, Sgn, Sgn1, Sgn2, IAlloc
      Integer :: I, J, P, Q, R, S, II, JJ, AA, BB, A, B

!=========================================================================!
!  Builds as much of the Hamiltonian as we can afford to save in core.    !
!=========================================================================!

! Allocate some local space
      Allocate(OccsA(NOccA), OccsB(NOccB), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in BuildHam"


! BuildHDiag
      Do IDetA = 1,NDetA; OccsA = OccListA(IDetA,:)
      Do IDetB = 1,NDetB; OccsB = OccListB(IDetB,:)
        Call GetHDiag(HElement,OccsA,OccsB,H1a,H1b,H2aa,H2ab,H2bb,NOccA,NOccB,NAO)
        IDet = IndAB2Full(IDetA,IDetB)
        HDiag(IDet) = HElement
      End Do
      End Do


! Build FockSSa.  This is the same-spin     part of Fa for each determinant
! Build FockOSb.  This is the opposite-spin part of Fb for each determinant
      Do IDet = 1,NDetA
        OccsA = OccListA(IDet,:)
        Do P = 1,NAO
        Do Q = 1,NAO
          HElement1 = H1a(P,Q)
          HElement2 = Zero
          Do J = 1,NOccA
            I = OccsA(J)
            HElement1 = HElement1 + H2aa(P,I,Q,I)
            HElement2 = HElement2 + H2ab(I,P,I,Q)
          End Do
          FockSSa(IDet,P,Q) = HElement1
          FockOSb(IDet,P,Q) = HElement2
        End Do
        End Do
      End Do


! Build FockSSb.  This is really the same-spin     part of Fb for each determinant
! Build FockOSa.  This is really the opposite-spin part of Fa for each determinant
      Do IDet = 1,NDetB
        OccsB = OccListB(IDet,:)
        Do P = 1,NAO
        Do Q = 1,NAO
          HElement1 = H1b(P,Q)
          HElement2 = Zero
          Do J = 1,NOccB
            I = OccsB(J)
            HElement1 = HElement1 + H2bb(P,I,Q,I)
            HElement2 = HElement2 + H2ab(P,I,Q,I)
          End Do
          FockSSb(IDet,P,Q) = HElement1
          FockOSa(IDet,P,Q) = HElement2
        End Do
        End Do
      End Do


! Build HDoubA
      HDoubA = Zero
      Do IDetA = 1,NDetA
        IDoubA = 1
        Do II = 1,NOccA
        Do JJ = II+1,NOccA
        Do AA = 1,NVrtA
        Do BB = AA+1,NVrtA
          I = OccListA(IDetA,II)
          J = OccListA(IDetA,JJ)
          A = VrtListA(IDetA,AA)
          B = VrtListA(IDetA,BB)
          JDetA = IndPQonA(A,I,IDetA)
          Sgn1 = Sign(1,JDetA)
          JDetA = JDetA*Sgn1
          KDetA = IndPQonA(B,J,JDetA)
          Sgn2 = Sign(1,KDetA)
          KDetA = KDetA*Sgn2
          HElement = Sgn1*Sgn2*H2aa(A,B,I,J)
          HDoubA(IDetA,IDoubA) = HElement
          IDoubA = IDoubA + 1
        End Do
        End Do
        End Do
        End Do
      End Do


! Build HDoubB
      HDoubB = Zero
      Do IDetB = 1,NDetB
        IDoubB = 1
        Do II = 1,NOccB
        Do JJ = II+1,NOccB
        Do AA = 1,NVrtB
        Do BB = AA+1,NVrtB
          I = OccListB(IDetB,II)
          J = OccListB(IDetB,JJ)
          A = VrtListB(IDetB,AA)
          B = VrtListB(IDetB,BB)
          JDetB = IndPQonB(A,I,IDetB)
          Sgn1 = Sign(1,JDetB)
          JDetB = JDetB*Sgn1
          KDetB = IndPQonB(B,J,JDetB)
          Sgn2 = Sign(1,KDetB)
          KDetB = KDetB*Sgn2
          HElement = Sgn1*Sgn2*H2bb(A,B,I,J)
          HDoubB(IDetB,IDoubB) = HElement
          IDoubB = IDoubB + 1
        End Do
        End Do
        End Do
        End Do
      End Do


! Deallocate and exit safely
      DeAllocate(OccsA, OccsB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in BuildHam"

      Return
      End Subroutine BuildHam






      Subroutine GetHDiag(HElement,OccsA,OccsB,H1a,H1b,H2aa,H2ab,H2bb,NOccA,NOccB,NAO)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO
      Integer,        Intent(In)  :: OccsA(NOccA), OccsB(NOccB)
      Real (Kind=pr), Intent(In)  :: H1a(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H1b(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2aa(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2ab(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2bb(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(Out) :: HElement
      Integer :: I, J, IOcc, JOcc

!========================================================================!
!  Get the energy of a determinant with occupied orbitals from OccList.  !
!========================================================================!

      HElement = Zero

! Get the alpha part of the energy
      Do IOcc = 1,NOccA
        I = OccsA(IOcc)
        HElement = HElement + H1a(I,I)
        Do JOcc = 1,NOccA
          J = OccsA(JOcc)
          HElement = HElement + F12*H2aa(I,J,I,J)
        End Do
      End Do


! Get the beta part of the energy
      Do IOcc = 1,NOccB
        I = OccsB(IOcc)
        HElement = HElement + H1b(I,I)
        Do JOcc = 1,NOccB
          J = OccsB(JOcc)
          HElement = HElement + F12*H2bb(I,J,I,J)
        End Do
      End Do


! Get the mixed spin part
      Do IOcc = 1,NOccA
      Do JOcc = 1,NOccB
        I = OccsA(IOcc)
        J = OccsB(JOcc)
        HElement = HElement + H2ab(I,J,I,J)
      End Do
      End Do

      Return
      End Subroutine GetHDiag






      Subroutine BuildHC(Psi,HPsi,H2aa,H2bb,H2ab,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: Psi(NDet)
      Real (Kind=pr), Intent(In)  :: H2aa(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2bb(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2ab(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(Out) :: HPsi(NDet)
      Real (Kind=pr) :: CICoeff, HElement
      Integer :: IDet, IDetA, IDetB, JDet, JDetA, JDetB, KDet, KDetA, KDetB
      Integer :: IDoubA, IDoubB
      Integer :: A, B, I, J, II, JJ, AA, BB, Sgn, Sgn1, Sgn2

!==================================================!
!  This takes |Psi> and returns |HPsi> = H |Psi>.  !
!==================================================!

! Step 1: Handle the diagonal part
      HPsi = HDiag*Psi


! Step 2: Loop over excitations IDet and read off IDetA, IDetB
      Do IDet = 1,NDet
        CICoeff = Psi(IDet)
        IDetA = IndFull2AB(IDet,1)
        IDetB = IndFull2AB(IDet,2)

! Build FockA and FockB for IDet
        FockA = FockSSa(IDetA,:,:) + FockOSa(IDetB,:,:)
        FockB = FockSSb(IDetB,:,:) + FockOSb(IDetA,:,:)

! Build the alpha single excitations
        Do II = 1,NOccA
        Do AA = 1,NVrtA
          I = OccListA(IDetA,II)
          A = VrtListA(IDetA,AA)
          JDetA = IndPQonA(A,I,IDetA)
          Sgn = Sign(1,JDetA)
          JDetA = JDetA*Sgn
          JDet = IndAB2Full(JDetA,IDetB)
          HPsi(JDet) = HPsi(JDet) + Sgn*CICoeff*FockA(A,I)

! We'll slip the mixed-spin doubles into here
          Do JJ = 1,NOccB
          Do BB = 1,NVrtB
            J = OccListB(IDetB,JJ)
            B = VrtListB(IDetB,BB)
            JDetB = IndPQonB(B,J,IDetB)
            Sgn2 = Sign(1,JDetB)
            JDetB = JDetB*Sgn2
            JDet = IndAB2Full(JDetA,JDetB)
            HPsi(JDet) = HPsi(JDet) + Sgn*Sgn2*CICoeff*H2ab(A,B,I,J)
          End Do
          End Do
        End Do
        End Do

! Build the beta single excitations
        Do II = 1,NOccB
        Do AA = 1,NVrtB
          I = OccListB(IDetB,II)
          A = VrtListB(IDetB,AA)
          JDetB = IndPQonB(A,I,IDetB)
          Sgn = Sign(1,JDetB)
          JDetB = JDetB*Sgn
          JDet = IndAB2Full(IDetA,JDetB)
          HPsi(JDet) = HPsi(JDet) + Sgn*CICoeff*FockB(A,I)
        End Do
        End Do
      End Do

! Handle the same-spin double excitations
      Do IDetA = 1,NDetA
        Do IDoubA = 1,NDoubA
          JDetA = IndDoubA(IDetA,IDoubA)
          HElement = HDoubA(IDetA,IDoubA)
          Do IDetB = 1,NDetB
            IDet = IndAB2Full(IDetA,IDetB)
            JDet = IndAB2Full(JDetA,IDetB)
            HPsi(JDet) = HPsi(JDet) + Psi(IDet)*HElement
          End Do
        End Do
      End Do
      Do IDetB = 1,NDetB
        Do IDoubB = 1,NDoubB
          JDetB = IndDoubB(IDetB,IDoubB)
          HElement = HDoubB(IDetB,IDoubB)
          Do IDetA = 1,NDetA
            IDet = IndAB2Full(IDetA,IDetB)
            JDet = IndAB2Full(IDetA,JDetB)
            HPsi(JDet) = HPSi(JDet) + Psi(IDet)*Helement
          End Do
        End Do
      End Do
      Return

! TMH: This is archive code that does the same-spin doubles without resort to IndDoub or HDoub
!      If, after testing, this seems fast enough, then we can eliminate those guys and just use this code
!   

! Step 3: Handle the aa double excitations
      Do IDetA = 1,NDetA
        Do II = 1,NOccA
        Do JJ = II+1,NOccA
        Do AA = 1,NVrtA
        Do BB = AA+1,NVrtA
          I = OccListA(IDetA,II)
          J = OccListA(IDetA,JJ)
          A = VrtListA(IDetA,AA)
          B = VrtListA(IDetA,BB)
          JDetA = IndPQonA(A,I,IDetA)
          Sgn1 = Sign(1,JDetA)
          JDetA = JDetA*Sgn1
          KDetA = IndPQonA(B,J,JDetA)
          Sgn2 = Sign(1,KDetA)
          KDetA = KDetA*Sgn2
          HElement = Sgn1*Sgn2*H2aa(A,B,I,J)
          Do IDetB = 1,NDetB
            IDet = IndAB2Full(IDetA,IDetB)
            KDet = IndAb2Full(KDetA,IDetB)
            HPsi(KDet) = HPsi(KDet) + Psi(IDet)*HElement
          End Do
        End Do
        End Do
        End Do
        End Do
      End Do


! Step 4: Handle the bb double excitations
      Do IDetB = 1,NDetB
        Do II = 1,NOccB
        Do JJ = II+1,NOccB
        Do AA = 1,NVrtB
        Do BB = AA+1,NVrtB
          I = OccListB(IDetB,II)
          J = OccListB(IDetB,JJ)
          A = VrtListB(IDetB,AA)
          B = VrtListB(IDetB,BB)
          JDetB = IndPQonB(A,I,IDetB)
          Sgn1 = Sign(1,JDetB)
          JDetB = JDetB*Sgn1
          KDetB = IndPQonB(B,J,JDetB)
          Sgn2 = Sign(1,KDetB)
          KDetB = KDetB*Sgn2
          HElement = Sgn1*Sgn2*H2bb(A,B,I,J)
          Do IDetA = 1,NDetA
            IDet = IndAB2Full(IDetA,IDetB)
            KDet = IndAB2Full(IDetA,KDetB)
            HPsi(KDet) = HPsi(KDet) + Psi(IDet)*HElement
          End Do
        End Do
        End Do
        End Do
        End Do
      End Do

      Return
      End Subroutine BuildHC






      Subroutine ShutDownHam
      Implicit None
      Integer :: IAlloc
      DeAllocate(HDiag, FockSSa, FockSSb, FockOSa, FockOSb, FockA, FockB, HDoubA, HDoubB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ShutDownHam!"
      Return
      End Subroutine ShutDownHam

   End Module Hamiltonian

