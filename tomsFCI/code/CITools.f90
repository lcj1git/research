
   Module CITools
   Use Precision
   Use Constants
   Use IndexCI
   Implicit None
   Logical :: ExcLevelDetReady, ExcLevelDetAReady, ExcLevelDetBReady
   Integer, Allocatable :: ExcLevelDet(:), ExcLevelDetA(:), ExcLevelDetB(:)
   Contains

      Subroutine ExtractTAmplitudes(CIVec,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: CIVec(NDet)
      Real (Kind=pr), Allocatable :: T(:), ExpTPsi(:)
      Integer :: IAlloc, ExcLevel, ExcLevelMaxA, ExcLevelMaxB, ExcLevelMax

!====================================================================!
!  This guy should be pulling out the T amplitudes that give CIVec.  !
!  The idea is to read off C1 = T1, apply Exp(-T1) |Psi>, then read  !
!  C2 = T2, apply Exp(-T2) |Psi>, and so on and so forth.            !
!====================================================================!

! We'll need some memory
      Allocate(T(NDet), ExpTPsi(NDet), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ExtractTAmplitudes"


! Fill in ExcLevelDet*
      Call SetUpExcLevelDet(NOccA,NOccB,NAO,NDetA,NDetB,NDet)


! Start with ExpTPsi = Psi, in intemediate normalization
      ExpTPsi = CIVec/CIVec(1)


! Get the maximum excitation level.
      ExcLevelMaxA = MaxVal(ExcLevelDetA)
      ExcLevelMaxB = MaxVal(ExcLevelDetB)
      ExcLevelMax  = MaxVal(ExcLevelDet)


! Loop over excitation levels.
! At each level, first read Tn = Cn, then apply Exp(-Tn) |Psi>
      Do ExcLevel = 1,ExcLevelMax
        Call GetT(ExpTPsi,T,ExcLevel,ExcLevelDet,NDet)                          ! Read Tn = Cn
        Call GetExpTnPsi(ExpTPsi,-T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet)  ! |ExpTnPsi> => Exp(-Tn) |ExpTnPsi>
      End Do


! Now dump amplitudes to disk
      Call PrintT(T,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Call PrintT2abab(T,NOccA,NOccB,NAO,NDet,"TAmps")


! Deallocate and exit safely
      Call ShutDownExcLevel
      DeAllocate(T, ExpTPsi, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ExtractTAmplitudes"

      Return
      End Subroutine ExtractTAmplitudes






      Subroutine SetUpExcLevelDet(NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer, Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Integer :: IAlloc

!==========================================!
!  This allocates and fills ExcLevelDet*.  !
!==========================================!

! Allocate space
      Allocate(ExcLevelDet(NDet), ExcLevelDetA(NDetA), ExcLevelDetB(NDetB),   &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in SetUpExcLevelDet"

! Fill them in
      Call FillExcLevel(ExcLevelDet,ExcLevelDetA,ExcLevelDetB,NOccA,NOccB,NAO,NDetA,NDetB,NDet)

      Return
      End Subroutine SetUpExcLevelDet






      Subroutine FillExcLevel(ExcLevelDet,ExcLevelDetA,ExcLevelDetB,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer, Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Integer, Intent(Out) :: ExcLevelDet(NDet)
      Integer, Intent(Out) :: ExcLevelDetA(NDetA)
      Integer, Intent(Out) :: ExcLevelDetB(NDetB)
      Integer :: IDetA, IDetB, IDet, IExc

!==========================================!
!  Fills in excitation level information.  !
!==========================================!

! Set ExcLevelDetA
      ExcLevelDetA = 0
      Do IDetA = 1,NDetA
        IExc = Sum(IAmOccA(IDetA,NOccA+1:NAO))
        ExcLevelDetA(IDetA) = IExc
      End Do
      ExcLevelDetAReady = .true.


! Set ExcLevelDetB
      ExcLevelDetB = 0
      Do IDetB = 1,NDetB
        IExc = Sum(IAmOccB(IDetB,NOccB+1:NAO))
        ExcLevelDetB(IDetB) = IExc
      End Do
      ExcLevelDetBReady = .true.


! Set ExcLevelDet
      ExcLevelDet = 0
      Do IDetA = 1,NDetA
      Do IDetB = 1,NDetB
        IDet = IndAB2Full(IDetA,IDetB)
        IExc = ExcLevelDetA(IDetA) + ExcLevelDetB(IDetB)
        ExcLevelDet(IDet) = IExc
      End Do
      End Do
      ExcLevelDetReady = .true.

      Return
      End Subroutine FillExcLevel






      Subroutine ShutDownExcLevel
      Implicit None
      Integer :: IAlloc

!==========================================================!
!  This deallocate ExcLevelDet and declares it not ready.  !
!==========================================================!

      DeAllocate(ExcLevelDetA, ExcLevelDetB, ExcLevelDet, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ShutDownExcLevel"

      ExcLevelDetReady  = .false.
      ExcLevelDetAReady = .false.
      ExcLevelDetBReady = .false.

      Return
      End Subroutine ShutDownExcLevel






      Subroutine GetT(CIVec,T,ExcLevel,ExcLevelDet,NDet)
      Implicit None
      Integer,        Intent(In)    :: NDet, ExcLevel
      Integer,        Intent(In)    :: ExcLevelDet(NDet)
      Real (Kind=pr), Intent(In)    :: CIVec(NDet)
      Real (Kind=pr), Intent(InOut) :: T(NDet)
      Integer :: IDet,  IExc

!========================================================!
!  This subroutine copies C(ExcLevel) into T(ExcLevel).  !
!========================================================!

      Do IDet = 1,NDet
        IExc = ExcLevelDet(IDet)
        If(IExc /= ExcLevel) Cycle
        T(IDet) = CIVec(IDet)
      End Do

      Return
      End Subroutine GetT






      Subroutine GetExpTnPsi(Psi,T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Integer,        Intent(In)    :: ExcLevel
      Real (Kind=pr), Intent(In)    :: T(NDet)
      Real (Kind=pr), Intent(InOut) :: Psi(NDet)
      Real (Kind=pr), Allocatable   :: TPsi(:), TTPsi(:)
      Integer :: ITimes, IAlloc, NExcMax

!============================================================!
!  Given T, this just takes |Psi> and returns Exp(T) |Psi>.  !
!  We do only T(ExcLevel).  So                               !
!    GetExpTnPsi(Psi,T,2,NOcc,NSO,NDet) = Exp(T2) |Psi>      !
!    GetExpTnPsi(Psi,T,8,NOcc,NSO,NDet) = Exp(T8) |Psi>      !
!  and so on.                                                !
!============================================================!

! Allocate space
      Allocate(TPsi(NDet), TTPsi(NDet), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not Allocate in GetExpTnPsi"
      If(.not. ExcLevelDetReady) Stop "Need ExcLevelDet in GetExpTnPsi"


! Set up limits on how many times we apply T
      NExcMax = MaxVal(ExcLevelDet)/ExcLevel


! Now, at the top of the loop on count I we want TPsi = T^(I-1) Psi
      TPsi = Psi
      Do ITimes = 1,NExcMax
        Call GetTnPsi(TPsi,TTPsi,T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet)  ! TTPsi = T TPsi
        TPsi = TTPsi/ITimes                                                     ! Update TPsi; normalize for the 1/n!
        Psi = Psi + TPsi
      End Do


! DeAllocate and exit safely
      DeAllocate(TPsi, TTPsi, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not DeAllocate in GetExpTnPsi"

      Return
      End Subroutine GetExpTnPsi






      Subroutine GetTnPsi(Psi,TPsi,T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Integer,        Intent(In)  :: ExcLevel
      Real (Kind=pr), Intent(In)  :: T(NDet), Psi(NDet)
      Real (Kind=pr), Intent(Out) :: TPsi(NDet)
      Integer,        Allocatable :: dOccsA(:), KDetAOccs(:), ExcListA(:,:)
      Integer,        Allocatable :: dOccsB(:), KDetBOccs(:), ExcListB(:,:)
      Real (Kind=pr) :: TAmp
      Integer :: ExcLevelMaxA, ExcLevelMaxB, IAlloc
      Integer :: IDet, IDetA, IDetB, JDet, JDetA, JDetb, KDet, KDetA, KDetB
      Integer :: RelSignA, RelSignB, ExcLevelA, ExcLevelB
      Integer :: JDetAExcLevel, JDetBExcLevel

!============================================================!
!  This evalutes Tn |Psi> where n is specified by ExcLevel.  !
!============================================================!

! Allocate necessary space
      Allocate(dOccsA(NAO), KDetAOccs(NAO), ExcListA(ExcLevel,2),   &
               dOccsB(NAO), KDetBOccs(NAO), ExcListB(ExcLevel,2),   &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not Allocate in GetTnPsi"
      If(.not. ExcLevelDetAReady) Stop "Need ExcLevelDetA in GetTnPsi"
      If(.not. ExcLevelDetBReady) Stop "Need ExcLevelDetB in GetTnPsi"
      If(.not. ExcLevelDetReady)  Stop "Need ExcLevelDet  in GetTnPsi"


! This is the maximum possible excitation level
      ExcLevelMaxA = MaxVal(ExcLevelDetA)
      ExcLevelMaxB = MaxVal(ExcLevelDetB)


! Loop over excitations.  Retain those with the right excitation level
      TPsi = Zero
      Do IDet = 1,NDet
        If(ExcLevelDet(IDet) /= ExcLevel) Cycle


! Okay.  This determinant has the right excitation level
! Read off the alpha and beta determinants and their individual excitation levels
        IDetA = IndFull2AB(IDet,1)
        IDetB = IndFull2AB(IDet,2)
        ExcLevelA = ExcLevelDetA(IDetA)
        ExcLevelB = ExcLevelDetB(IDetB)


! Get the list of occupied annihilation operators and virtual creation operators for each spin
!   * If exciting out of orbital i, dOccs(i) = -1, else dOccs(i) = 0
!   * If exciting into orbital a, dOccs(a) = 1, else dOccs(a) = 0
        dOccsA = IAmOccA(IDetA,:)
        dOccsA(1:NOccA) = dOccsA(1:NOccA) - 1 
        dOccsB = IAmOccB(IDetB,:)
        dOccsB(1:NOccB) = dOccsB(1:NOccB) - 1 


! Convert into an excitation list and get the sign
! For example, suppose IDetA is |1 3 5 6>.  This is
!     |1 3 5 6> = - 5! 2 6! 4 |1 2 3 4>
! We'll read the sign and store that this is (5! 2) (6! 4)
        Call GetExcList(dOccsA,ExcListA,ExcLevelA,ExcLevel,NOccA,NAO,RelSignA,1)
        Call GetExcList(dOccsB,ExcListB,ExcLevelB,ExcLevel,NOccB,NAO,RelSignB,2)


! I have the alpha and beta signs.  Just save the value of the T amplitude
        TAmp = T(IDet)*RelSignA*RelSignB


! Now I know the value of the T amplitude itself and I know the operator string
! Loop over determinants and apply T(IDet) to each determinant
        Do JDetA = 1,NDetA
          JDetAExcLevel = ExcLevelDetA(JDetA)
          If(JDetAExcLevel + ExcLevelA > ExcLevelMaxA) Cycle


! Get the occupation numbers for KDetA = IDetAth excitation on JDetA
! Skip those not permitted
          KDetAOccs = IAmOccA(JDetA,:) + dOccsA
          If(Any(KDetAOccs > 1) .or. Any(KDetAOccs < 0)) Cycle


! Okay.  This excitation is allowed.  Calculate KDetA and its sign
! given that we're doing the excitation specified in ExcListA on JDetA
          Call Excite(JDetA,KDetA,RelSignA,ExcListA,ExcLevelA,ExcLevel,1)


! Do the same thing for B
          Do JDetB = 1,NDetB
            JDetBExcLevel = ExcLevelDetB(JDetB)
            If(JDetBExcLevel + ExcLevelB > ExcLevelMaxB) Cycle
            KDetBOccs = IAmOccB(JDetB,:) + dOccsB
            If(Any(KDetBOccs > 1) .or. Any(KDetBOccs < 0)) Cycle
            Call Excite(JDetB,KDetB,RelSignB,ExcListB,ExcLevelB,ExcLevel,2)


! Alright.  We turn |JDetA JDetB> into |KDetA KDetB> with an amplitude and a bunch of signs...
            JDet = IndAB2Full(JDetA,JDetB)
            KDet = IndAB2Full(KDetA,KDetB)
            TPsi(KDet) = TPsi(KDet) + TAmp*Psi(JDet)*RelSignA*RelSignB
          End Do
        End Do
      End Do


! Deallocate and done
      DeAllocate(dOccsA, KDetAOccs, ExcListA, dOccsb, KDetBOccs, ExcListB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in GetTnPsi"

      Return
      End Subroutine GetTnPsi






      Subroutine GetExcList(dOccs,ExcList,ExcLevelSpin,ExcLevel,NOcc,NAO,RelSign,ISpin)
      Implicit None
      Integer, Intent(In)  :: NOcc, NAO, ExcLevelSpin, ExcLevel, ISpin
      Integer, Intent(In)  :: dOccs(NAO)
      Integer, Intent(Out) :: ExcList(ExcLevel,2), RelSign
      Integer :: Ind, IDet, JDet, I, A

!===================================================================!
!  This takes how the occupation numbers change from the reference  !
!  and works out the list of creation and annihilation operators    !
!  and the sign needed to create the state with this change.        !
!                                                                   !
!  Do alpha or beta depending on ISpin.                             !
!===================================================================!

! Step 1: Loop over dOccs from 1 to NOcc.  Each time we hit -1, this is an occupied annihilation operator
      Ind = 0
      Do I = 1,NOcc
        If(dOccs(I) == -1) Then
          Ind = Ind + 1
          ExcList(Ind,1) = I
        End If
      End Do


! Step 2: Loop over dOccs from NOcc+1 to NAO.  Each time we get 1, this is a virtual creation operator
      Ind = 0
      Do A = NOcc+1,NAO
        If(dOccs(A) == 1) Then
          Ind = Ind + 1
          ExcList(Ind,2) = A
        End If
      End Do


! Step 3: get the sign
      RelSign = 1
      IDet = 1
      Do Ind = 1,ExcLevelSpin
        I = ExcList(Ind,1)
        A = ExcList(Ind,2)
        If(ISpin == 1) Then
          JDet = IndPQonA(A,I,IDet)   ! Look up sign of A! I |IDet>
        Else If(ISpin == 2) Then
          JDet = IndPQonB(A,I,IDet)   ! Look up sign of A! I |IDet>
        Else
          Stop "Incorrect ISpin in GetExcList"
        End If
        If(JDet < 0) RelSign = -RelSign
        IDet = Abs(JDet)
      End Do

      Return
      End Subroutine GetExcList






      Subroutine Excite(JDet,KDet,RelSign,ExcList,ExcLevelSpin,ExcLevel,ISpin)
      Integer, Intent(In)  :: ExcLevelSpin, ExcLevel, ISpin
      Integer, Intent(In)  :: JDet, ExcList(ExcLevel,2)
      Integer, Intent(Out) :: KDet, RelSign
      Integer :: Ind, I, A, Sgn, IDet

!=======================================================!
!  There's a list of excitations specified in ExcList.  !
!  We apply those excitations on JDet to produce KDet.  !
!  There's a possible sign in RelSign.                  !
!                                                       !
!  Do alpha or beta depending on ISpin.                 !
!=======================================================!

      IDet = JDet
      RelSign = 1
      Do Ind = 1,ExcLevelSpin
! Look up the indices for A! I
        I = ExcList(Ind,1)
        A = ExcList(Ind,2)

! Look up KDet = A! I IDet
        If(ISpin == 1) Then
          KDet = IndPQonA(A,I,IDet)
        Else If(ISpin == 2) Then
          KDet = IndPQonB(A,I,IDet)
        Else
          Stop "Incorrect ISpin in Excite!"
        End If

! Take care of sign and so on
        If(KDet < 0) RelSign = -RelSign
        IDet = Abs(KDet)
      End Do
      KDet = IDet

      Return
      End Subroutine Excite






      Subroutine GetExpTPsi(Psi,T,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)    :: T(NDet)
      Real (Kind=pr), Intent(InOut) :: Psi(NDet)
      Integer :: IExc

!====================================!
!  This build |Psi> => Exp(T) |Psi>  !
!====================================!

      If(.not. ExcLevelDetReady) Stop "Need ExcLevelDet in GetExpTPsi"
      Do IExc = 1,MaxVal(ExcLevelDet)
        Call GetExpTnPsi(Psi,T,IExc,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      End Do

      Return
      End Subroutine GetExpTPsi






      Subroutine GetTPsi(Psi,TPsi,T,NOcCA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: T(NDet), Psi(NDet)
      Real (Kind=pr), Intent(Out) :: TPsi(NDet)
      Real (Kind=pr), Allocatable :: TnPsi(:)
      Integer :: IExc, IAlloc

!===================================================!
!  This evalutes T |Psi> where it runs over all T.  !
!===================================================!

      Allocate(TnPsi(NDet), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in GetTPsi"
      If(.not. ExcLevelDetReady) Stop "Need ExcLevelDet in GetTPsi"

      TPsi = Zero
      Do IExc = 1,MaxVal(ExcLevelDet)
        Call GetTnPsi(Psi,TnPsi,T,IExc,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
        TPsi = TPsi + TnPsi
      End Do

      DeAllocate(TnPsi, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in GetTPsi"

      Return
      End Subroutine GetTPsi






      Subroutine PrintT(T,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: T(NDet)
      Real (Kind=pr), Allocatable :: NrmT(:)
      Integer :: ExcLevel, nMax

!=================================!
!  Print out the whole T vector.  !
!=================================!

      If(.not. ExcLevelDetReady) Stop "Need ExcLevelDet in PrintT"
      nMax = MaxVal(ExcLevelDet)

! Allocate space and print T for each excitation level
      Allocate(NrmT(nMax))
      Do ExcLevel = 1,nMax
        Call PrintTn(T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NrmT(ExcLevel))
      End Do
      Write(8,1000) NrmT
1000  Format(F15.10,15(3x,F15.10))

      Return
      End Subroutine PrintT






      Subroutine PrintTn(T,ExcLevel,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NrmT)
      Implicit None
      Integer,        Intent(In)  :: ExcLevel, NOccA, NOccB, NAO, NDetA, NDetB, NDet
      Real (Kind=pr), Intent(In)  :: T(NDet)
      Real (Kind=pr), Intent(Out) :: NrmT
!     Integer,        Allocatable :: dOccs(:), Create(:), Destroy(:)
!     Character (Len=120) :: ExcString, OccString, VrtString
!     Character (Len=4)   :: FName
!     Character (Len=3)   :: SpinOrbInd
!     Integer :: IDet, IDetExcLevel
!     Integer :: IOccA, IOccB, IVrtA, IVrtB
!     Integer :: IOcc, IAmNumber, RelSign
!     Integer :: IAlloc, NStates

!===========================================!
!  This prints the nth-order T amplitudes.  !
!  I need to check sign of T.               !
!===========================================!

!     Allocate(dOccs(NSO), Create(ExcLevel), Destroy(ExcLevel), Stat=IAlloc)
!     If(IAlloc /= 0) Stop "Could not Allocate in PrintTn"
!
!     If(ExcLevel < 10) Then
!       Write(FName,1001) ExcLevel
!     Else
!       Write(FName,1002) ExcLevel
!     End If
!1001  Format("T_",I1)
!1002  Format("T_",I2)
!
!      Open(7,File=Trim(AdjustL(FName)),Status="Replace")
!        NrmT = Zero
!        NStates = 0
!        Do IDet = 1,NDet
!! Are we writing this guy?
!          IDetExcLevel = Sum(IAmOcc(IDet,NOcc+1:NSO))
!          If(IDetExcLevel /= ExcLevel) Cycle
!          NStates = NStates + 1
!
!! I guess we are.  Now we need to make ExcString
!! Start by finding the alpha occupieds that are empty
!          OccString = ""
!          Do IOccA = 1,NOccA
!            IOcc = IOccA
!            If(IAmOcC(IDet,IOcc) == 0) Then
!              Write(SpinOrbInd,"(I3)") IOccA
!              OccString = Trim(AdjustL(OccString)) // " " // Trim(AdjustL(SpinOrbInd)) // "a "
!            End If
!          End Do
!
!! Now get the beta occupieds that are empty
!          Do IOccB = 1,NOccB
!            IOcc = IOccB + NOccA
!            If(IAmOcc(IDet,IOcc) == 0) Then
!              Write(SpinOrbInd,"(I3)") IOccB
!              OccString = Trim(AdjustL(OccString)) // " " // Trim(AdjustL(SpinOrbInd)) // "b "
!            End If
!          End Do
!
!! Now get the alpha virtuals that are not empty
!          VrtString = ""
!          Do IVrtA = NOccA+1,NAO
!            IOcc = IVrtA + NOccB
!            If(IAmOcc(IDet,IOcc) == 1) Then
!              Write(SpinOrbInd,"(I3)") IVrtA
!              VrtString = Trim(AdjustL(VrtString)) // " " // Trim(AdjustL(SpinOrbInd)) // "a "
!            End If
!          End Do
!
!! Now get the beta virtuals that are not empty
!          Do IVrtB = NOccB+1,NAO
!            IOcc = IVrtB + NAO
!            If(IAmOcc(IDet,IOcc) == 1) Then
!              Write(SpinOrbInd,"(I3)") IVrtB
!              VrtString = Trim(AdjustL(VrtString)) // " " // Trim(AdjustL(SpinOrbInd)) // "b "
!            End If
!          End Do
!          ExcString = "T(" // Trim(AdjustL(OccString)) // " => " // Trim(AdjustL(VrtString)) // ") = "
!
!
!! Quite separately, I need to get the sign of T.  That excitation string I've
!! written could take me to determinant IDet or -IDet.  Get it right.  Take this
!! straight from GetTnPsi
!          dOccs = IAmOcc(IDet,:)
!          dOccs(1:NOcc) = dOccs(1:NOcc) - 1
!          Destroy = 0
!          IAmNumber = 0
!          Do IOcc = 1,NOcc
!            If(dOccs(IOcc) == -1) Then
!              IAmNumber = IAmNumber + 1
!              Destroy(IAmNumber) = IOcc
!            End If
!          End Do
!          Create = 0
!          IAmNumber = 0
!          Do IOcc = NOcc+1,NSO
!            If(dOccs(IOcc) == 1) Then
!              IAmNumber = IAmNumber + 1
!              Create(IAmNumber) = IOcc
!            End If
!          End Do
!          Call GetRelSign(1,IDet,Create,Destroy,NOcc,ExcLevel,RelSign)
!
!          Write(7,2000) Trim(AdjustL(ExcString)), T(IDet)*RelSign
!          NrmT = NrmT + T(IDet)*T(IDet)
!        End Do
!        NrmT = Sqrt(NrmT/NStates)
!      Close(7)
!2000  Format(A,3x,ES20.12)
!
!      DeAllocate(dOccs, Create, Destroy, Stat=IAlloc)
!      If(IAlloc /= 0) Stop "Could not deallocate in PrintTn"
      Return
      End Subroutine PrintTn






      Subroutine PrintT2abab(T,NOccA,NOccB,NAO,NDet,FName)
      Implicit None
      Integer,           Intent(In) :: NOccA, NOccB, NAO, NDet
      Real (Kind=pr),    Intent(In) :: T(NDet)
      Character (Len=5), Intent(In) :: FName
      Integer :: I, J, A, B
      Integer :: IDetA, IDetB, IDet
      Integer :: SgnA, SgnB

!=======================!
!  Prints T2abab only.  !
!=======================!

      Open(7,File=FName,Status="Replace")
! Loop over single alpha excitations and single beta excitations
        Do I = 1,NOccA
        Do A = NOccA+1,NAO
        Do B = NOccB+1,NAO
        Do J = 1,NOccB

! Look up the alpha and beta excited determinants with sign
          IDetA = IndPQonA(A,I,1)
          IDetB = IndPQonB(B,J,1)
          SgnA = Sign(1,IDetA)
          SgnB = Sign(1,IDetB)


! Get the ab excited determinant index
          IDetA = SgnA*IDetA
          IDetB = SgnB*IDetB
          IDet = IndAB2Full(IDetA,IDetB)


! Write to disk
          Write(7,1000) I, J, A, B, T(IDet)*SgnA*SgnB
        End Do
        End Do
        End Do
        End Do
      Close(7)
1000  Format(4(I3,3x),ES20.12)

      Return
      End Subroutine PrintT2abab






      Subroutine WriteBasis(NOccA,NOccB,NDet)
      Implicit None
      Integer, Intent(In)  :: NOccA, NOccB, NDet
      Integer, Allocatable :: OccsA(:), OccsB(:)
      Integer :: IOcc, IDetA, IDetB, IDet, IAmOrb, IAlloc
      Character (Len=3)   :: ISpinOrb
      Character (Len=120) :: IDetString

!===========================================!
! Write determinant information to string.  !
!===========================================!

! We'll need some space
      Allocate(OccsA(NOccA), OccsB(NOccB), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in WriteBasis"


! Loop over determinants...
      Do IDet = 1,NDet
        IDetA = IndFull2AB(IDet,1)
        IDetB = IndFull2AB(IDet,2)
        OccsA = OccListA(IDetA,:)
        OccsB = OccListB(IDetB,:)
        IDetString = ""

! For each alpha occupied, get its index and update IDetString
        Do IOcc = 1,NOccA
          IAmOrb = OccsA(IOcc)
          Write(ISpinOrb,1000) IAmOrb, "a"
          IDetString = Trim(AdjustL(IDetString)) // " " // Trim(AdjustL(ISpinOrb))
        End Do

! Do the same for beta
        Do IOcc = 1,NOccB
          IAmOrb = OccsB(IOcc)
          Write(ISpinOrb,1000) IAmOrb, "b"
          IDetString = Trim(AdjustL(IDetString)) // " " // Trim(AdjustL(ISpinOrb))
        End Do
        Write(6,2000) IDet,Trim(AdjustL(IDetString))
      End Do
1000  Format(I2,A1)
2000  Format("Det #",I7," is  ",A)

      DeAllocate(OccsA, OccsB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in WriteBasis"

      Return
      End Subroutine WriteBasis

   End Module CITools

