
   Module ObaraSaika
   Use Precision
   Use Constants
   Implicit None
   Real (Kind=pr), Parameter :: Thresh = 1.0e-13_pr
   Contains

      Function Olap2c(ZA,CA,A,LA,ZB,CB,B,LB)
      Implicit None
      Real (Kind=pr), Intent(In)  :: ZA,    ZB       ! Exponents
      Real (Kind=pr), Intent(In)  :: CA,    CB       ! Norm Constants
      Real (Kind=pr), Intent(In)  :: A(3),   B(3)    ! Centers
      Integer,        Intent(In)  :: LA(3), LB(3)    ! Angular Momentum
      Real (Kind=pr)              :: Olap2C
! Local variables
      Real (Kind=pr), Allocatable :: OTmp(:,:)
      Real (Kind=pr) :: Zeta, DNrm, Xi, Olap
      Real (Kind=pr) :: P(3)
      Real (Kind=pr) :: Rab2
      Integer :: NA, NB
      Integer :: IA, IB, I

!================================================!
!  Two-center overlap integrals evaluated using  !
!  the Obara-Saika recursion scheme.             !
!                                                !
!  The functions we're integrating have:         !
!    ZA, ZB = Exponents                          !
!    CA, CB = Normalization constants            !
!    LA, LB = Angular momentum information       !
!    A,  B  = Centers                            !
!================================================!

      DNrm = CA * CB
      Zeta = ZA + ZB
      Xi   = ZA*ZB/Zeta
      P    = (ZA*A + ZB*B)/Zeta
      Rab2 = Dot_Product(A-B,A-B)
      NA   = MaxVal(LA)
      NB   = MaxVal(LB)


!===============================================!
!  Now, we initialize the <s|s> integral and    !
!  set up the recursion for the higher angular  !
!  momenta.  We're incrementing each Cartesian  !
!  component in turn.                           !
!===============================================!

      Olap = (Pi/Zeta)**F32 * Exp(-Xi*Rab2)
! TMH: New!
      If(Olap < Thresh) Then
        Olap2C = Zero
        Return
      End If
! TMH: Done!
      Allocate(OTmp(-1:NA,-1:NB))

      Do I = 1,3
       OTmp = Zero
       OTmp(0,0) = Olap

       Do IB = 0,LB(I)-1
        OTmp(0,IB+1) = (P(I) - B(I))*OTmp(0,IB)                &
                     + Real(IB)/(Two*Zeta)*OTmp(0,IB-1)
       End Do

       Do IB = 0,LB(I)
        Do IA = 0,LA(I)-1
         OTmp(IA+1,IB) = (P(I) - A(I))*OTmp(IA,IB)             &
                       + Real(IB)/(Two*Zeta)*OTmp(IA,IB-1)     &
                       + Real(IA)/(Two*Zeta)*OTmp(IA-1,IB)
        End Do
       End Do
             
       Olap = OTmp(LA(I),LB(I))
      End Do
      DeAllocate(OTmp)

      Olap2c = Olap*DNrm

      Return
      End Function Olap2c






      Function TInt(ZA,CA,A,LA,ZB,CB,B,LB)
      Implicit None
      Real (Kind=pr), Intent(In)  :: ZA,    ZB       ! Exponents
      Real (Kind=pr), Intent(In)  :: CA,    CB       ! Norm Constants
      Real (Kind=pr), Intent(In)  :: A(3),   B(3)    ! Centers
      Integer,        Intent(In)  :: LA(3), LB(3)    ! Angular Momentum
      Real (Kind=pr)              :: TInt
! Local variables
      Real (Kind=pr), Allocatable :: OTmp(:,:), TTmp(:,:)
      Real (Kind=pr) :: Zeta, DNrm, Xi, Olap
      Real (Kind=pr) :: P(3)
      Real (Kind=pr) :: Rab2
      Integer :: NA, NB
      Integer :: IA, IB, I

!================================================!
!  Kinetic energy integrals evaluated using the  !
!  Obara-Saika recursion scheme.                 !
!                                                !
!  The functions we're integrating have:         !
!    ZA, ZB = Exponents                          !
!    CA, CB = Normalization constants            !
!    LA, LB = Angular momentum information       !
!    A,  B  = Centers                            !
!================================================!

      DNrm = CA * CB
      Zeta = ZA + ZB
      Xi   = ZA*ZB/Zeta
      P    = (ZA*A + ZB*B)/Zeta
      Rab2 = Dot_Product(A-B,A-B)
      NA   = MaxVal(LA)
      NB   = MaxVal(LB)


!===============================================!
!  Now, we initialize the <s|s> integrals and   !
!  set up the recursion for the higher angular  !
!  momenta.  We're incrementing each Cartesian  !
!  component in turn.                           !
!===============================================!

      Olap = (Pi/Zeta)**F32 * Exp(-Xi*Rab2) 
      TInt = Xi*(Three - Two*Xi*Rab2)*Olap
! TMH: New!
      If(Olap < Thresh) Then
        TInt = Zero
        Return
      End If
! TMH: Done!
      Allocate(OTmp(-1:NA,-1:NB), TTmp(-1:NA,-1:NB))

      Do I = 1,3
       OTmp = Zero
       TTmp = Zero
       OTmp(0,0) = Olap
       TTmp(0,0) = TInt

       Do IB = 0,LB(I)-1
        OTmp(0,IB+1) = (P(I) - B(I))*OTmp(0,IB)                &
                     + Real(IB)/(Two*Zeta)*OTmp(0,IB-1)
        TTmp(0,IB+1) = (P(I) - B(I))*TTmp(0,IB)                &
                     + Two*Xi*OTmp(0,IB+1)                     &
                     + Real(IB)*(TTmp(0,IB-1)/(Two*Zeta)       &
                               - OTmp(0,IB-1)*Xi/ZB)
       End Do

       Do IB = 0,LB(I)
        Do IA = 0,LA(I)-1
         OTmp(IA+1,IB) = (P(I) - A(I))*OTmp(IA,IB)             &
                       + Real(IB)/(Two*Zeta)*OTmp(IA,IB-1)     &
                       + Real(IA)/(Two*Zeta)*OTmp(IA-1,IB)
         TTmp(IA+1,IB) = (P(I) - A(I))*TTmp(IA,IB)             &
                       + Two*Xi*OTmp(IA+1,IB)                  &
                       + Real(IA)*(TTmp(IA-1,IB)/(Two*Zeta)    &
                                 - OTmp(IA-1,IB)*Xi/ZA)        &
                       + Real(IB)*TTmp(IA,IB-1)/(Two*Zeta)
        End Do
       End Do
       Olap = OTmp(LA(I),LB(I))
       TInt = TTmp(LA(I),LB(I))
      End Do
      DeAllocate(OTmp,TTmp)

      TInt = TInt*DNrm
       
      Return
      End Function TInt






      Function VInt(ZA,CA,A,LA,ZB,CB,B,LB,Q,C)
      Implicit None
      Real (Kind=pr), Intent(In)  :: ZA,    ZB       ! Exponents
      Real (Kind=pr), Intent(In)  :: CA,    CB       ! Norm Constants
      Real (Kind=pr), Intent(In)  :: A(3),   B(3)    ! Centers
      Integer,        Intent(In)  :: LA(3), LB(3)    ! Angular Momentum
      Real (Kind=pr), Intent(In)  :: Q, C(3)         ! Charge info
      Real (Kind=pr)              :: VInt
! Local variables
      Real (Kind=pr), Allocatable :: VTmp(:,:,:), V(:)
      Real (Kind=pr) :: Zeta, DNrm, Xi, Olap
      Real (Kind=pr) :: P(3)
      Real (Kind=pr) :: Rab2, U
      Integer :: NA, NB, LTot
      Integer :: IA, IB, I, M

!================================================!
!  Coulomb integrals evaluated using the         !
!  Obara-Saika recursion scheme.                 !
!                                                !
!  The functions we're integrating have:         !
!    ZA, ZB = Exponents                          !
!    CA, CB = Normalization constants            !
!    LA, LB = Angular momentum information       !
!    A,  B  = Centers                            !
!  The Coulomb center has:                       !
!    Q      = Charge                             !
!    C      = Center                             !
!  If Q is positive, the integral is positive.   !
!================================================!

      DNrm = CA * CB * Q
      Zeta = ZA + ZB
      Xi   = ZA*ZB/Zeta
      P    = (ZA*A + ZB*B)/Zeta
      Rab2 = Dot_Product(A-B,A-B)
      LTot = Sum(LA+LB)
      NA   = MaxVal(LA)
      NB   = MaxVal(LB)


!=========================================================!
!  After initializing the <s|s> integral, we have a loop  !
!  over the nuclear centers, and we evaluate the nuclear  !
!  attraction integrals a nucleus at a time.  In each     !
!  integral, we increment each Cartesian component in     !
!  turn.  Since increasing the angular momentum requires  !
!  us to know auxiliary integrals indexed by m, we have   !
!  that loop to worry about as well.                      !
!=========================================================!

      Olap = (Pi/Zeta)**F32 * Exp(-Xi*Rab2)
      VInt = Zero

      Allocate(V(0:LTot), VTmp(-1:NA,-1:NB,0:LTot))

      U = Zeta*Dot_Product(P-C,P-C)
      Do M = 0, LTot
       V(M) = Two*Sqrt(Zeta/Pi)*Olap*FBoys(M,U)
      End Do
! TMH: New!
      If(MaxVal(Abs(V(0:LTot))) < Thresh) Then
        VInt = Zero
        Deallocate(V, VTmp)
        Return
      End If
! TMH: Done!

      Do I = 1,3
       VTmp = Zero
       VTmp(0,0,:) = V

       Do IB = 0,LB(I)-1
        Do M = 0,LTot-1
         VTmp(0,IB+1,M) = (P(I) - B(I))*VTmp(0,IB,M)                 &
                        - (P(I) - C(I))*VTmp(0,IB,M+1)               &
                        + Real(IB)/(Two*Zeta)*(VTmp(0,IB-1,M)        &
                                             - VTmp(0,IB-1,M+1))
        End Do
       End Do

       Do IB = 0,LB(I)
        Do IA = 0,LA(I)-1
         Do M = 0,LTot-1
          VTmp(IA+1,IB,M) = (P(I) - A(I))*VTmp(IA,IB,M)              &
                          - (P(I) - C(I))*VTmp(IA,IB,M+1)            &
                          + Real(IA)/(Two*Zeta)*(VTmp(IA-1,IB,M)     &
                                               - VTmp(IA-1,IB,M+1))  &
                          + Real(IB)/(Two*Zeta)*(VTmp(IA,IB-1,M)     &
                                               - VTmp(IA,IB-1,M+1))
         End Do
        End Do
       End Do

       V = VTmp(LA(I),LB(I),:)
      End Do

      VInt = V(0)*DNrm
      DeAllocate(V,VTmp)

      Return
      End Function VInt






      Function ERInt(ZA,CA,A,LA,ZB,CB,B,LB,ZC,CC,C,LC,ZD,CD,D,LD)
      Implicit None
      Real (Kind=pr), Intent(In)  :: ZA,   ZB,     ZC,    ZD     ! Exponents
      Real (Kind=pr), Intent(In)  :: CA,   CB,     CC,    CD     ! Norm Constants
      Real (Kind=pr), Intent(In)  :: A(3),   B(3),  C(3),  D(3)  ! Centers
      Integer,        Intent(In)  :: LA(3), LB(3), LC(3), LD(3)  ! Angular Momentum
      Real (Kind=pr)              :: ERInt
! Local variables
      Real (Kind=pr), Allocatable :: ETmp(:,:,:,:,:), E(:)
      Real (Kind=pr) :: Zeta, Eta, Rho, DNrm, Oab, Ocd
      Real (Kind=pr) :: P(3), Q(3), W(3)
      Real (Kind=pr) :: Rab2, Rcd2, Rpq2, T
      Integer :: NA, NB, NC, ND, LTot
      Integer :: IA, IB, IC, ID, M, I
   
!====================================================!
!  Electron repulsion integrals evaluated using the  !
!  Obara-Saika recursion scheme.                     !    
!                                                    !
!  The functions we're integrating have:             !
!    ZK = Exponent of function K                     !
!    CK = Normalization constant of function K       !
!    LK = Angular momentum for function K            !
!    K  = Vector pointing to center K                !
!  Note that this integral is in Mulliken notation,  !
!  that is, with functions A, B, C, and D we have    !
!    ERI(A,B,C,D) = (A B | 1/r12 | C D)              !
!                 = <A C | 1/r12 | B D>.             !
!====================================================!

      Zeta = ZA + ZB 
      Eta  = ZC + ZD
      Rho  = Eta*Zeta/(Eta+Zeta)
      DNrm = CA * CB * CC * CD
      P    = (ZA*A + ZB*B)/Zeta
      Q    = (ZC*C + ZD*D)/Eta
      W    = (Zeta*P + Eta*Q)/(Zeta+Eta)
      Rab2 = Dot_Product(A-B,A-B)
      Rcd2 = Dot_Product(C-D,C-D)
      Rpq2 = Dot_Product(P-Q,P-Q)
      T    = Rho*Rpq2
      LTot = Sum(LA+LB+LC+LD)
      NA   = MaxVal(LA)
      NB   = MaxVal(LB)
      NC   = MaxVal(LC)
      ND   = MaxVal(LD)


!=====================================================!
!  After initializing the <s|s> integrals, we start   !
!  incrementing the Cartesian components of each      !
!  function in turn.  Since increasing the angular    !
!  momentum requires auxiliary integrals indexed by   !
!  m, we have to loop over m as well.                 !
!=====================================================!

      Oab = (Pi/Zeta)**F32 * Exp(-ZA*ZB*Rab2/Zeta)
      Ocd = (Pi/ Eta)**F32 * Exp(-ZC*ZD*Rcd2/Eta)
      Allocate(E(0:LTot), ETmp(-1:NA,-1:NB,-1:NC,-1:ND,0:LTot))

      Do M = 0,LTot
       E(M) = Two*Sqrt(Rho/Pi)*Oab*Ocd*FBoys(M,T)
      End Do
! TMH: New!
      If(MaxVal(Abs(E(0:LTot))) < Thresh) Then
        ERInt = Zero
        Deallocate(E, ETmp)
        Return
      End If
! TMH: Done!

      Do I = 1,3
       ETmp = Zero
       ETmp(0,0,0,0,:) = E

       Do ID = 0,LD(I)-1
        Do M = 0,LTot-1
         ETmp(0,0,0,ID+1,M) = (Q(I) - D(I))*ETmp(0,0,0,ID,M)             &
                            + (W(I) - Q(I))*ETmp(0,0,0,ID,M+1)           &
               + Real(ID)/(Two*Eta)*(ETmp(0,0,0,ID-1,M)                  &
                                   - ETmp(0,0,0,ID-1,M+1)*Rho/Eta)
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)-1
         Do M = 0,LTot-1
          ETmp(0,0,IC+1,ID,M) = (Q(I) - C(I))*ETmp(0,0,IC,ID,M)          &
                              + (W(I) - Q(I))*ETmp(0,0,IC,ID,M+1)        &
               + Real(ID)/(Two*Eta)*(ETmp(0,0,IC,ID-1,M)                 &
                                   - ETmp(0,0,IC,ID-1,M+1)*Rho/Eta)      &
               + Real(IC)/(Two*Eta)*(ETmp(0,0,IC-1,ID,M)                 &
                                   - ETmp(0,0,IC-1,ID,M+1)*Rho/Eta)
         End Do
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)
         Do IB = 0,LB(I)-1
          Do M = 0,LTot-1
           ETmp(0,IB+1,IC,ID,M) = (P(I) - B(I))*ETmp(0,IB,IC,ID,M)       &
                                + (W(I) - P(I))*ETmp(0,IB,IC,ID,M+1)     &
               + Real(ID)/(Two*Eta+Two*Zeta)*ETmp(0,IB,IC,ID-1,M+1)      &
               + Real(IC)/(Two*Eta+Two*Zeta)*ETmp(0,IB,IC-1,ID,M+1)      &
               + Real(IB)/(Two*Zeta)*(ETmp(0,IB-1,IC,ID,M)               &
                                    - ETmp(0,IB-1,IC,ID,M+1)*Rho/Zeta)
          End Do
         End Do
        End Do
       End Do

       Do ID = 0,LD(I)
        Do IC = 0,LC(I)
         Do IB = 0,LB(I)
          Do IA = 0,LA(I)-1
           Do M = 0,LTot-1
            ETmp(IA+1,IB,IC,ID,M) = (P(I) - A(I))*ETmp(IA,IB,IC,ID,M)    &
                                  + (W(I) - P(I))*ETmp(IA,IB,IC,ID,M+1)  &
               + Real(ID)/(Two*Eta+Two*Zeta)*ETmp(IA,IB,IC,ID-1,M+1)     &
               + Real(IC)/(Two*Eta+Two*Zeta)*ETmp(IA,IB,IC-1,ID,M+1)     &
               + Real(IB)/(Two*Zeta)*(ETmp(IA,IB-1,IC,ID,M)              &
                                    - ETmp(IA,IB-1,IC,ID,M+1)*Rho/Zeta)  &
               + Real(IA)/(Two*Zeta)*(ETmp(IA-1,IB,IC,ID,M)              &
                                    - ETmp(IA-1,IB,IC,ID,M+1)*Rho/Zeta)
           End Do
          End Do
         End Do
        End Do
       End Do
       E = ETmp(LA(I),LB(I),LC(I),LD(I),:)
      End Do

      ERInt = E(0)*DNrm
      DeAllocate(E,ETmp)

      Return
      End Function ERInt






      Function FBoys(M,T)
      Implicit None
      Integer :: M
      Real (Kind=pr) :: T, FBoys
      Real (Kind=pr) :: GammP, A
      Real (Kind=pr), Parameter :: Thresh = 1.0e-4_pr
      Real (Kind=pr) :: R, X, S, Fac
      Integer :: K

      A = M + F12
      If(T < Zero) Stop "T must be non-negative in FBoys!"
      If(A < Zero) Stop "M must be non-negative in FBoys!"

! For very small T, there's a problem with GSer
      If(T < Thresh) Then
        R = Sqrt(T)
        X = Two*T
        FBoys = One/(2*M+1)
        Do K = 1,6
          FBoys = FBoys + (-T)**K/Fac(K) *One/(Two*M+Two*K+One)
        End Do
      Else If(T < A+1) Then
        Call GSer(GammP,A,T)    ! Compute from taylor series
        FBoys = GammP/(2*Sqrt(T)*T**M)
      Else
        Call GCF(GammP,A,T)     ! Compute from continued fraction
        FBoys = GammP/(2*Sqrt(T)*T**M)
      End If

      Return

   Contains

        Subroutine GSer(GammP,A,T)
        Real (Kind=pr), Intent(In)  :: A, T
        Real (Kind=pr), Intent(Out) :: GammP
        Real (Kind=pr), Parameter   :: Tol = 1.0e-12_pr
        Integer,        Parameter   :: MaxIter = 200
        Real (Kind=pr) :: AP, Term, Sum
        Integer :: N

        AP   = A
        Term = One/AP
        Sum  = Term
        Do N = 1,MaxIter
          AP = AP + 1
          Term = Term*T/AP
          Sum = Sum + Term
          If(Abs(Term) < Abs(Sum)*Tol) Exit
        End Do
        If(Abs(Term) > Abs(Sum)*Tol) Stop "GSer needs more terms!"
        GammP = Sum*Exp(-T+a*log(T))

        Return
        End Subroutine GSer




        Subroutine GCF(GammP,A,T)
        Real (Kind=pr), Intent(In)  :: A, T
        Real (Kind=pr), Intent(Out) :: GammP
        Integer,        Parameter   :: MaxIter = 200
        Real (Kind=pr), Parameter   :: Tol = 1.0e-12_pr, FPMin = 1.0e-30_pr
        Real (Kind=pr) :: B, C, D, H, AN, Term, GamA, DFac
        Integer :: M, I

        B = T + One-A
        C = One/FPMin
        D = One/B
        H = D
        Do I = 1,MaxIter
          AN = -I*(I-A)
          B  = B + 2
          D  = AN*D + B
          If(Abs(D) < FPMin) D = FPMin
          C  = B + AN/C
          If(Abs(C) < FPMin) C = FPMin
          D = One/D
          Term = D*C
          H = H*Term
          If(Abs(Term-One) < Tol) Exit
        End Do
        If(Abs(Term-One) > Tol) Stop "GCF did not converge"
        GammP = Exp(-T+A*Log(T))*H

! Now we just need to subtract this from Gamma(A) = Gamma(M+1/2)
        M = NInt(A-F12)
        GamA = DFac(M)/(2**M)*Sqrt(Pi)

        GammP = GamA-GammP
        Return
        End Subroutine GCF

      End Function FBoys

   End Module ObaraSaika

