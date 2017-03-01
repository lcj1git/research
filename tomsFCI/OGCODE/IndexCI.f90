
   Module IndexCI
   Implicit None
   Integer, Allocatable :: IndAB2Full(:,:)
   Integer, Allocatable :: IndFull2AB(:,:)
   Integer, Allocatable :: OccListA(:,:)
   Integer, Allocatable :: OccListB(:,:)
   Integer, Allocatable :: VrtListA(:,:)
   Integer, Allocatable :: VrtListB(:,:)
   Integer, Allocatable :: IAmOccA(:,:) 
   Integer, Allocatable :: IAmOccB(:,:)
   Integer, Allocatable :: IndDoubA(:,:)
   Integer, Allocatable :: IndDoubB(:,:)
   Integer, Allocatable :: IndPQonA(:,:,:)
   Integer, Allocatable :: IndPQonB(:,:,:)
   Integer :: NSingA, NSingB, NDoubA, NDoubB, NVrtA, NVrtB
   Contains

!==========================================================!
!  This module handles all sorts of indexing I'll need.    !
!                                                          !
!  I write all determinants as                             !
!     |IDet> = |IDetA IDetB>                               !
!  I thus need to be able to go back and forth between     !
!  IDet on the one hand and IDetA, IDetB on the other.     !
!  To do this, I use                                       !
!     IndAB2Full(IDetA,IDetB) = IDet                       !
!     IndFull2AB(IDet,:) = (/IDetA, IDetB /).              !
!                                                          !
!  To specify the alpha and beta determinants, I use the   !
!  arrays OccList* and IAmOcc*.  These are                 !
!     OccList(IDet,:) = list of orbitals occupied in IDet  !
!     IAmOcc(IDet,:) = list of occupation numbers in IDet  !
!  For example, given 6 AOs, the determinant |1 4 5> has   !
!     OccList = (1, 4, 5)                                  !
!     IAmOcc  = (1, 0, 0, 1, 1, 0)                         !
!                                                          !
!  These arrays are filled in SetUpIndexCI.                !
!----------------------------------------------------------!
!  In order to make the FCI more efficient, I also have    !
!  arrays which index excitations.  These are filled in    !
!  SetUpExcitationsIndices.                                !
!                                                          !
!  First, we have                                          !
!     IndDoubA(IDeTA,I) = JDetA                            !
!  which means that the Ith double excitation out of       !
!  IDetA is into JDetA.  Note that these don't say what    !
!  particle-hole excitation I corresponds to.  This will   !
!  probably be replaced.                                   !
!                                                          !
!  We remedy this flaw in IndPQonA and IndPQonB.  These    !
!  return                                                  !
!     IndPQonA(P,Q,IDetA) = +- JDetA                       !
!  where this means that                                   !
!     P! Q |IDetA> =  |JDetA> is sign is positive          !
!     P! Q |IDetA> = -|JDetA> if sign is negative          !
!     P! Q |IDetA> = 0        if value is zero.            !
!  So if I know which excitation operator I am using, I    !
!  can work out what state it excites me to.               !
!==========================================================!

      Subroutine SetUpIndexCI(NOccA,NOccB,NAO,NDetA,NDetB,NDet)
      Implicit None
      Integer, Intent(In)  :: NOccA, NOccB, NAO
      Integer, Intent(Out) :: NDetA, NDetB, NDet
      Integer :: IAlloc(8), Binom
      Integer :: IDetA, IDetB, IDet
      Integer :: IOcc, I, IVrt

!=====================================================!
!  This subroutine sets up the list of determinants.  !
!=====================================================!

! Allocate space
      NDetA = Binom(NAO,NOccA)
      NDetB = Binom(NAO,NOccB)
      NDet  = NDetA*NDetB
      NVrtA = NAO - NOccA
      NVrtB = NAO - NOccB
      IAlloc = 0
      Allocate(OccListA(NDetA,NOccA),   Stat=IAlloc(1))
      Allocate(OccListB(NDetB,NOccB),   Stat=IAlloc(2))
      Allocate(VrtListA(NDetA,NVrtA),   Stat=IAlloc(3))
      Allocate(VrtListB(NDetB,NVrtB),   Stat=IAlloc(4))
      Allocate(IAmOccA(NDetA,NAO),      Stat=IAlloc(5))
      Allocate(IAmOccB(NDetB,NAO),      Stat=IAlloc(6))
      Allocate(IndAB2Full(NDetA,NDetB), Stat=IAlloc(7))
      Allocate(IndFull2AB(NDet,2),      Stat=IAlloc(8))
      If(Any(IAlloc /= 0)) Stop "Could not allocate in SetUpIndexCI"


! Build OccListA, VrtListA, and IAmOccA
      OccListA = 0
      IAmOccA  = 0
      Do IDetA = 1,NDetA
        Call GetDet(NAO,NOccA,IDetA,OccListA(IDetA,:))
        Do IOcc = 1,NOccA
          I = OccListA(IDetA,IOcc)
          IAmOccA(IDetA,I) = 1
        End Do
        IVrt = 1
        Do I = 1,NAO
          If(IAmOccA(IDetA,I) == 1) Cycle
          VrtListA(IDetA,IVrt) = I
          IVrt = IVrt + 1
        End Do
      End Do


! Build OccListB, VrtListB, and IAmOccB
      OccListB = 0
      Do IDetB = 1,NDetB
        Call GetDet(NAO,NOccB,IDetB,OccListB(IDetB,:))
        Do IOcc = 1,NOccB
          I = OccListB(IDetB,IOcc)
          IAmOccB(IDetB,I) = 1
        End Do
        IVrt = 1
        Do I = 1,NAO
          If(IAmOCcB(IDetB,I) == 1) Cycle
          VrtListB(IdetB,IVrt) = I
          IVrt = IVrt + 1
        End Do
      End Do


! Now we've built the alpha and beta determinants
! Build IndAB2Full and IndFull2AB
      IndAB2Full = 0
      IndFull2AB = 0
      IDet = 0
      Do IDetA = 1,NDetA
      Do IDetB = 1,NDetB
        IDet = IDet + 1
        IndAB2Full(IDetA,IDetB) = IDet
        IndFull2AB(IDet,:) = (/IDetA, IDetB/)
      End Do
      End Do


! Set up the excitation indexing
      Call SetUpExcitationIndices(NOccA,NOccB,NAO,NDetA,NDetB)

      Return
      End Subroutine SetUpIndexCI






      Subroutine SetUpExcitationIndices(NOccA,NOccB,NAO,NDetA,NDetB)
      Implicit None
      Integer, Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB
      Integer, Allocatable :: OccsA(:), OccsB(:), VrtsA(:), VrtsB(:), OccsA2(:), OccsB2(:)
      Integer :: IAlloc, IDet, JDet, KDet, ISing, IDoub, Sgn
      Integer :: I, J, A, B, II, JJ, AA, BB

!==========================================================!
!  This sets up all the indexing for excitations.          !
!                                                          !
!  Recall                                                  !
!     IndPQonA(P,Q,IDetA) = +- JDetA                       !
!  where this means that                                   !
!     P! Q |IDetA> =  |JDetA> is sign is positive          !
!     P! Q |IDetA> = -|JDetA> if sign is negative          !
!     P! Q |IDetA> = 0        if value is zero.            !
!  Generically, to find IndPQ we do the following: for a   !
!  given valid P/Q/IDet combination, we                    !
!    * Occupy P the position where Q was                   !
!    * Sort the new list to get the sign                   !
!    * Look up the determinant in the new list             !
!                                                          !
!  For doubles, we need only unique doubles.  Algorithm:   !
!    * For a given determinant IDet,i loop over occupied   !
!      pairs I > J and virtual pairs A > B.                !
!    * Get |JDet> = A! I |IDet>                            !
!    * Get |KDet> = B! J |JDet>                            !
!    * Then this double excitation from IDet is KDet.      !
!  We're not fussed with sign here, although we could be.  !
!==========================================================!

! Compute the number of excitations of each spin
      NSingA = NOccA*NVrtA
      NSingB = NOccB*NVrtB
      NDoubA = NSingA*(NOccA-1)*(NVrtA-1)/4
      NDoubB = NSingB*(NOccB-1)*(NVrtB-1)/4


! Allocate global and local variables
      Allocate(IndDoubA(NDetA,NDoubA),       &
               IndDoubB(NDetB,NDoubB),       &
               IndPQonA(NAO,NAO,NDetA),      &
               IndPQonB(NAO,NAO,NDetB),      &
               OccsA(NOccA), OccsA2(NOccA),  &
               OccsB(NOccB), OccsB2(NOccB),  &
               VrtsA(NVrtA), VrtsB(NVrtB),   &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not Allocate in SetUpExcitationIndices"


! Build IndPQonA
      IndPQonA = 0
      Do IDet = 1,NDetA
        OccsA = OccListA(IDet,:)
        VrtsA = VrtListA(IDet,:)
        Do II = 1,NOccA
          I = OccsA(II)
          IndPQonA(I,I,IDet) = IDet
          Do AA = 1,NVrtA
            A = VrtsA(AA)
            OccsA2 = OccsA
            OccsA2(II) = A
            Call Sort(OccsA2,NOccA,Sgn)
            Call GetDetInv(NAO,NOccA,JDet,OccsA2)
            IndPQonA(A,I,IDet) = Sgn*JDet
          End Do
        End Do
      End Do


! Build IndPQonB
      IndPQonB = 0
      Do IDet = 1,NDetB
        OccsB  = OccListB(IDet,:)
        VrtsB  = VrtListB(IDet,:)
        Do II = 1,NOccB
          I = OccsB(II)
          IndPQonB(I,I,IDet) = IDet
          Do AA = 1,NVrtB
            A = VrtsB(AA)
            OccsB2 = OccsB
            OccsB2(II) = A
            Call Sort(OccsB2,NOccB,Sgn)
            Call GetDetInv(NAO,NOccB,JDet,OccsB2)
            IndPQonB(A,I,IDet) = Sgn*JDet
          End Do
        End Do
      End Do


! Double excitations are more annoying, but basically follow the same logic
      IndDoubA = 0
      Do IDet = 1,NDetA
        OccsA = OccListA(IDet,:)
        VrtsA = VrtListA(IDet,:)
        IDoub = 1
! Loop over occupied pairs PQ and virtual pairs RS
        Do II = 1,NOccA
        Do JJ = II+1,NOccA
        Do AA = 1,NVrtA
        Do BB = AA+1,NVrtA
          I = OccsA(II)
          J = OccsA(JJ)
          A = VrtsA(AA)
          B = VrtsA(BB)
          JDet = Abs(IndPQonA(A,I,IDet))
          KDet = Abs(IndPQonA(B,J,JDet))
          IndDoubA(IDet,IDoub) = KDet
          IDoub = IDoub + 1
        End Do 
        End Do 
        End Do 
        End Do 
      End Do 


! And IndDoubB
      IndDoubB = 0
      Do IDet = 1,NDetB
        OccsB = OccListB(IDet,:)
        VrtsB = VrtListB(IDet,:)
        IDoub = 1
! Loop over occupied pairs PQ and virtual pairs RS
        Do II = 1,NOccB
        Do JJ = II+1,NOccB
        Do AA = 1,NVrtB
        Do BB = AA+1,NVrtB
          I = OccsB(II)
          J = OccsB(JJ)
          A = VrtsB(AA)
          B = VrtsB(BB)
          JDet = Abs(IndPQonB(A,I,IDet))
          KDet = Abs(IndPQonB(B,J,JDet))
          IndDoubB(IDet,IDoub) = KDet
          IDoub = IDoub + 1
        End Do 
        End Do 
        End Do 
        End Do 
      End Do 


! Deallocate and done
      DeAllocate(OccsA, OccsB, OccsA2, OccsB2, VrtsA, VrtsB, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in SetUpExcitationIndices"

      Return
      End Subroutine SetUpExcitationIndices






      Subroutine GetDet(NAO,NOcc,IDet,OccList)
      Implicit None
      Integer, Intent(In)  :: NAO, NOcc, IDet
      Integer, Intent(Out) :: OccList(NOcc)
! Local variables
      Integer K, I, R, Binom

!==============================================================!
!  This subroutine creates an occupation number list OccList   !
!  for NOcc electrons in NBF orbitals.  IComb specifies which  ! 
!  combination of the Binom(NAO,NOcc) we're after.             !
!                                                              !
!  ACM Algorithm 515, B. P. Buckles and M. Lybanon             !
!                                                              !
!  Restrictions: 1 <= IComb <= Binom(NAO,NOcc)                 !
!                NOcc <= NAO                                   !
!==============================================================!

! Handle the case NOcc = 0
      OccList = 0
      If(NOcc == 0) Return


! Handle the case NOcc = 1
      If(NOcc == 1) Then
        OccList(1) = IDet
        Return
      End If


! Begin executable by initializing lower bound index
      K = 0
      Do I = 1,NOcc-1
        If(I /= 1) OccList(I) = OccList(I-1)
! Loop to check validity of each element
        Do
          OccList(I) = OccList(I) + 1
          R = Binom(NAO-OccList(I),NOcc-I)
          K = K+R
          If(K >= IDet) Exit
        End Do
        K = K-R
      End Do
      OccList(NOcc) = OccList(NOcc-1) + IDet-K

      Return
      End Subroutine GetDet






      Subroutine GetDetInv(NBF,NOcc,IComb,OccList)
      Implicit None
      Integer, Intent(In)  :: NBF, NOcc, OccList(NOcc)
      Integer, Intent(Out) :: IComb
! Local variables
      Integer :: K, I, R, Binom
      Integer :: OccListTemp(NOcc)

!==========================================================!
!  This is the inverse of Buckles and Lybanon, taken from  !
!       http://saliu.com/bbs/messages/348.html             !
!==========================================================!

! Handle the case NOcc = 0
      If(NOcc == 0) Return


! Handle the case NOcc = 1
      If(NOcc == 1) Then
        IComb = OccList(1)
        Return
      End If


! Begin executable by initializing lower bound index
      IComb = 0
      K = 0
      OccListTemp = 0
      Do I = 1,NOcc-1
        If(I /= 1) OccListTemp(I) = OccListTemp(I-1)
        Do
          OccListTemp(I) = OccListTemp(I) + 1
          R = Binom(NBF-OccListTemp(I),NOcc-I)
          K = K+R
          If(OccListTemp(I) >= OccList(I)) Exit
        End Do
        K = K-R
      End Do
      IComb = K + OccList(NOcc) - OccList(NOcc-1)

      Return
      End Subroutine GetDetInv






      Subroutine Sort(List,LenList,Sgn)
      Implicit None
      Integer, Intent(In)    :: LenList
      Integer, Intent(InOut) :: List(LenList)
      Integer, Intent(Out)   :: Sgn
      Integer :: I, J, Tmp

!==================================================================!
!  Sorts a list in place and returns the sign of the permutation.  !
!  We sort in ascending order.                                     !
!==================================================================!

      Sgn = 1
      Do I = 2,LenList
      Do J = 1,I-1
        If(List(I) < List(J)) Then
          Sgn = -Sgn
          Tmp = List(J)
          List(J) = List(I)
          List(I) = Tmp
        End If
      End Do
      End Do

      Return
      End Subroutine Sort






      Subroutine ShutDownIndexCI
      Implicit None
      Integer :: IAlloc(12)
      IAlloc = 0
      DeAllocate(IndAB2Full, Stat=IAlloc(1))
      DeAllocate(IndFull2AB, Stat=IAlloc(2))
      DeAllocate(OccListA,   Stat=IAlloc(3))
      DeAllocate(OccListB,   Stat=IAlloc(4))
      DeAllocate(VrtListA,   Stat=IAlloc(5))
      DeAllocate(VrtListB,   Stat=IAlloc(6))
      DeAllocate(IAmOccA,    Stat=IAlloc(7))
      DeAllocate(IAmOccB,    Stat=IAlloc(8))
      DeAllocate(IndDoubA,   Stat=IAlloc(9))
      DeAllocate(IndDoubB,   Stat=IAlloc(10))
      DeAllocate(IndPQonA,   Stat=IAlloc(11))
      DeAllocate(IndPQonB,   Stat=IAlloc(12))
      If(Any(IAlloc /= 0)) Stop "Could not DeAllocate in ShutDownIndexCI"
      Return
      End Subroutine ShutDownIndexCI
     
   End Module IndexCI

