
   Module IntTrans
   Use Precision
   Use Constants
   Contains

      Subroutine IntTran2(Matrix,Evec,NBF)
      Implicit None
      Integer,        Intent(In)    :: NBF
      Real (Kind=pr), Intent(In)    :: Evec(NBF,NBF)
      Real (Kind=pr), Intent(InOut) :: Matrix(NBF,NBF)
      Real (Kind=pr), Allocatable   :: Scr(:,:)
      Integer :: IAlloc

!=============================================!
!  Does the 2-index integral transformation.  !
!                                             !
!  There ia an annoying subtlety:             !
!  In the CC code, we store                   !
!    T_i^a = T(I,A)                           !
!  which mean lower indices come first.  But  !
!  lower indices are associated with the ket  !
!  and my SCF code has put ket indices last.  !
!  Thus, we need a transpose.                 !
!=============================================!

      Allocate(Scr(NBF,NBF), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran2"

      Scr = MatMul(Matrix,Evec)
      Matrix = MatMul(Transpose(Evec),Scr)
      Matrix = Transpose(Matrix)

      DeAllocate(Scr, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran2"

      Return
      End Subroutine IntTran2






      Subroutine IntTran4(ERI,Evec1,Evec2,NBF,AntiSymm)
      Implicit None
      Logical,        Intent(In)    :: AntiSymm
      Integer,        Intent(In)    :: NBF
      Real (Kind=pr), Intent(In)    :: Evec1(NBF,NBF)
      Real (Kind=pr), Intent(In)    :: Evec2(NBF,NBF)
      Real (Kind=pr), Intent(InOut) :: ERI(NBF,NBF,NBF,NBF)
      Real (Kind=pr), Allocatable   :: ERI2(:,:,:,:)
      Integer :: IAlloc
      Integer :: Mu, Nu, Lam, Sig, P, Q, R, S

!=============================================!
!  Does the 4-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!                                             !
!  Also important note: see transpose above.  !
!=============================================!

      Allocate(ERI2(NBF,NBF,NBF,NBF), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran4"

! Transform index 4 - in Dirac notation, we're builidng <Mu Lam | Nu P>
! This uses Evec2 because it's electron 2
      Do Mu = 1,NBF
      Do Nu = 1,NBF
      Do Lam = 1,NBF
      Do P = 1,NBF
        ERI2(Mu,Nu,Lam,P) = Dot_Product(Evec2(:,P),ERI(Mu,Nu,Lam,:))
      End Do
      End Do
      End Do
      End Do


! Transform index 3 - in Dirac notation, we're building <Mu Q | Nu P>
! This uses Evec2 because it's electron 2 again
      Do Mu = 1,NBF
      Do Nu = 1,NBF
      Do P  = 1,NBF
      Do Q = 1,NBF
        ERI(Mu,Nu,Q,P) = Dot_Product(Evec2(:,Q),ERI2(Mu,Nu,:,P))
      End Do
      End Do
      End Do
      End Do


! Transform index 2 - in Dirac notation, we're building <Mu Q | R P>
! This time we use Evec1
      Do Mu = 1,NBF
      Do P = 1,NBF
      Do Q = 1,NBF
      Do R = 1,NBF
        ERI2(Mu,R,Q,P) = Dot_Product(Evec1(:,R),ERI(Mu,:,Q,P))
      End Do
      End Do
      End Do
      End Do


! Transform index 1 and drop to Dirac notation
! Once again, we use Evec1
      Do P = 1,NBF
      Do Q = 1,NBF
      Do R = 1,NBF
      Do S = 1,NBF
        ERI(S,Q,R,P) = Dot_Product(Evec1(:,S),ERI2(:,R,Q,P))
      End Do
      End Do
      End Do
      End Do


! Here we take our transpose and possibly antisymmetrize
      If(AntiSymm) Then
        Do P = 1,NBF
        Do Q = 1,NBF
        Do R = 1,NBF
        Do S = 1,NBF
          ERI2(P,Q,R,S) = ERI(R,S,P,Q) - ERI(S,R,P,Q)
        End Do
        End Do
        End Do
        End Do
      Else
        Do P = 1,NBF
        Do Q = 1,NBF
        Do R = 1,NBF
        Do S = 1,NBF
          ERI2(P,Q,R,S) = ERI(R,S,P,Q)
        End Do
        End Do
        End Do
        End Do
      End If
      ERI = ERI2


! Deallocate and exit safely
      DeAllocate(ERI2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran4"

      Return
      End Subroutine IntTran4

   End Module IntTrans

