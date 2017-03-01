
   Module DIIS
   Use Precision
   Use Constants
   Implicit None
   Integer,        Parameter   :: NDIIS = 5
   Integer,        Parameter   :: LenV  = NDIIS+1
   Integer,        Parameter   :: StartDIIS = NDIIS + 5
   Real (Kind=pr), Allocatable :: ErrVecs(:,:)
   Real (Kind=pr) :: M(LenV,LenV)
   Real (Kind=pr) :: V(LenV)
   Contains

      Subroutine SetUpDIIS(LenVec)
      Implicit None
      Integer, Intent(In) :: LenVec
      Allocate(ErrVecs(LenVec,NDIIS))
      Return
      End Subroutine SetUpDIIS






      Subroutine DoDIIS(Coeffs)
      Implicit None
      Real (Kind=pr), Intent(Out) :: Coeffs(NDIIS)
      Integer :: Info, IPiv(LenV)
      Integer :: I, J

!========================================================!
!  Given error vectors ErrVecs, set up the DIIS problem  !
!  and get the coefficient array Coeffs.                 !
!========================================================!

! Initialize V and M
      V = Zero
      M = One
      V(LenV) = One
      M(LenV,LenV) = Zero


! Build the B-matrix part of M
      Do I = 1,NDIIS
      Do J = 1,NDIIS
          M(I,J) = Dot_Product(ErrVecs(:,I),ErrVecs(:,J))
      End Do
      End Do


! Solve M.C = V, placing the solution in V
      Call DGeSV(LenV,1,M,LenV,IPiv,V,LenV,Info)
      If(Info == 0) Coeffs = V(1:NDIIS)

      Return
      End Subroutine DoDIIS






      Subroutine ShutDownDIIS
      Implicit None
      DeAllocate(ErrVecs)
      Return
      End Subroutine ShutDownDIIS

   End Module DIIS

