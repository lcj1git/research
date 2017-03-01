
   Module CI
   Use Precision
   Use Constants
   Use Hamiltonian
   Implicit None
   Contains

      Subroutine SolveCI(H2aa,H2bb,H2ab,CIVecs,ESCF,ECI,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NRoots)
      Implicit None
      Integer,        Intent(In)  :: NOccA, NOccB, NAO, NDetA, NDetB, NDet, NRoots
      Real (Kind=pr), Intent(In)  :: ESCF
      Real (Kind=pr), Intent(In)  :: H2aa(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2bb(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)  :: H2ab(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(Out) :: CIVecs(NDet,NRoots)
      Real (Kind=pr), Intent(Out) :: ECI(NRoots)
! Local variables
      Real (Kind=pr) :: Tmp, ECISorted(NRoots)
      Logical :: Converged, Fail
      Integer :: iSweep, I, J, iRoot, NFound
      Integer, Parameter :: NSweeps = 50

!====================================================================!
!  Solve HC = EC by Lanczos diagonalization.                         !
!  We exit from the subroutine Lanczos on one of three conditions:   !
!    1) We have converged HC = EC to tolerance                       !
!    2) The Krylov space is linearly dependent                       !
!    3) The metric has a negative eigenvalue                         !
!                                                                    !
!  If condition one happens, we're done.                             !
!                                                                    !
!  If condition two happens, we take our best guess for CI and use   !
!  it as an initial guess for another Lanczos diagonalization        !
!                                                                    !
!  If condition three happens, the Krylov space can not be improved  !
!  further and we exit with an error.                                !
!====================================================================!

! Initialize the CI vectors.
      Call InitializeCIVecs(CIVecs,NDet,NRoots)


! Write the header to the output file
      Open(7,File="Output",Position="Append")
        Write(7,1000)
      Close(7)



! Do the Lanczos sweeps
      Do NFound = 0,NRoots-1
        Do iSweep = 1,NSweeps
! Write the header for this sweep/root combo
          Open(7,File="Output",Position="Append")
            Write(7,2000) iSweep, NFound +1
            Write(7,3000)
          Close(7)

! Do the Lanczos sweep and write the result
          Call Lanczos(CIVecs,ECI,H2aa,H2bb,H2ab,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NFound,NFound+1,Converged,Fail)
          Open(7,File="Output",Position="Append")
            Write(7,3000)
          Close(7)

! If the root converged, no more sweeps are needed
! If the Lanczos failed utterly, no more sweeps are possible
          If(Converged) Exit
          If(Fail) Stop "Lanczos diagonalization failed."
        End Do


! If after NSweeps sweeps we didn't converge, stop here
        If(.not. Converged) Exit


! Okay.  The root is converged.  Write to output and go to the next root
        Open(7,File="Output",Position="Append")
          Write(7,4000) NFound +1, ECI(NFound+1)
          Write(7,5000)
        Close(7)
      End Do


! Final output -- sort the roots we did find
      ECISorted = ECI
      Do I = 1,NFound
      Do J = 1,I-1
        If(ECISorted(J) > ECISorted(I)) Then
          Tmp = ECISorted(J)
          ECISorted(J) = ECISorted(I)
          ECISorted(I) = Tmp
        End If
      End Do
      End Do


! Write to disk
      Open(7,File="Output",Position="Append")
        Do iRoot = 1,NRoots
          Write(7,6000) iRoot, ECI(iRoot)
        End Do
        Write(7,6001)
      Close(7)

1000  Format(/,14x,50('*'),/,                               &
             14x,'*',15x,'FCI summary follows',14x,'*',/,   &
             14x,'*',48('='),'*')
2000  Format(14x,'*',8x,'Lanczos Sweep #',I3," for Root #",I3,8x,'*')
3000  Format(14x,'*',48('-'),'*')
4000  Format(14x,'*  Root #',I2,' has converged to ',F16.10,'    *')
5000  Format(14x,'*',48('='),'*')
6000  Format(14x,'*  FCI Energy of state ',I2,' is ',F18.12,'  *')
6001  Format(14x,50('*'))

      Return
      End Subroutine SolveCI






      Subroutine InitializeCIVecs(CIVecs,NDet,NRoots)
      Implicit None
      Integer,        Intent(In)  :: NDet, NRoots
      Real (Kind=pr), Intent(Out) :: CIVecs(NDet,NRoots)
      Real (Kind=pr), Allocatable :: SortedHDiag(:)
      Integer,        Allocatable :: StateOrder(:)
      Real (Kind=pr) :: FTmp
      Integer        :: ITmp, I, J, IAlloc

!=====================================================================!
!  This initializes CIVecs to the NRoots lowest energy determinants.  !
!=====================================================================!

! Allocate
      Allocate(SortedHDiag(NDet), StateOrder(NDet), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in InitializeCIVecs"


! Initialize SortedHDiag to HDiag, set StateOrder
      SortedHDiag = HDiag
      Do I = 1,NDet
        StateOrder(I) = I
      End Do


! Sort states from lowest to highest
      Do I = 1,NDet
      Do J = 1,I-1
        If(SortedHDiag(J) > SortedHDiag(I)) Then
          FTmp = SortedHDiag(J)
          SortedHDiag(J) = SortedHDiag(I)
          SortedHDiag(I) = FTmp
          ITmp = StateOrder(J)
          StateOrder(J) = StateOrder(I)
          StateOrder(I) = ITmp
        End If
      End Do
      End Do


! Write my guess
      CIVecs = Zero
      Do I = 1,NRoots
        CIVecs(I,I) = One
      End Do


! Deallocate and exit safely
      DeAllocate(SortedHDiag, StateOrder, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in InitializeCIVecs"

      Return
      End Subroutine InitializeCIVecs






      Subroutine Lanczos(CIVecs,ECI,H2aa,H2bb,H2ab,NOccA,NOccB,NAO,NDetA,NDetB,NDet,NFound,NRoots,Converged,Fail)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NOccB, NAO, NDetA, NDetB, NDet, NFound, NRoots
      Real (Kind=pr), Intent(In)    :: H2aa(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)    :: H2bb(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(In)    :: H2ab(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(InOut) :: CIVecs(NDet,NRoots)
      Real (Kind=pr), Intent(Out)   :: ECI(NRoots)
      Logical,        Intent(Out)   :: Converged, Fail
! Krylov space variables
      Real (Kind=pr), Allocatable :: q(:,:),  Hq(:,:)
      Real (Kind=pr), Allocatable :: Matrix(:,:), Metric(:,:)
      Real (Kind=pr), Allocatable :: Evecs(:,:),  Evals(:), HC(:,:)
      Real (Kind=pr), Allocatable :: X(:), Y(:), W(:)
      Real (Kind=pr), Allocatable :: Res(:)
! Miscellaneous local variables
      Real (Kind=pr) :: ResNrm, Alp, E, S
      Integer        :: I, J, IFlag, iRoot
      Integer        :: IAlloc(11)
      Logical        :: OpenClose
! Miscellaneous parameters
      Integer,        Parameter :: MaxDim = 50     ! Max size of Krylov space
      Real (Kind=pr), Parameter :: Tol = 1.0e-8_pr ! Tolerance on HC-EC

!======================================================!
!  Does Lanczos diagonalization                        !
!                                                      !
!  Usage: Call this routine over and over again until  !
!  Converged = .true.  If we exit with Fail = .true.,  !
!  we have negative eigenvalues of the metric of the   !
!  Krylov subspace and can go no further.              !
!======================================================!

! Step 1: Allocate space
      Allocate(q(NDet,MaxDim),          Stat=IAlloc(1))
      Allocate(Hq(NDet,MaxDim),         Stat=IAlloc(2))
      Allocate(Matrix(MaxDim,MaxDim),   Stat=IAlloc(3))
      Allocate(Metric(MaxDim,MaxDim),   Stat=IAlloc(4))
      Allocate(Evecs(MaxDim,MaxDim),    Stat=IAlloc(5))
      Allocate(Evals(MaxDim),           Stat=IAlloc(6))
      Allocate(HC(NDet,NRoots),         Stat=IAlloc(7))
      Allocate(W(NDet),                 Stat=IAlloc(8))
      Allocate(X(NDet),                 Stat=IAlloc(9))
      Allocate(Y(NDet),                 Stat=IAlloc(10))
      Allocate(Res(NRoots),             Stat=IAlloc(11))
      If(Any(IAlloc /= 0)) Stop "Could not Allocate in Lanczos"


!==================================!
!  This drives the Lanczos sweep.  !
!==================================!

! First step: Orthogonalize
      Call Orthogonalize(CIVecs,NDet,NRoots)


! Initialize logical flags
      Converged = .false.
      Fail = .false.
      OpenClose = NDet > 100000
      If(.not. OpenClose) Open(7,File="Output",Status="Old",Position="Append")

! Initialize the first NRoots elements in q to CIVecs, normalized
! Initialize the first NRoots elements in Hq to H.CIVecs
      q = Zero
      Hq = Zero
      Do iRoot = 1,NRoots
        W = CIVecs(:,iRoot)
        X = W/Sqrt(Dot_Product(W,W))
        Call BuildHC(X,Y,H2aa,H2bb,H2ab,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
        q(:,iRoot)  = X
        Hq(:,iRoot) = Y
      End Do


! Build the first block of Matrix and Metric
      Metric(1:NRoots,1:NRoots) = MatMul(Transpose(q(:,1:NRoots)),q(:,1:NRoots))
      Matrix(1:NRoots,1:NRoots) = MatMul(Transpose(q(:,1:NRoots)),Hq(:,1:NRoots))


! Increase the Krylov space until convergence or something
      Do I = NRoots+1,MaxDim
! Initialize q(I) to Hq(I-1)
        W = Hq(:,I-1)
        X = W

! Orthogonalize q(I) against all previous q, update q, build Hq
        Do J = 1,I-1
          Y = q(:,J)
          Alp = Dot_Product(X,Y)
          W = W - Alp*Y
        End Do
        X = W/Sqrt(Dot_Product(W,W))
        Call BuildHC(X,Y,H2aa,H2bb,H2ab,NOccA,NOccB,NAO,NDetA,NDetB,NDet)
        q(:,I)  = X
        Hq(:,I) = Y


! Augment Matrix and Metric
        Do J = 1,I
          Metric(J,I) = Dot_Product(q(:,J),q(:,I))
          Metric(I,J) = Metric(J,I)
          Matrix(J,I) = Dot_Product(q(:,J),Hq(:,I))
          Matrix(I,J) = Matrix(J,I)
        End Do


! Solve K.U = S.U.E
        Call EigenSolve(Matrix,Metric,Evecs,Evals,I,MaxDim,IFlag)


! Pick up the approriate case
        Select Case(IFlag)
          Case(1)       ! Linear dependency (small eigenvalue of Metric)
            Print *, "Linear dependency in Krylov space - ending sweep early!"
            Exit
          Case(2)       ! Linear dependency (negative eigenvalue of Metric)
            Fail = .true.
            Exit
          Case Default  ! Linearly independent, so update ECI and CIVecs
            ECI = Evals(1:NRoots)
            CIVecs  = MatMul( q(:,1:I),Evecs(1:I,1:NRoots))
            HC = MatMul(Hq(:,1:I),Evecs(1:I,1:NRoots))
        End Select



! Check for convergence on H.CIVecs - ECI*CIVecs
        Do iRoot = 1,NRoots
          X = HC(:,iRoot)
          Y = CIVecs(:,iRoot)
          W = X - ECI(iRoot)*Y
          Res(iRoot) = MaxVal(Abs(W))/MaxVal(Abs(X))
        End Do
        ResNrm = MaxVal(Res)
        Converged = (ResNrm < Tol)
        If(OpenClose) Then
          Open(7,File="Output",Status="Old",Position="Append")
          Write(7,2000) ECI(NRoots), I, Res(NRoots)
          Close(7)
        Else
          Write(7,2000) ECI(NRoots), I, Res(NRoots)
        End If
        If(Converged) Exit
      End Do
      If(.not. OpenClose) Close(7)
2000  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')


!=========================================================!
!  At this point, either we have converged, or we have    !
!  not converged, or the Krylov space has somehow picked  !
!  up a negative eigenvalue.                              !
!=========================================================!

! Deallocate
      DeAllocate(q,        Stat=IAlloc(1))
      DeAllocate(Hq,       Stat=IAlloc(2))
      DeAllocate(Matrix,   Stat=IAlloc(3))
      DeAllocate(Metric,   Stat=IAlloc(4))
      DeAllocate(Evecs,    Stat=IAlloc(5))
      DeAllocate(Evals,    Stat=IAlloc(6))
      DeAllocate(HC,       Stat=IAlloc(7))
      DeAllocate(W,        Stat=IAlloc(8))
      DeAllocate(X,        Stat=IAlloc(9))
      DeAllocate(Y,        Stat=IAlloc(10))
      DeAllocate(Res,      Stat=IAlloc(11))
      If(Any(IAlloc /= 0)) Stop "Could not DeAllocate in Lanczos"

! Distinguish cases
      If(Converged) Return
      If(Fail)      Stop "Negative eigenvalue of Krylov space overlap!"

      Return
      End Subroutine Lanczos






      Subroutine Orthogonalize(V,NDim,N)
      Implicit None
      Integer,        Intent(In)    :: NDim, N
      Real (Kind=pr), Intent(InOut) :: V(NDim,N)
! Local junk
      Real (Kind=pr), Allocatable :: Olap(:,:)
      Real (Kind=pr), Allocatable :: UMat(:,:)
      Real (Kind=pr), Allocatable :: Scr(:,:)
      Real (Kind=pr), Allocatable :: Evals(:)
      Real (Kind=pr), Allocatable :: V0(:,:)
      Integer, Parameter :: Thresh = 1.0E-10_pr
      Integer :: IAlloc, I

!====================!
!  Orthogonalize V.  !
!====================!

! We'll need some space
      Allocate(Olap(N,N),     &
               UMat(N,N),     &
               Scr(N,N),      &
               EVals(N),      &
               V0(NDim,N),    &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in Orthogonalize"


! Build the overlap matrix
      Olap = MatMul(Transpose(V),V)


! Diagonalize it
      Call DiagR(Olap,Evals,UMat,N)


! Stop with an error if the overlap matrix is ill-conditioned
      If(MinVal(Evals) < Zero) Stop "Negative overlap element in initial orthogonalization"
      If(MinVal(Evals)/MaxVal(Evals) < Thresh) Stop "Initial matrix is linearly independent to tolerance"


! Okay.  Build U.s^{-1/2}.U!
      Olap = Zero
      Do I = 1,N
        Olap(I,I) = One/Sqrt(Evals(I))
      End Do
      Scr = MatMul(UMat,Olap) ! U.s^{-1/2}
      Olap = Transpose(UMat)  ! U!
      UMat = MatMul(SCr,Olap) ! U.s^{-1/2}.U!


! Fill in V = V.UMat
      V0 = MatMul(V,UMat)
      V = V0


! Deallocate and exit safely
      DeAllocate(Olap, UMat, Scr, EVals, V0, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in Orthogonalize"

      Return
      End Subroutine Orthogonalize






      Subroutine EigenSolve(Matrix,Metric,EVecs,Evals,N,MaxDim,IFlag)
      Implicit None
      Integer,        Intent(In)  :: N, MaxDim
      Integer,        Intent(Out) :: IFlag
      Real (Kind=pr), Intent(Out) :: EVals(MaxDim)
      Real (Kind=pr), Intent(Out) :: EVecs(MaxDim,MaxDim)
      Real (Kind=pr), Intent(In)  :: Matrix(MaxDim,MaxDim)
      Real (Kind=pr), Intent(In)  :: Metric(MaxDim,MaxDim)
! Local variables
      Real (Kind=pr), Parameter :: LinDepTol = 1.0e-8_pr
      Real (Kind=pr), Allocatable, Dimension(:,:) :: Mat, Olap, UMat, X
      Integer :: I, IAlloc(4)

!================================================!
!  Solves H.C = S.C.E without destryong H or S.  !
!================================================!

! Allocate space
      Allocate(Mat(N,N),    Stat=IAlloc(1))
      Allocate(Olap(N,N),   Stat=IAlloc(2))
      Allocate(UMat(N,N),   Stat=IAlloc(3))
      Allocate(X(N,N),      Stat=IAlloc(4))
      If(Any(IAlloc /= 0)) Stop "Could not Allocate in EigenSolve"


! Build the symmetric orthogonalization matrix X = U.s^(-1/2).U+
      Olap = Metric(1:N,1:N)
      Call DiagR(Olap,Evals(1:N),UMat,N)
      Olap = Zero
      Do I = 1,N
        If(Evals(I) < Zero) Then
          IFlag = 2
          GoTo 20
        Else If(Evals(I) < LinDepTol) Then
          IFlag = 1
          GoTo 20
        End If
        Olap(I,I) = 1.0_pr/Sqrt(Evals(I))
      End Do
      Mat = MatMul(UMat,Olap)  ! Mat  = U.s^{-1/2}
      Olap = Transpose(UMat)   ! Olap = U+
      X = MatMul(Mat,Olap)     ! X = U.s^{-1/2}.U+


! Transform Mat
      Mat  = Matrix(1:N,1:N)
      UMat = MatMul(X,Mat)    ! UMat = X.H
      Mat  = MatMul(UMat,X)   ! Mat  = X.H.X


! Diagonalize Mat
      Call DiagR(Mat,Evals(1:N),UMat,N)
      Olap = MatMul(X,UMat)
      Evecs(1:N,1:N) = Olap
      IFlag = 0


! Exit
20    DeAllocate(Mat,   Stat=IAlloc(1))
      DeAllocate(Olap,  Stat=IAlloc(2))
      DeAllocate(UMat,  Stat=IAlloc(3))
      DeAllocate(X,     Stat=IAlloc(4))
      If(Any(IAlloc /= 0)) Stop "Could not DeAllocate in EigenSolve" 

      Return
      End Subroutine EigenSolve

   End Module CI

