
!============================================!
!  Fac(k)       = k!                         !
!  DFac(k)      = (2k-1)!!                   !
!  Binom(p,q)   = p!/[q! (p-q)!]             !
!  F1(L,eta)    = Int(x^L exp(-eta*x^2) dx)  !
!============================================!


      Function Fac(K)
      Use Precision
      Implicit None
      Integer :: I, K
      Real (Kind=pr) :: Fac
      Fac = 1.0_pr
      If(K==0) Then
       Continue
      Else
       Do I = 1,K
        Fac = Fac*Real(I)
       End Do
      End If
      End Function Fac






      Function DFac(K)
      Use Precision
      Implicit None
      Integer :: K
      Real (Kind=pr) :: Fac, DFac
      DFac = Fac(2*K)/(Fac(K)*2.0_pr**K)
      End Function DFac






      Function Binom(N,K)
      Implicit None
      Integer :: Binom, N, K, I
      Binom = 1
      Do I = 1,K
        Binom = Binom*(N+1-I)/I
      End Do
      Return
      End Function Binom






      Function F1(L,Eta)
      Use Precision
      Use Constants
      Implicit None
      Integer :: L, L2
      Real (Kind=pr) :: Eta, Tmp, F1, DFac
      L2 = L/2
      If(2*L2 /= L) Then
       Tmp = Zero
      Else
       Tmp = Sqrt(Pi/Eta)*DFac(L2)/(Two*Eta)**L2
      End If
      F1 = Tmp
      End Function F1






      Subroutine Timer(OldTime,String)
      Implicit None
      Integer, Intent(InOut) :: OldTime(8)
      Character (len=*)  :: String
      Character (len=8)  :: Date
      Character (len=10) :: Time
      Character (len=5)  :: Zone
      Integer            :: NewTime(8), TimeDiff(8)
      Integer            :: Hours, Minutes
      Real               :: Seconds

!===========================!
!  Writes out timing data.  !
!===========================!

      Call Date_and_Time(Date,Time,Zone,NewTime)
      TimeDiff = NewTime - OldTime

      If(TimeDiff(8) < 0) Then
       TimeDiff(8) = TimeDiff(8) + 1000
       TimeDiff(7) = TimeDiff(7) - 1
      End If

      If(TimeDiff(7) < 0) Then
       TimeDiff(7) = TimeDiff(7) + 60
       TimeDiff(6) = TimeDiff(6) - 1
      End If

      If(TimeDiff(6) < 0) Then
       TimeDiff(6) = TimeDiff(6) + 60
       TimeDiff(5) = TimeDiff(5) - 1
      End If

      TimeDiff(5) = TimeDiff(5)+24*TimeDiff(3)   !  Adding days

      Hours   = TimeDiff(5)
      Minutes = TimeDiff(6)
      Seconds = TimeDiff(7) + Real(TimeDiff(8))/1000
      Write(6,1000) Hours, Minutes, Seconds, String
      Call Date_and_Time(Date,Time,Zone,OldTime)

1000  Format(I3,' hr, ',I2,' min, ',F6.3,' sec for ',A)

      Return
      End Subroutine Timer

