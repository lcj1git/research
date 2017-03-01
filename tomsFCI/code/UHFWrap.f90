
   Module UHFWrap
   Use Precision
   Use HamiltonianData
   Use UHF
   Use DIIS
   Implicit None
   Contains

      Subroutine DoUHF(Olap,OneH,ERI,EvecA,EvecB,NAO,NOccA,NOccB,ENuc,ESCF)
      Implicit None
      Integer,        Intent(In)  :: NAO, NOccA, NOccB
      Real (Kind=pr), Intent(In)  :: ENuc
      Real (Kind=pr), Intent(In)  :: Olap(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: OneH(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr), Intent(Out) :: EvecA(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: EvecB(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: ESCF
      Integer :: IAlloc, NA, NB, I

!================================================!
!  Hides the UHF section of the code from main.  !
!================================================!

! Set up DIIS
      Call SetUpDIIS(2*NAO*NAO)


!==========================================================================================!
!  What we do here depends on the Hamiltonian type.                                        !
!                                                                                          !
!  For the pairing Hamiltonian, we know the seniority-preserving solution already so we    !
!  build that directly.  If it's repulsive pairing, we could mix orbitals and try again.   !
!                                                                                          !
!  For the molecular Hamiltonian, we do a UHF starting from the core guess, then mix and   !
!  try again.                                                                              !
!                                                                                          !
!  For Hubbard, we do my fancy guess at half filling, occupy the half filling orbitals     !
!  and do another UHF, then mix and try one more time.  It's not ideal!                    !
!==========================================================================================!

      Select Case(HamiltonianType)
        Case(0) ! Molecular Hamiltonian
! First do the core guess...
          Call DrvUHF(Olap,OneH,ERI,EvecA,EvecB,ENuc,ESCF,NOccA,NOccB,NAO,.true.,0)

! Mix orbitals and try again
          Call UHFMix(EvecA,EvecB,NAO,NOccA,NOccB)
          Call DrvUHF(Olap,OneH,ERI,EvecA,EvecB,ENuc,ESCF,NOccA,NOccB,NAO,.true.,1)


        Case(1) ! Pairing Hamiltonian
! Just build the seniority-preserving solution and get the energy.
! Later, we could mix orbitals and try again, but let's skip for now
          EvecA = Zero
          Do I = 1,NAO
            EvecA(I,I) = One
          End Do
          EvecB = EvecA
          Call SimpleUHF(OneH,ERI,EvecA,EvecB,ENuc,ESCF,NOccA,NOccB,NAO)


        Case(2) ! Hubbard Hamiltonian
! Do my fancy guess at half-filling...
          NA = NAO/2
          NB = NAO/2
          Call DrvUHF(Olap,OneH,ERI,EvecA,EvecB,ENuc,ESCF,NA,NB,NAO,.true.,2)

! Occupy the orbitals and try again
          Call DrvUHF(Olap,OneH,ERI,EvecA,EvecB,ENuc,ESCF,NOccA,NOccB,NAO,.true.,1)

! Lastly, mix orbitals and try once more
          Call UHFMix(EvecA,EvecB,NAO,NOccA,NOccB)
          Call DrvUHF(Olap,OneH,ERI,EvecA,EvecB,ENuc,ESCF,NOccA,NOccB,NAO,.true.,1)
      End Select


!=================================!
!  Turn off DIIS and we're done.  !
!=================================!========================================!

      Call ShutDownDIIS

      Return
      End Subroutine DoUHF

   End Module UHFWrap

