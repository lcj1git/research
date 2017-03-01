
   Module IO
   Use Precision
   Use Constants
   Use HamiltonianData
   Implicit None
   Integer :: IUseAngstroms = 0
   Private
   Public  :: ReadInput
   Contains

      Subroutine ReadInput(NOccA,NoccB,NAO,NRoots,ENuc)
      Implicit None
! I/O variables
      Integer,        Intent(Out) :: NOccA, NOccB, NAO, NRoots
      Real (Kind=pr), Intent(Out) :: ENuc
! Length of a line and of a possible keyword
      Integer,               Parameter   :: LLine   = 79
      Integer,               Parameter   :: LName   = 19
! List of possible keywords, and are they set?
      Character (len=LName), Allocatable :: ParamName(:)
      Logical,               Allocatable :: SetOnce(:), SetTwice(:)
! Format statement for line/keyword/value, plus...  line/keyword/value
      Character (len=5)     :: FormatString
      Character (len=LLine) :: Line, KeyWord, Value
! Miscellaneous variables
      Logical :: Error, Exists
      Integer :: NParams, IAlloc, I, ExStatus, LineNumber, NOcc
! Things we need to set
      Character (Len=9) :: Units
      Integer :: Charge, Multiplicity

!=====================================================================!
!  This is a pretty complicated subroutine for me, as I know little   !
!  about string-handling.  But here's what everything does.           !
!                                                                     !
!  We open Input and read it one line at a time.  Lines are assumed   !
!  to be either of the general form                                   !
!       Keyword: Value                                                !
!       Keyword:                                                      !
!  We want each of the possible keywords to be set exactly once, so   !
!  we keep track of which are set and of which are set twice.  We     !
!  also check to see if the user tried to set an invalid keyword,     !
!  telling him on which line he tried to do so.                       !
!                                                                     !
!  The parameters NParams, LName, and LLine mean:                     !
!     NParams: The number of keywords we expect to be set.            !
!     LName:   The length of the string naming each keyword.          !
!     LLine:   The length of a line.                                  !
!  Other variables are:                                               !
!     Charge:        The total charge of the system, read in.         !
!     Multiplicity:  The multiplicity of the system, read in.         !
!     CorFunc:       The correlation functional to be used, read in.  !
!     CorPot:        The correlation potential to be used, read in.   !
!     ExFunc:        The exchange functional to be used, read in.     !
!     ExPot:         The exchange potential to be used, read in.      !
!     Exists:        .True. if Input file exists.                     !
!     Error:         .True. if there is an error in the routine.      !
!     FormatString:  If LLine = xyz, FormatString = (Axyz).           !
!     Line:          The line we read in.                             !
!     Keyword:       The keyword we read in.                          !
!     Value:         The value corresponding to that keyword.         !
!     ParamName:     Array of acceptable keyword names.               !
!     SetOnce:       Array, .True. when the right keyword is set.     !
!     SetTwice:      Array, .True. if the keyword is multiply set.    !
!                                                                     !
!  After creating FormatString, we set the defaults for charge,       !
!  multiplicity, and Kohn-Sham stuff.                                 !
!=====================================================================!

      Write(FormatString,'(a2,i2,a1)') '(a',LLine,')'


!==============================================================!
!  If the input file exists, open it and read the first line.  !
!==============================================================!

      Inquire(File='Input',Exist=Exists)
      If(.not. Exists) Stop 'No Input file'

! Read the Hamiltonian Type
      Open(4,File='Input',Status='Old')
      Read(4,FormatString,End=10) Line
      Call ParseLine(Line,KeyWord,Value)
      If(Trim(AdjustL(KeyWord)) /= "Hamiltonian Type") Stop "Hamiltonian Type must be set in the first line of Input!"
      Read(Value,*) HamiltonianType

! Read the number of roots
      Read(4,FormatString,End=10) Line
      Call ParseLine(Line,Keyword,Value)
      If(Trim(AdjustL(Keyword)) /= "NRoots") Stop "NRoots must be set in the second line of Input!"
      Read(Value,*) NRoots


!===================================================================!
!  Given HamiltonianType, we know how many parameters we can have.  !
!  Set NParams and allocate ParamName, SetOnce, SetTwice.           !
!===================================================================!

      Select Case(HamiltonianType)
        Case(0); NParams = 4                          ! Molecular
        Case(1); NParams = 4                          ! Pairing
        Case(2); NParams = 8                          ! Hubbard
        Case Default
          Stop "Unknown Hamiltonian Type in Input!"
      End Select

      Allocate(ParamName(NParams), SetOnce(NParams), SetTwice(NParams), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ReadInput!"


!===========================================!
!  Set SetOnce, SetTwice, and ParamName.    !
!  Also set default values if we have any.  !
!===========================================!

      SetOnce  = .False.
      SetTwice = .False.
      Select Case(HamiltonianType)
        Case(0)
          Charge = 0
          Multiplicity = 1
          Units = "Bohr"
          ParamName = (/'Charge             ',     &
                        'Multiplicity       ',     &
                        'Molecular Geometry ',     &
                        'Units              '/)
        Case(1)
          ENuc = Zero
          ParamName = (/'# of Alpha         ',     &
                        '# of Beta          ',     &
                        '# of Levels        ',     &
                        'G                  '/)
        Case(2)
          ENuc = Zero
          ParamName = (/'# of X Sites       ',     &
                        '# of Y Sites       ',     &
                        '# of Alpha         ',     &
                        '# of Beta          ',     &
                        'T                  ',     &
                        'U                  ',     &
                        'Periodic in X      ',     &
                        'Periodic in Y      '/)
      End Select


!==========================================!
!  Start reading lines and setting stuff.  !
!==========================================!

      LineNumber = 1
      ExStatus = 0
      Do
! Read the line and update the line number
        LineNumber = LineNumber + 1
        Read(4,FormatString,End=10) Line
        Call ParseLine(Line,KeyWord,Value)


! Check if we're setting a parameter for the first time or for a second time
        Do I = 1,NParams
          If(Trim(ParamName(I)) == Trim(AdjustL(KeyWord))) Then
            If(SetOnce(I))  SetTwice(I) = .True.
            If(SetTwice(I)) Write(6,1020) Trim(ParamName(I))
            SetOnce(I) = .True.
          End If
        End Do


! Given the keyword, set the appropriate varaible
        Select Case (Trim(AdjustL(Keyword)))
          Case ('SKIP')

          Case ('Charge')
            Read(Value,*) Charge

          Case ('Multiplicity')
            Read(Value,*) Multiplicity
!           If(Multiplicity /= 1) Stop "Cannot do non-singlets yet!"

          Case ('Molecular Geometry')
            Read(4,FormatString,End=10) Line
            Call ReadGeometry(NAtom,ENuc,NAO)
            LineNumber = LineNumber + NAtom + 2

          Case ('# of X Sites')
            Read(Value,*) NSitesX

          Case ('# of Y Sites')
            Read(Value,*) NSitesY

          Case ('G')
            Read(Value,*) PairingG

          Case ('# of Levels')
            Read(Value,*) NAO

          Case ('# of Alpha')
            Read(Value,*) NOccA

          Case ('# of Beta')
            Read(Value,*) NOccB

          Case ('T')
            Read(Value,*) HubbardT

          Case ('U')
            Read(Value,*) HubbardU

          Case('Periodic in X')
            Read(Value,*) DoPBCX

          Case('Periodic in Y')
            Read(Value,*) DoPBCY

          Case('Units')
            Read(Value,*) Units
            Select Case(Trim(AdjustL(Units)))
              Case('Angstroms'); IUseAngstroms = 1
              Case('Bohr');      IUseAngstroms = 0
              Case Default
                Stop "Unknown units -- Angstroms or Bohr only!"
            End Select

          Case Default   !  Unknown keyword
            ExStatus=1
            Exit
        End Select
      End Do
10    Close(4,Status='Keep')


!==================================================================!
!  If the user tried to set an undefined keyword, tell him where.  !
!  Then check if the required keywords are set.                    !
!  Then check if any are multiply set.                             !
!  Finally, if there was an error, we stop.                        !
!==================================================================!

! Tell the user where he tried to set an undefined keyword
      If(ExStatus == 1) Write(6,1000) LineNumber

! Check if each parameter is set once.
! For the molecular Hamiltonian, we have defaults for Charge and Multiplicty
! If the user failed to set a required parameter, tell him
      If(HamiltonianType == 0) SetOnce(1:2) = .true.
      If(HamiltonianType == 0) SetOnce(5) = .true.
      Do I = 1,NParams
       If(.not. SetOnce(I)) Write(6,1010) Trim(ParamName(I))
      End Do

! Exit with an error if needed
      Error = ExStatus == 1 .or. Any(SetTwice) .or. .not. All(SetOnce)
      If(Error) Stop 'Error in Input'


!=============================================================!
!  Do any Hamiltonian-specific work that still must be done.  !
!=============================================================!

      Select Case(HamiltonianType)
        Case(0)   ! For the molecular Hamiltonian, set NOccA and NOccB
          NOcc = NInt(Sum(AtomCharge)) - Charge
          If(NOcc < 0)         Stop 'Must have positive number of electrons'
          NOccA = (NOcc + Multiplicity - 1)/2
          NOccB = (NOcc - Multiplicity + 1)/2
          If(NOccA + NOccB /= NOcc) Stop "Impossible multiplicity!"

        Case(2)   ! For the Hubbard Hamiltonian, set NAO
          NAO = NSitesX*NSitesY
      End Select


!=========================================!
!  Finally, tell the user what he's got.  !
!=========================================!

      Select Case(HamiltonianType)
        Case(0)
          Call BasicDataOutput(NAO,NAtom,NOcc)
        Case(1)
          Open(7,File="Output",Status="Replace")
            Write(7,2000) NOccA, NOccB, NAO
            Write(7,2001) PairingG
          Close(7)
        Case(2)
          Open(7,File='Output',Status='Replace')
            Select Case(DoPBCX)
              Case(.true.);  Write(7,2100) NSitesX
              Case(.false.); Write(7,2110) NSitesX
            End Select
            Select Case(DoPBCY)
              Case(.true.);  Write(7,2200) NSitesY
              Case(.false.); Write(7,2210) NSitesY
            End Select
            Write(7,2300) NOccA, NOccB, HubbardT, HubbardU
          Close(7)
      End Select

! Generic format statements
1000  Format('Error: Line number ',I4,' of Input not recognized; ',  &
             'subsequent lines not read.')
1010  Format('Error: Parameter ',A,' not set.')
1020  Format('Error: Parameter ',A,' multiply set.')

! Format statements for the pairing Hamiltonian
2000  Format("There are ",I3," alpha electrons and ",I3," beta electrons in ",I3," levels.")
2001  Format("G = ",F8.5)

! Format statements for the Hubbard Hamiltonian
2100  Format(3x,"The system is periodic     in x with ",I3," sites.")
2110  Format(3x,"The system is non-periodic in x with ",I3," sites.")
2200  Format(3x,"The system is periodic     in y with ",I3," sites.")
2210  Format(3x,"The system is non-periodic in y with ",I3," sites.")
2300  Format(3x,"There are ",I3," spin-up electrons and ",    &
                             I3," spin-down electrons.",/,    &
             3x,"T = ",F18.10,/,3x,"U = ",F18.10,/)

      Return
      End Subroutine ReadInput






      Subroutine ParseLine(Line,KeyWord,Value)
      Implicit None
      Character (Len=*),         Intent(In)  :: Line
      Character (Len=Len(Line)), Intent(Out) :: KeyWord, Value
      Character (Len=Len(Line))              :: WorkingLine
      Integer :: I, CommentPos, ColonPos

!==================================================!
!  Swiped from Paul, reads keywords and values.    !
!--------------------------------------------------!
!  Check for comments and convert them to spaces,  !
!  then turns control characters into spaces       !
!==================================================!

      WorkingLine = Line
      CommentPos  = Index(WorkingLine,'!')
      If(CommentPos /= 0) WorkingLine = WorkingLine(1:CommentPos-1)
      Do I = 1,Len(Line)
       If(IAChar(WorkingLine(I:I)) <= 31) WorkingLine(I:I) = ' '
      End Do


!=================================================================!
!  If we have a blank line at this point, we skip this line.      !
!  Otherwise, we presumably have Keyword: Value and return them.  !
!=================================================================!

      If(Len_Trim(WorkingLine) == 0) Then
       KeyWord = 'SKIP'
       Return
      End If
      ColonPos = Index(WorkingLine,':',back=.False.)
      KeyWord = WorkingLine(1:ColonPos-1)
      Value   = WorkingLine(ColonPos+1:)

      Return
      End Subroutine ParseLine






      Subroutine ReadGeometry(NAtom,ENuc,NAO)
      Implicit None
      Integer :: NAtom, NAO
      Real (Kind=pr) :: ENuc
      Integer, Parameter           :: LLine = 79
      Character (len=LLine)        :: Line,Keyword,Value
      Character (len=2)            :: Element
      Character (len=5)            :: FormatString
      Character (len=12)           :: Basis
      Character (len=12), External :: AllCase
      Real (Kind=pr), Dimension(3) :: R
      Integer,        Dimension(3) :: IAlloc
      Integer :: I, J
      Real (Kind=pr) :: Rij

!===================================!
!  Read in the molecular geometry.  !
!===================================!

      Write(FormatString,'(a2,i2,a1)') '(a',LLine,')'

      NAtom = 0
      Open(4,File='Input',Position='AsIs')
      Do
       Read(4,FormatString,End=10) Line
       Call ParseLine(Line,Keyword,Value)
       If(KeyWord == 'SKIP') GoTo 10 
       NAtom = NAtom + 1
      End Do
10    Continue

      Do I = 1,NAtom+1
       BackSpace 4
      End Do

!=================================================!
!  Now we know how many atoms we have, so we can  !
!  allocate and fill in the atomic data arrays.   !
!=================================================!

      Call SetUpGeometry(NAtom)
      Do I = 1,NAtom
       Read(4,*) Element, R, Basis
       AtomName(I)     = Element
       AtomCenter(I,:) = R
       AtomBasis(I)    = Basis
       Select Case(Element)
        Case('H','h')
         AtomCharge(I) = 1.0_pr
        Case('HE','He','he')
         AtomCharge(I) = 2.0_pr
        Case('LI','Li','li')
         AtomCharge(I) = 3.0_pr
        Case('BE','Be','be')
         AtomCharge(I) = 4.0_pr
        Case('B','b')
         AtomCharge(I) = 5.0_pr
        Case('C','c')
         AtomCharge(I) = 6.0_pr
        Case('N','n')
         AtomCharge(I) = 7.0_pr
        Case('O','o')
         AtomCharge(I) = 8.0_pr
        Case('F','f')
         AtomCharge(I) = 9.0_pr
        Case('NE','Ne','ne')
         AtomCharge(I) = 10.0_pr
        Case('NA','Na','na')
         AtomCharge(I) = 11.0_pr
        Case('MG','Mg','mg')
         AtomCharge(I) = 12.0_pr
        Case('AL','Al','al')
         AtomCharge(I) = 13.0_pr
        Case('SI','Si','si')
         AtomCharge(I) = 14.0_pr
        Case('P','p')
         AtomCharge(I) = 15.0_pr
        Case('S','s')
         AtomCharge(I) = 16.0_pr
        Case('CL','Cl','cl')
         AtomCharge(I) = 17.0_pr
        Case('AR','Ar','ar')
         AtomCharge(I) = 18.0_pr
        Case('X','x','GH','Gh','gh')
         AtomCharge(I) = 0.0_pr
        Case Default
         Stop 'Unknown element in Input; we only know up to Argon'
       End Select
      End Do

! If units are in angstroms, convert to bohr here
      If(IUseAngstroms == 1) AtomCenter = AtomCenter*1.88972613_pr

      ENuc  = Zero
      Do I = 1,NAtom
       Do J = I+1,NAtom
        Rij = Sqrt(Dot_Product(AtomCenter(I,:) - AtomCenter(J,:),       &
                               AtomCenter(I,:) - AtomCenter(J,:)))
        ENuc = ENuc + AtomCharge(I)*AtomCharge(J)/Rij
       End Do
      End Do


!=======================================================================!
!  Now that we've got the geometry and know WHICH basis each atom has,  !
!  we can pick up the information ABOUT that basis, then set the rest   !
!  of our variables and normalize the contracted AOs.                   !
!=======================================================================!

      Call GetBasisInfo(NAtom,NAO)
      Call Normalize(NAO)
      Do I = 1,NAO
       AOCenter(I,:) = AtomCenter(IAOAtom(I),:)
      End Do

      Return
      End Subroutine ReadGeometry






      Subroutine GetBasisInfo(NAtom,NAO)
      Implicit None
      Integer :: NAtom, NAO
      Integer, Parameter    :: LLine = 79
      Character (len=LLine) :: Line,Keyword,Value
      Character (len=5)     :: FormatString
      Logical               :: SetVars, FoundBas, RightAtom, RightBasis
      Logical               :: Exists
      Character (len=12), External :: AllCase
      Integer :: JPrim, JCont, IAtom, IDummy

!============================!
!  Reads in the basis info.  !
!============================!

      Write(FormatString,'(a2,i2,a1)') '(a',LLine,')'
      Inquire(File='Basis',Exist=Exists)
      If(.not. Exists) Stop 'Basis file does not exist'
      Open(55,File='Basis',Status='Old')


!===============================================!
!  For each atom, read through the file Basis   !
!  until we come to the correct atom:basis      !
!  combination.  Pick up the number of AOs for  !
!  that atom.  If we don't find the pair, we    !
!  exit with an error.                          !
!===============================================!

      NAO   = 0
      IDummy = 0
      SetVars = .false.

      Do IAtom = 1,NAtom
       FoundBas = .false.
       Rewind(55)
       Do
        Read(55,FormatString,End=10) Line
        Call ParseLine(Line,KeyWord,Value)
        RightAtom  = Trim(AdjustL(AllCase(Keyword))) ==                 &
                     Trim(AdjustL(AllCase(AtomName(IAtom))))
        RightBasis = Trim(AdjustL(AllCase(Value)))   ==                 &
                     Trim(AdjustL(AllCase(AtomBasis(IAtom))))
        If(RightAtom .and. RightBasis) Then
         FoundBas = .true.
         Call ReadBas(SetVars,NAO,IDummy,IAtom)
         Exit
        End If
       End Do
10     Continue
       If(.not. FoundBas) Then
        Write(6,1000) AtomName(IAtom), AtomBasis(IAtom)
        Stop 'Could not find basis set'
       End If
      End Do
      NPrim = IDummy


!================================================!
!  Now we have the number of basis functions so  !
!  we read in the actual data after allocation.  !
!================================================!

      JPrim = 0
      JCont = 0
      SetVars = .true.

      Call SetUpBasis(NAO)
      Do IAtom = 1,NAtom
       Rewind(55)
       Do
        Read(55,FormatString) Line
        Call ParseLine(Line,KeyWord,Value)
        RightAtom  = Trim(AdjustL(AllCase(Keyword))) ==                 &
                     Trim(AdjustL(AllCase(AtomName(IAtom))))
        RightBasis = Trim(AdjustL(AllCase(Value)))   ==                 &
                     Trim(AdjustL(AllCase(AtomBasis(IAtom))))
        If(RightAtom .and. RightBasis) Then
         Call ReadBas(SetVars,JCont,JPrim,IAtom)
         Exit
        End If
       End Do
      End Do
      Close(55)

1000  Format('Could not find basis ',A,': ',A)

      Return
      End Subroutine GetBasisInfo






      Subroutine ReadBas(SetVars,JCont,JPrim,IAtom)
      Implicit None
      Logical :: SetVars
      Integer :: JPrim, JCont, IAtom, LowerBound, UpperBound
      Integer,        Dimension(:), Allocatable :: LShell, NContShell
      Real (Kind=pr), Dimension(:), Allocatable :: Alpha, Coeff
      Real (Kind=pr) :: A, DNrm, F1
      Integer :: NShells, IShell, ICont, LTot, IPrim, NPrim_ICont
      Integer :: K, L, M, N

!================================================================!
!  GetBasisInfo reads the Basis file.  This routine, once we've  !
!  got the right basis set for a given atom, fills in the right  !
!  section of the various arrays in BasisData.  We also use it   !
!  to find out how many primitive and contracted AOs we have.    !
!================================================================!

      Open(55,Position='AsIs')
      Read(55,*) NShells
      Allocate(LShell(NShells), NContShell(NShells))
      Read(55,*) LShell
      Read(55,*) NContShell

      Do IShell = 1,NShells
       LTot = LShell(IShell)
       Do ICont = 1,NContShell(IShell)
        Read(55,*) NPrim_ICont
        Allocate(Alpha(NPrim_ICont), Coeff(NPrim_ICont))
        Read(55,*) Alpha, Coeff
        K = 0
        Do L = LTot,0,-1
         Do M = LTot-L,0,-1
          N = LTot - L - M
          K = K + 1
          JCont = JCont + 1
          If(SetVars) Then
           IAOAtom(JCont) = IAtom
           LMNAO(JCont,:) = (/L,M,N/)
           LowerBound     = JPrim + (K-1)*NPrim_ICont + 1
           UpperBound     = JPrim + K*NPrim_ICont
           FirstPrimInCont(JCont) = LowerBound
           LastPrimInCont(JCont)  = UpperBound
           AOExp  (LowerBound:UpperBound) = Alpha
           AOCoeff(LowerBound:UpperBound) = Coeff
           Do IPrim = LowerBound,UpperBound
            A = AOExp(IPrim)
            DNrm = F1(2*L,2*A)*F1(2*M,2*A)*F1(2*N,2*A)
            AOCoeff(IPrim) = AOCoeff(IPrim)/Sqrt(DNrm)
           End Do
          End If
         End Do
        End Do
        JPrim = JPrim + NPrim_ICont*(LTot+1)*(LTot+2)/2
        DeAllocate(Alpha,Coeff)
       End Do
      End Do
      DeAllocate(LShell,NContShell)

      Return
      End Subroutine ReadBas






      Subroutine Normalize(NAO)
      Implicit None
      Integer, Intent(In) :: NAO
      Integer :: IBas, L, M, N, I, J, IStart, IEnd
      Real (Kind=pr) :: DNorm, A, C, F1

!=======================!
!  Normalizes the AOs.  !
!=======================!

      Do IBas = 1,NAO
       DNorm = Zero
       IStart = FirstPrimInCont(IBas)
       IEnd   = LastPrimInCont(IBas)
       L      = LMNAO(IBas,1)
       M      = LMNAO(IBas,2)
       N      = LMNAO(IBas,3)
       Do I = IStart,IEnd
       Do J = IStart,IEnd
        A = (AOExp(I)   + AOExp(J))/Two
        C = AOCoeff(I) * AOCoeff(J)
        DNorm = DNorm + C*F1(2*L,2*A)*F1(2*M,2*A)*F1(2*N,2*A)
       End Do
       End Do
       AONorm(IBas) = One/Sqrt(DNorm)
      End Do

      Return
      End Subroutine Normalize






      Subroutine BasicDataOutput(NAO,NAtom,NOcc)
      Implicit None
      Integer, Intent(In) :: NAO, NAtom, NOcc
      Integer, Parameter :: LLine = 79
      Character (len=5)  :: FormatString
      Character (len=50) :: String, Tmp
      Character (len=1)  :: AngMom
      Character (len=2)  :: N
      Logical            :: RightAtom, RightBasis, BasEqual, AtomEqual
      Integer, Dimension(:), Allocatable :: LShell, NContShell
      Character (len=12), External       :: AllCase
      Character (len=LLine)              :: Line,Keyword,Value
      Integer :: I, J, NShells

!=============================================================!
!  Here, we're just displaying the input data in a nice way.  !
!=============================================================!

      Write(FormatString,'(a2,i2,a1)') '(a',LLine,')'
      Open(55,File='Basis',Status='Old')
      Rewind(55)

      Open(7,File='Output',Status='Replace')
      Write(7,1000) NAO,NOcc,NAtom
      Write(7,*)
      Write(7,1010)

      Do I = 1,NAtom
       Write(7,1020) AtomName(I), AtomCenter(I,:), AtomBasis(I)
      End Do

      Write(7,*)
      Do I = 1,NAtom
       Rewind(55)
       BasEqual = .false.
       AtomEqual = .false.
       Do J = 1,I-1
        If(Trim(AdjustL(AllCase(AtomBasis(J)))) ==    &
           Trim(AdjustL(AllCase(AtomBasis(I))))) BasEqual = .true.
        If(Trim(AdjustL(AllCase(AtomName(J)))) ==     &
           Trim(AdjustL(AllCase(AtomName(I))))) AtomEqual = .true.
       End Do
       If(BasEqual .and. AtomEqual) Cycle
       Do
        Read(55,FormatString,End=10) Line
        Call ParseLine(Line,KeyWord,Value)
        RightAtom  = Trim(AdjustL(AllCase(Keyword))) ==                 &
                     Trim(AdjustL(AllCase(AtomName(I))))
        RightBasis = Trim(AdjustL(AllCase(Value)))   ==                 &
                     Trim(AdjustL(AllCase(AtomBasis(I))))
        If(.not. (RightAtom .and. RightBasis)) Cycle
        Read(55,*) NShells
        Allocate(LShell(NShells), NContShell(NShells))
        Read(55,*) LShell
        Read(55,*) NContShell
        String = Repeat(' ',50)
        Do J = 1,NShells
         Write(N,'(I2)') NContShell(J)
         If(LShell(J) == 0) AngMom = 's'
         If(LShell(J) == 1) AngMom = 'p'
         If(LShell(J) == 2) AngMom = 'd'
         If(LShell(J) >= 3) AngMom = AChar(99 + LShell(J))
         Tmp = Trim(AdjustL(String))//' '//Trim(N)//AngMom
         String = AdjustL(Tmp)
        End Do
        DeAllocate(LShell,NContShell)
        Exit
       End Do
       Write(7,1030) AtomName(I), AtomBasis(I), Trim(String)
      End Do
10    Close(55,Status='Keep')
      Write(7,*)
      Close(7)


1000  Format(3x,"There are ",I3," basis functions for ",I3,         &
                " electrons on ",I3," centers.")
1010  Format("Atom      X        Y       Z       Basis")
1020  Format(A2,5x,F7.3,2x,F7.3,2x,F7.3,3x,A)
1030  Format(A,':',A,4x,A)

      Return
      End Subroutine BasicDataOutput

   End Module IO

