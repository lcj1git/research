#T RHF/6-31G(d) Pop=Reg Test

single point E calculation for acetone

0 1
C1
O  C1 CO 
C2 C1 CC O   CCO
C3 C1 CC O  -CCO C2  0.0
H1 C2 CH C1  CCH O   0.0
H2 C2 CH C1  CCH O   DA
H3 C2 CH C1  CCH O  -DA
H4 C3 CH C1  CCH O   0.0
H4 C3 CH C1  CCH O   DA
H5 C3 CH C1  CCH O  -DA
  Variables:
CO  1.07
CH  1.09
CC  1.52
CCO 120.0
CCH 109.5
DA  60.0

