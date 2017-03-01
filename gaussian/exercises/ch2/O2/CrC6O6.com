#T RHF/3-21G  Test

geometric optimization of chromium hexacarbonyl

0 1
Cr
C1 Cr CCr
O1 C1 CO  Cr 180.0
C2 Cr CCr C1 180.0 O1  0.0
O2 C2 CO  Cr 180.0 C1  0.0
C3 Cr CCr C1  90.0 C2  0.0
O3 C3 CO  Cr 180.0 C1  0.0
C4 Cr CCr C3 180.0 C2  0.0
O4 C4 CO  Cr 180.0 C1  0.0
C5 Cr CCr C3  90.0 C4 90.0
  Variables:
CCr 1.94
CO  1.14

