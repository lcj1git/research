%chk=acrolein.chk
#P uhf/sto-3g pop=naturalorbitals test

preliminary uhf calculation to determine active space for cas first excited state calculation

0 3
 C              
 H                  1    1.07000000
 O                  1    1.25840000    2  119.88652650
 C                  1    1.53999999    3  120.22694677    2  180.00000000    0
 H                  4    1.07000000    1  119.88652673    3    0.00000000    0
 C                  4    1.35520000    1  120.22694671    3  180.00000000    0
 H                  6    1.07000001    4  120.22694581    1    0.00000000    0
 H                  6    1.07000000    4  119.88652756    1  180.00000000    0

