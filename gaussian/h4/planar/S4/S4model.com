%chk=S4model.chk
#T CASSCF(4,4)/sto-3g scan test

S4 model PES scan with varying H2 bond length using UMP4 with mbs

0 1
 H              
 H                  1    R
 H                  1    R   2   90.00000000
 H                  2    R   1   90.00000000 3 0.0

R 1.0 60 0.05

