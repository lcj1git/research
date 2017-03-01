%chk=S4model.chk
# UMP4/sto-3g scan test

P4 model PES scan over varying intermolecular H2 distances in the 3-21g basis

0 1
 H              
 H                  1    R
 H                  1    1.6   2   90.00000000
 H                  2    R     1   90.00000000 3 0.0

R 1.0 60 0.05

