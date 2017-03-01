%chk=H4model.chk
#T ccsd/3-21g scan test

pes of trapezoidal model of four H atoms scanned over varying values of angle delta, using ccsd in 3-21g basis

0 1
 H              
 H  1  a
 H  1  a  2  deltaPlus2pi
 H  2  a  1  deltaPlus2pi  3  0.0

deltaPLus2pi 90.0 18 5.0
a            1.6

