 %chk=V4model.chk
 #T ccsd/sto-3g scan test

 V4 model (maximally twisted T4 with varying r) PES scan with ccsd in the sto-3g basis

 0 1
 H  
 H 1 a  
 H 2 alpha 1 90.0
 H 3 a     2 90.0 1 d

 alpha    1.0 60 0.05
 d    90.000000
 a    1.6000000    

