 %chk=T4model.chk
 #T ccsd/sto-3g scan test

 T4 model (P4 with one H2 twisted oop) SEP calculation with RHF in the sto-3g basis

 0 1
 H  
 H 1 a  
 H 2 alpha 1 90.0
 H 3 a     2 90.0 1 d

 d    0.0000000              18  5.0
 a    1.6000000    
 alpha    1.6000000    

