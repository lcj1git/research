 %chk=h10.chk
 #uhf/gen scf(tight,novaracc,noincfock,fulldiag)        stable(crhf) guess=read geom=nocrowd test
 
 hf stability analysis of h10 triangular lattice, y=  0.47690247052288248      alpha=   21.501000000000001     
 
 0 1
 H
 H 1   1.3011725299533654     
 H 1   1.3011725299533654       2    136.99799999999999     
 H 2   1.3011725299533654       1    43.002000000000002       3    0.0000000000000000     
 H 3   1.3011725299533654       1    180.00000000000000       2    0.0000000000000000     
 H 4   1.3011725299533654       2    180.00000000000000       1    0.0000000000000000     
 H 5   1.3011725299533654       3    180.00000000000000       2    0.0000000000000000     
 H 6   1.3011725299533654       4    180.00000000000000       1    0.0000000000000000     
 H 7   1.3011725299533654       5    180.00000000000000       2    0.0000000000000000     
 H 8   1.3011725299533654       6    180.00000000000000       1    0.0000000000000000     
 
@H.gbs
 
