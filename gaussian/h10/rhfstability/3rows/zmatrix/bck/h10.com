 %chk=h10.chk
 #rhf/gen scf(symm,tight,novaracc,noincfock,fulldiag)   stable(crhf) geom=nocrowd test
 
 hf stability analysis of h10 triangular lattice, y=  0.24942049080490034     
 
 0 1
 H
 H 1    2.3281657757612129     
 H 1    2.3281657757612129       2    167.69999999999999     
 H 2    2.3281657757612129       1    12.300000000000001       3    0.0000000000000000     
 H 3    2.3281657757612129       4    180.00000000000000       1    0.0000000000000000     
 H 3    2.3281657757612129       1    180.00000000000000       2    0.0000000000000000     
 H 3    4.6295335826376514       2    180.00000000000000       1    0.0000000000000000     
 H 7    2.3281657757612129       6    180.00000000000000       1    0.0000000000000000     
 H 7    2.3281657757612129       5    180.00000000000000       2    0.0000000000000000     
 H 7    4.6295335826376514       3    180.00000000000000       1    0.0000000000000000     
 
@H.gbs
 
