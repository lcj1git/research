 %chk=h12.chk
 #rhf/gen scf(symm,tight,novaracc,noincfock,fulldiag)   stable(crhf) guess=read geom=nocrowd test
 
 hf stability analysis of h10 triangular lattice, y=  0.24633708657984518     
 
 0 1
 H
 H 1   2.3566507346348495     
 H 1   2.3566507346348495       2    168.00000000000000     
 H 2   2.3566507346348495       1    12.000000000000000       3    0.0000000000000000     
 H 4   2.3566507346348495       2    180.00000000000000       3    0.0000000000000000     
 H 5   2.3566507346348495       4    180.00000000000000       3    0.0000000000000000     
 H 3   2.3566507346348495       1    180.00000000000000       4    0.0000000000000000     
 H 7   2.3566507346348495       3    180.00000000000000       4    0.0000000000000000     
 H 3   2.3566507346348495       4    180.00000000000000       1    0.0000000000000000     
 H 7   2.3566507346348495       5    180.00000000000000       3    0.0000000000000000     
 H 4   2.3566507346348495       3    180.00000000000000       1    0.0000000000000000     
 H 5   2.3566507346348495       7    180.00000000000000       3    0.0000000000000000     
 
@H.gbs
 
