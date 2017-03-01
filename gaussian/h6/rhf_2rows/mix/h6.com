 %chk=h6.chk
 #rhf/gen scf(tight,novaracc,noincfock,fulldiag)        stable(crhf) guess=mix geom=nocrowd test
 
 hf stability analysis of h6 triangular lattice, y=   2.5066432657105051      alpha=   84.750000000000000     
 
 0 1
 H
 H 1   2.5172031055669808     
 H 1   2.5172031055669808       2    10.500000000000000     
 H 2   2.5172031055669808       1    169.50000000000000       3    0.0000000000000000     
 H 3   2.5172031055669808       1    180.00000000000000       2    0.0000000000000000     
 H 4   2.5172031055669808       2    180.00000000000000       1    0.0000000000000000     
 
@H.gbs
 
