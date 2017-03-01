%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) geom=nocrowd test
 
scf stability calculation for a triangular configuration of 8 H atoms
 
0 1
H .193093 0.0 5.98
H 0.0 0.0 2.99
H .386186 0.0 2.99
H -.193093 0.0 0.0
H .193093 0.0 0.0
H -.386186 0.0 -2.99
H 0.0 0.0 -2.99
H -.193093 0.0 -5.98
 
@H.gbs
 
