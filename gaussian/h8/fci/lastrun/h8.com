%chk=h8.chk
 
#t casscf(8,8,nroot=3)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a triangular lattice of 8 H atoms (made up of 6 iscosceles triangles)
 
0 1
H 0.0 0.0 2.99
H .193093 0.0 0.0
H -.193093 0.0 0.0
H 0.0 0.0 -2.99
H .386186 0.0 -2.99
H -.386186 0.0 -2.99
H .193093 0.0 -5.98
H -.193093 0.0 -5.98
 
@H.gbs
 
