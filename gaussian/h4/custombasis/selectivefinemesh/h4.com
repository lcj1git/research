h4.chk
 
#t casscf(4,4,nroot=5)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a rhombic configuration of 4 H atoms (made up of 2 iscosceles triangles)
 
0 3
H 0.0 0.0 2.99
H .193093 0.0 0.0
H -.193093 0.0 0.0
H 0.0 0.0 -2.99
 
@H.gbs
 
