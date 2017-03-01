%chk=h8.chk
 
#t casscf(8,8,nroot=2,fulldiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=2.560b=.451054)
 
0 1
H .225527 0.0 5.120
H 0.0 0.0 2.560
H .451054 0.0 2.560
H -.225527 0.0 0.0
H .225527 0.0 0.0
H -.451054 0.0 -2.560
H 0.0 0.0 -2.560
H -.225527 0.0 -5.120
 
@H.gbs
 
