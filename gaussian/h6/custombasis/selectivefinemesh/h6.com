h6.chk
 
# casscf(6,6,nroot=3,fulldiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a triangular configuration of 6 H atoms (made up of 4 iscosceles triangles)
 
0 1
H 0.0 0.0 2.995
H .192771 0.0 0.0
H -.192771 0.0 0.0
H 0.0 0.0 -2.995
H .385542 0.0 -2.995
H -.385542 0.0 -2.995
 
@H.gbs
 
