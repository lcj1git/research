%chk=h8.chk
 
#t casscf(8,8,nroot=3,fulldiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=2.999b=.385028)
 
0 1
H -.385028 0.0 2.999
H 0.0 0.0 2.999
H .385028 0.0 2.999
H -.192514 0.0 0.0
H .192514 0.0 0.0
H -.385028 0.0 -2.999
H 0.0 0.0 -2.999
H .385028 0.0 -2.999
 
@H.gbs
 
