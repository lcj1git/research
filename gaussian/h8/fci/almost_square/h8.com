%chk=h8.chk
 
#t casscf(8,8,nroot=3,fulldiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=2.9995b=.384964)
 
0 1
H 0.0 0.0 2.9995
H .192482 0.0 0.0
H -.192482 0.0 0.0
H 0.0 0.0 -2.9995
H .384964 0.0 -2.9995
H -.384964 0.0 -2.9995
H .192482 0.0 -5.9990
H -.192482 0.0 -5.9990
 
@H.gbs
 
