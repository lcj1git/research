%chk=h10.chk
 
#t casscf(10,10,nroot=3,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=0.76b=1.519342)
 
0 1
H .759671 0.0 1.52
H 2.279013 0.0 1.52
H 0.0 0.0 0.76
H 1.519342 0.0 0.76
H -.759671 0.0 0.0
H .759671 0.0 0.0
H -1.519342 0.0 -0.76
H 0.0 0.0 -0.76
H -2.279013 0.0 -1.52
H -.759671 0.0 -1.52
 
@H.gbs
 
