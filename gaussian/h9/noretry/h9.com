%chk=h9.chk
 
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 9 H atoms (made up of 8 iscosceles triangles h=0.25b=4.618802)
 
0 2
H 0.0 0.0 .50
H 2.309401 0.0 0.25
H -2.309401 0.0 0.25
H 0.0 0.0 0.0
H 4.618802 0.0 0.0
H -4.618802 0.0 0.0
H 0.0 0.0 -.50
H 2.309401 0.0 -0.25
H -2.309401 0.0 -0.25
 
@H.gbs
 
