%chk=h9n3.chk
 
#t casscf(9,9,nroot=3)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a rhombic configuration of 9 H atoms (made up of 8 iscosceles triangles), param y=0.25 try #1
 
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
 
