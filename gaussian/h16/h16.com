%chk=h16.chk
%mem=175MW
 
#t casscf(16,16,nroot=3)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a triangular lattice of 16 H atoms (made up of 18 iscosceles triangles)
 
0 1
H 0.0 0.0 .87
H -1.990863 0.0 .58
H 1.990863 0.0 .58
H -3.981726 0.0 .29
H 0.0 0.0 .29
H 3.981726 0.0 .29
H -5.972589 0.0 0.0
H -1.990863 0.0 0.0
H 1.990863 0.0 0.0
H 5.972589 0.0 0.0
H -3.981726 0.0 -.29
H 0.0 0.0 -.29
H 3.981726 0.0 -.29
H -1.990863 0.0 -.58
H 1.990863 0.0 -.58
H 0.0 0.0 -.87
 
@H.gbs
 
