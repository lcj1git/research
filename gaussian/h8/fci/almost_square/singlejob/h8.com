%chk=h8.chk
 
#t casscf(8,8,nroot=2, fulldiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a triangular lattice of 8 H atoms height=y=0.2900,base=2/rt3y=3.981726
 
0 1
H 0.0 0.0 0.2900
H 1.990863 0.0 0.0
H -1.990863 0.0 0.0
H 0.0 0.0 -0.2900
H 3.981726 0.0 -0.2900
H -3.981726 0.0 -0.2900
H 1.990863 0.0 -.5800
H -1.990863 0.0 -.5800
 
@H.gbs
 
