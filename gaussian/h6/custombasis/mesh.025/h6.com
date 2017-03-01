h6.chk
 
#t casscf(6,6,nroot=2)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
mbs cascf calculation for a triangular configuration of 6 H atoms (made up of 4 iscosceles triangles)
 
0 3
H 0.0 0.0 .525
H 1.099714 0.0 0.0
H -1.099714 0.0 0.0
H 0.0 0.0 -.525
H -2.199428 0.0 -.525
H 2.199428 0.0 -.525
 
@H.gbs
 
