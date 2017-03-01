%chk=h10n3.chk
 
#t casscf(10,10,nroot=3,davidsondiag)/gen scf(tight,novaracc,noincfock) guess=read geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=2.185b=.528466)
 
0 1
H -.528466 0.0 2.185
H 0.0 0.0 2.185
H .528466 0.0 2.185
H -.792699 0.0 0.0
H -.264233 0.0 0.0
H .264233 0.0 0.0
H .792699 0.0 0.0
H -.528466 0.0 -2.185
H 0.0 0.0 -2.185
H .528466 0.0 -2.185
 
@H.gbs
 
