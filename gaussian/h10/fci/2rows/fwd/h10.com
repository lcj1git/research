%chk=h10n1.chk
 
#t casscf(10,10,davidsondiag)/gen scf(tight,novaracc,noincfock) guess=read geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=1.735b=.665532)
 
0 1
H .332766 0.0 3.470
H .998298 0.0 3.470
H 0.0 0.0 1.735
H .665532 0.0 1.735
H -.332766 0.0 0.0
H .332766 0.0 0.0
H -.665532 0.0 -1.735
H 0.0 0.0 -1.735
H -.998298 0.0 -3.470
H -.332766 0.0 -3.470
 
@H.gbs
 
