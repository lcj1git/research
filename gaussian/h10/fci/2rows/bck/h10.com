%chk=h10n1.chk
 
#t casscf(10,10,davidsondiag)/gen scf(tight,novaracc,noincfock) guess=read geom=nocrowd test
 
mbs casscf calculation for a triangular lattice of 8 H atoms (made up of 10 iscosceles triangles h=.335b=3.446866)
 
0 1
H 1.723433 0.0 .670
H 5.170299 0.0 .670
H 0.0 0.0 .335
H 3.446866 0.0 .335
H -1.723433 0.0 0.0
H 1.723433 0.0 0.0
H -3.446866 0.0 -.335
H 0.0 0.0 -.335
H -5.170299 0.0 -.670
H -1.723433 0.0 -.670
 
@H.gbs
 
