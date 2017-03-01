%chk=h9.chk
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
preliminary casscf(9,9) calculation to improve doublet convergence at slightly different configuration
y=0.33

0 2
H 0.0 0.0 .66
H -1.749546 0.0 0.33
H 1.749546 0.0 0.33
H -3.499092 0.0 0.0
H 0.0 0.0 0.0
H 3.499092 0.0 0.0
H -1.749546 0.0 -0.33
H 1.749546 0.0 -0.33
H 0.0 0.0 -.66
 
@H.gbs
 
