%chk=h9.chk
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) guess=read geom=nocrowd test
 
preliminary casscf(9,9) doublet calculation with initial guess from nearby converged pt
y=0.25
 
0 2
H 0.0 0.0 .50
H -2.309401 0.0 0.25
H 2.309401 0.0 0.25
H -4.618802 0.0 0.0
H 0.0 0.0 0.0
H 4.618802 0.0 0.0
H -2.309401 0.0 -0.25
H 2.309401 0.0 -0.25
H 0.0 0.0 -.50
 
@H.gbs
 
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
preliminary casscf(9,9) calculation to improve doublet convergence at slightly different configuration
y=0.32

0 2
H 0.0 0.0 .64
H -1.804219 0.0 0.32
H 1.804219 0.0 0.32
H -3.608438 0.0 0.0
H 0.0 0.0 0.0
H 3.608438 0.0 0.0
H -1.804219 0.0 -0.32
H 1.804219 0.0 -0.32
H 0.0 0.0 -.64
 
@H.gbs
 
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
preliminary casscf(9,9) calculation to improve doublet convergence at slightly different configuration
y=0.32

0 2
H 0.0 0.0 .64
H -1.804219 0.0 0.32
H 1.804219 0.0 0.32
H -3.608438 0.0 0.0
H 0.0 0.0 0.0
H 3.608438 0.0 0.0
H -1.804219 0.0 -0.32
H 1.804219 0.0 -0.32
H 0.0 0.0 -.64
 
@H.gbs
 
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
preliminary casscf(9,9) calculation to improve doublet convergence at slightly different configuration
y=0.32

0 2
H 0.0 0.0 .64
H -1.804219 0.0 0.32
H 1.804219 0.0 0.32
H -3.608438 0.0 0.0
H 0.0 0.0 0.0
H 3.608438 0.0 0.0
H -1.804219 0.0 -0.32
H 1.804219 0.0 -0.32
H 0.0 0.0 -.64
 
@H.gbs
 
#t casscf(9,9,davidsondiag)/gen scf(tight,novaracc,noincfock) geom=nocrowd test
 
preliminary casscf(9,9) calculation to improve doublet convergence at slightly different configuration
y=0.32

0 2
H 0.0 0.0 .64
H -1.804219 0.0 0.32
H 1.804219 0.0 0.32
H -3.608438 0.0 0.0
H 0.0 0.0 0.0
H 3.608438 0.0 0.0
H -1.804219 0.0 -0.32
H 1.804219 0.0 -0.32
H 0.0 0.0 -.64
 
@H.gbs
 
