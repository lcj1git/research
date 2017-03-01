%chk=h10.chk
 
#t casscf(10,10)/sto-3g scf(tight,novaracc,noincfock) guess=read geom=nocrowd test
 
mbs cascf calculation for a triangular lattice of 8 H atoms height=y=0.3624,base=2/rt3y=3.186258
 
0 1
H -3.186258 0.0 0.3624
H 0.0 0.0 0.3624
H 3.186258 0.0 0.3624
H -4.779387 0.0 0.0
H -1.593129 0.0 0.0
H 1.593129 0.0 0.0
H 4.779387 0.0 0.0
H -3.186258 0.0 -0.3624
H 0.0 0.0 -0.3624
H 3.186258 0.0 -0.3624
 
