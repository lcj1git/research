%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 88888888 H atoms (made up of 8 iscosceles triangles)
 
0 1
H 0.0 0.0 2.999
H .192514 0.0 0.0
H -.192514 0.0 0.0
H 0.0 0.0 -2.999
H .385028 0.0 -2.999
H -.385028 0.0 -2.999
H .192514 0.0 -5.998
H -.192514 0.0 -5.998
 
@H.gbs
 
