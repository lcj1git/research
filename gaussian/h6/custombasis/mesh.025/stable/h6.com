h6.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 6 H atoms (made up of 4 iscosceles triangles)
 
0 1
H 0.0 0.0 3.475
H .166143 0.0 0.0
H -.166143 0.0 0.0
H 0.0 0.0 -3.475
H -.332286 0.0 -3.475
H .332286 0.0 -3.475
 
@H.gbs
 
