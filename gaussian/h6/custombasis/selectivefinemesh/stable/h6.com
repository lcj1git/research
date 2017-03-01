h6.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 6 H atoms (made up of 4 iscosceles triangles)
 
0 1
H 0.0 0.0 3.49
H .165429 0.0 0.0
H -.165429 0.0 0.0
H 0.0 0.0 -3.49
H -.330858 0.0 -3.49
H .330858 0.0 -3.49
 
@H.gbs
 
