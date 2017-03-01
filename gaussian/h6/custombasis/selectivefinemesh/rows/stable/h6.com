h6.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for 6 H atoms arranged in two rows (made up of 4 iscosceles triangles)
 
0 1
H 0.0 0.0 2.995
H .192771 0.0 0.0
H -.192771 0.0 0.0
H 0.0 0.0 -2.995
H -.385542 0.0 -2.995
H -.192771 0.0 -5.990
 
@H.gbs
 
