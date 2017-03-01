%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 88888888 H atoms (made up of 8 iscosceles triangles)
 
0 1
H -4.528236 0.0 .255
H 0.0 0.0 .255
H 4.528236 0.0 .255
H -2.264118 0.0 0.0
H 2.264118 0.0 0.0
H -4.528236 0.0 -.255
H 0.0 0.0 -.255
H 4.528236 0.0 -.255
 
@H.gbs
 
