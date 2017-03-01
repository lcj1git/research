%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 88888888 H atoms (made up of 8 iscosceles triangles)
 
0 1
H -1.360070 0.0 .849
H 0.0 0.0 .849
H 1.360070 0.0 .849
H -.680035 0.0 0.0
H .680035 0.0 0.0
H -1.360070 0.0 -.849
H 0.0 0.0 -.849
H 1.360070 0.0 -.849
 
@H.gbs
 
