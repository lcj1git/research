%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 8 H atoms (made up of 8 iscosceles triangles)
 
0 1
H 0.0 0.0 0.880
H .656079 0.0 0.0
H -.656079 0.0 0.0
H 0.0 0.0 -0.880
H 1.312158 0.0 -0.880
H -1.312158 0.0 -0.880
H .656079 0.0 -1.760
H -.656079 0.0 -1.760
 
@H.gbs
 
