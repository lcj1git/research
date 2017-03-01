%chk=h8.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=read guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 8 H atoms (made up of 8 iscosceles triangles)
 
0 1
H .208279 0.0 5.544
H 0.0 0.0 2.772
H .416558 0.0 2.772
H -.208279 0.0 0.0
H .208279 0.0 0.0
H -.416558 0.0 -2.772
H 0.0 0.0 -2.772
H -.208279 0.0 -5.544
 
@H.gbs
 
