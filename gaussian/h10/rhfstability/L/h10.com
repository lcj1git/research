%chk=h10.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 10 H atoms (made up of 8 iscosceles triangles)
 
0 1
H -.231867 0.0 4.98
H -.463734 0.0 2.49
H 0.0 0.0 2.49
H .463734 0.0 2.49
H -.231867 0.0 0.0
H .231867 0.0 0.0
H .695601 0.0 0.0
H 0.0 0.0 -2.49
H .463734 0.0 -2.49
H .231867 0.0 -4.98
 
@H.gbs
 
