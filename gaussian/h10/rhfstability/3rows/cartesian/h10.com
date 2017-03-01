%chk=h10.chk
 
# rhf/gen scf(tight,novaracc,noincfock,fulldiag) stable(crhf) guess=mix geom=nocrowd test
 
scf stability calculation for a triangular configuration of 10 H atoms (made up of 8 iscosceles triangles)
 
0 1
H -.463734 0.0 2.49
H 0.0 0.0 2.49
H .463734 0.0 2.49
H -.695601 0.0 0.0
H -.231867 0.0 0.0
H .231867 0.0 0.0
H .695601 0.0 0.0
H -.463734 0.0 -2.49
H 0.0 0.0 -2.49
H .463734 0.0 -2.49
 
@H.gbs
 
