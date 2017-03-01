h4.chk
 
#t rhf/gen scf(tight,novaracc,noincfock) geom=nocrowd stable(crhf) test
 
scf stability calculation for a rhombic configuration of 4 H atoms (made up of 2 iscosceles triangles)
 
0 1
H 0.0 0.0 3.475
H .166143 0.0 0.0
H -.166143 0.0 0.0
H 0.0 0.0 -3.475
 
@H.gbs
 
