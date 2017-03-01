 %chk=rhombich4.chk
 #P casscf(4,4,nroot=3)/sto-3g test

 mbs casscf(4,4,nroot=3) spe calculation for rhombic
 H4 with side length =    1.7320510    
 short diagonal =  0.88823128    
 and long diagonal =   7.8000002    

           0           1
 H   0.0000000       0.0000000       3.9000001    
 H   0.0000000       0.0000000      -3.9000001    
 H  0.44411564       0.0000000       0.0000000    
 H -0.44411564       0.0000000       0.0000000    

