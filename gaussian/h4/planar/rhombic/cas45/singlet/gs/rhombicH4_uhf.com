 %chk=rhombicH4.chk
 #P uhf/sto-3g pop=naturalorbitals test

 preliminary uhf calculation to explore orbital symmetries for
 mbs CASSCF(4,5) spe calculation for rhombic H4 with
 side length =    1.7320510    
 short diagonal =   2.0000000    
 and long diagonal =   3.4641020    

           0           1
 H   0.0000000       0.0000000       1.7320510    
 H   0.0000000       0.0000000      -1.7320510    
 H   1.0000000       0.0000000       0.0000000    
 H  -1.0000000       0.0000000       0.0000000    

