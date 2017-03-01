 %chk=h3.chk
 #P CASSCF(3,3)/sto-3g pop=naturalorbitals test

 mbs CASSCF(3,3) spe calculation for isosceles H3 with
 the height of the triangle =   1.7300000    
 and the length of the third side  =   2.0023711    

           0           2
 H   0.0000000       0.0000000       1.7300000    
 H   1.0011855       0.0000000       0.0000000    
 H  -1.0011855       0.0000000       0.0000000    

