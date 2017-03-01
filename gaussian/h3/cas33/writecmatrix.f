      program writecmatrix
c writes a cmatrix input file for the rhombic configuration of 4 H atoms
      implicit real*4 (a-h,o-z)
      parameter (rt3 = 1.732051)

      read(5,*) y
      open(12, file = 'h3.com')
      write(12,*) '%chk=h3.chk'
      write(12,*) '#P CASSCF(3,3)/sto-3g pop=naturalorbitals test'
      write(12,*)
      write(12,*)'mbs CASSCF(3,3) spe calculation for isosceles H3 with'
      write(12,*) 'the height of the triangle =',y
      write(12,*) 'and the length of the third side  =', 2.*rt3/y
      write(12,*)
      write(12,*) 0, 2
      write(12,*) 'H', 0.0, 0.0, y
      write(12,*) 'H', rt3/y, 0.0, 0.0
      write(12,*) 'H', -rt3/y, 0.0, 0.0
      write(12,*)
      end
