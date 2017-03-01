      program writecmatrix
c writes a cmatrix input file for the rhombic configuration of 4 H atoms
      implicit real*4 (a-h,o-z)
      parameter (rt3 = 1.732051)

      read(5,*) y
      open(6, file = 'rhombich4.com')
      write(6,*) '%chk=rhombich4.chk'
      write(6,*) '#P casscf(4,4,nroot=5)/sto-3g test'
      write(6,*)
      write(6,*) 'mbs casscf(4,4,nroot=5) spe calculation for rhombic' 
      write(6,*) 'H4 with side length = ', rt3 
      write(6,*) 'short diagonal =', 2.*rt3/y
      write(6,*) 'and long diagonal =', 2.*y
      write(6,*)
      write(6,*) 0, 1
      write(6,*) 'H', 0.0, 0.0, y
      write(6,*) 'H', 0.0, 0.0, -y
      write(6,*) 'H', rt3/y, 0.0, 0.0
      write(6,*) 'H', -rt3/y, 0.0, 0.0
      write(6,*)

      end
