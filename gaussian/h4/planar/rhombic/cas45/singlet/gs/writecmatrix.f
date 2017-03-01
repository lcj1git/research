      program writecmatrix
c writes a cmatrix input file for the rhombic configuration of 4 H atoms
      implicit real*4 (a-h,o-z)
      parameter (rt3 = 1.732051)

c      read(5,*) y
      y = rt3
      open(12, file = 'rhombicH4.com')
      write(12,*) '%chk=rhombicH4.chk'
      write(12,*) '#T CASSCF(4,5)/sto-3g test'
      write(12,*)
      write(12,*) 'mbs CASSCF(4,5) spe calculation for rhombic H4 with'
      write(12,*) 'side length = ', rt3 
      write(12,*) 'short diagonal =', 2.*rt3/y
      write(12,*) 'and long diagonal =', 2.*y
      write(12,*)
      write(12,*) 0, 1
      write(12,*) 'H', 0.0, 0.0, y
      write(12,*) 'H', 0.0, 0.0, -y
      write(12,*) 'H', rt3/y, 0.0, 0.0
      write(12,*) 'H', -rt3/y, 0.0, 0.0
      write(12,*)
      end
