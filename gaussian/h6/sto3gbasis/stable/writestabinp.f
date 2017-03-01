      program writecmatrix
c writes a cmatrix input file for the rhombic configuration of 4 H atoms
      implicit real*4 (a-h,o-z)
      parameter (rt3 = 1.732051)
      character*72 chkfile

      read(5,*) y
      read(5,*) n
      read(5,*) mult
      read(5,*) chkfile
      write(6,*) '%chk=',chkfile
      write(6,*) '#p stable(crhf)/sto-3g guess=mix test'
      write(6,*)
      write(6,*)' mbs rhf crhf stability calculation for a rhombic',
     &          ' configuration of 9 H atoms (made up of',
     &          ' 8 isosceles triangles ) with the height'
      write(6,*)' of each triangle =',y
      write(6,*)' and the length of the third side  =', 2.*rt3/y
      write(6,*)
      write(6,*) 0, mult
      write(6,*) 'H',   0.0,     0.0,  2.*y
      write(6,*) 'H',  rt3/y,    0.0,  y
      write(6,*) 'H', -rt3/y,    0.0,  y
      write(6,*) 'H',  2.*rt3/y, 0.0, 0.0
      write(6,*) 'H', -2.*rt3/y, 0.0, 0.0
      write(6,*) 'H',   0.0,     0.0, 0.0
      write(6,*)
      end
