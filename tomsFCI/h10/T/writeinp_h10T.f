c23456------------------------------------------------------------------
      program writeinp
      implicit real*8 (a-h,o-z)

      ntyp = 0
      nrt = 5
      m = 1
      rrt3=1.0/sqrt(3.0)
      read(5,*) y
      x = rrt3 / y
      open (6,file='Input')
      write(6,*) 'Hamiltonian Type:', ntyp
      write(6,*) 'NRoots:', nrt
      write(6,*) 'Multiplicity:', m
      write(6,*) 'Units:', 'Angstroms'
      write(6,*) ' '
      write(6,*) 'Molecular Geometry:'
      write(6,*) 'Atom  ', 'X   ', 'Y   ', 'Z   ', 'Basis'
      write(6,*) 'H', x, 0.d0, y*2.0, 'CustomH'
      write(6,*) 'H', 0.d0, 0.d0, y, 'CustomH'
      write(6,*) 'H', 2.d0*x, 0.d0, y, 'CustomH'
      write(6,*) 'H', -x, 0.d0, 0.d0, 'CustomH'
      write(6,*) 'H', x, 0.d0, 0.d0, 'CustomH'
      write(6,*) 'H', -2.d0*x, 0.d0, -y, 'CustomH'
      write(6,*) 'H', 0.d0, 0.d0, -y, 'CustomH'
      write(6,*) 'H', 2.d0*x, 0.d0, -y, 'CustomH'
      write(6,*) 'H', -x, 0.d0, -2.0*y, 'CustomH'
      write(6,*) 'H', x, 0.d0, -2.0*y, 'CustomH'
      close(6)

      end
