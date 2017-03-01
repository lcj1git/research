      program writeT4zmatrix
      implicit real*4 (a-h, o-z)
      parameter (pi = 3.14159)



      read(5,*) a
      read(5,*) alpha
      read(5,*) theta_i
      read(5,*) theta_f
      theta_i = theta_i*pi/180.
      theta_f = theta_f*pi/180.
      di = a*sin(theta_i)/2.
      df = a*sin(theta_f)/2.

      open(10, file = 'T4model.com')

      write(10,*) '%chk=T4model.chk'
      write(10,*) '#T ccsd/sto-3g scan test'
      write(10,*)
      write(10,*) 'T4 model (P4 with one H2 twisted oop) SEP calculation
     - with RHF in the sto-3g basis'
      write(10,*)
      write(10,*) '0 1'
      write(10,*)'H ', 0.0, 0.0, 0.0
      write(10,*)'H ', 0.0, 'a', 0.0
      write(10,*)'H ', 'alpha ', 0.0, 'd'
      write(10,*)'H ', 'alpha ', 'a ', '-d'
      write(10,*)
      nstep = 18 
      stepsize = (df-di)/dfloat(nstep) 
      write(10,*) 'd ', di, nstep, stepsize
      write(10,*) 'a ', a
      write(10,*) 'alpha ', alpha
      write(10,*)
      
      close(10)

      end

