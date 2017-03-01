c23456------------------------------------------------------------------
      program comboutput
      implicit real*8 (a-h,o-z)
      parameter (nlinefwd = 190)
      parameter (nlinebck = 130)
c alpha, y describing triangle geometry
      real*8 alpha(nlinefwd+nlinebck), y(nlinefwd+nlinebck)
c rrhf and crhf hessian eigenvalues
      real*8 uhfi(nlinefwd+nlinebck), cuhf(nlinefwd+nlinebck,2)
c ruhf hessian eigen value, e of 2nd highest occ mo
      real*8 ghf(nlinefwd+nlinebck)
c uhf e
      real*8 euhf(nlinefwd+nlinebck) 

      nline= nlinefwd + nlinebck      

      open (5,file='fwd/h10hess')
      do i=nlinefwd,1,-1
       read(5,*) alpha(i), y(i), uhfi(i), cuhf(i,1), cuhf(i,2), ghf(i)
      enddo
      close(5)
      open (5,file='bck/h10hess')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), uhfi(i), cuhf(i,1), cuhf(i,2), ghf(i)
      enddo
      close(5)
      open(6,file='h10hess')
      do i=1,nline
       write(6,*) alpha(i), y(i), uhfi(i), cuhf(i,1), cuhf(i,2), ghf(i)
      enddo
      close(6)

      open (5,file='fwd/h10UHFE')
      do i=nlinefwd,1,-1
       read(5,*) alpha(i), y(i), euhf(i)
      enddo
      close(5)
      open (5,file='bck/h10UHFE')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), euhf(i)
      enddo
      close(5)
      open(6,file='h10UHFE')
      do i=1,nline
       write(6,*) alpha(i), y(i), euhf(i)
      enddo
      close(6)

      end
