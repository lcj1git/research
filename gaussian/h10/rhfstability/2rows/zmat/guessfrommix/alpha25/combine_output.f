c23456------------------------------------------------------------------
      program comboutput
      implicit real*8 (a-h,o-z)
      parameter (nlinefwd = 240)
      parameter (nlinebck = 77)
c alpha, y describing triangle geometry
      real*8 alpha(nlinefwd+nlinebck), y(nlinefwd+nlinebck)
c rrhf and crhf hessian eigenvalues
      real*8 rrhf(nlinefwd+nlinebck), crhf(nlinefwd+nlinebck)
c ruhf hessian eigen value
      real*8 ruhf(nlinefwd+nlinebck)
c e of mos
      real*8 emo(nlinefwd+nlinebck,10), egap(nlinefwd+nlinebck)
c rhf e
      real*8 erhf(nlinefwd+nlinebck) 

      nline= nlinefwd + nlinebck      

      open (5,file='fwd/h10hess')
      do i=nlinefwd,1,-1
       read(5,*) alpha(i), y(i), rrhf(i), crhf(i), ruhf(i)
      enddo
      close(5)
      open (5,file='bck/h10hess')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), rrhf(i), crhf(i), ruhf(i)
      enddo
      close(5)
      open(6,file='h10hess')
      do i=1,nline
       write(6,*) alpha(i), y(i), rrhf(i), crhf(i), ruhf(i)
      enddo
      close(6)

      open (5,file='fwd/h10MOE')
      do i=nlinefwd,1,-1
       read(5,*) alpha(i), y(i), (emo(i,j),j=1,10), egap(i)
      enddo
      close(5)
      open (5,file='bck/h10MOE')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), (emo(i,j),j=1,10), egap(i)
      enddo
      close(5)
      open(6,file='h10MOE')
      do i=1,nline
       write(6,*) alpha(i), y(i), (emo(i,j),j=1,10), egap(i)
      enddo
      close(6)

      open (5,file='fwd/h10RHFE')
      do i=nlinefwd,1,-1
       read(5,*) alpha(i), y(i), erhf(i)
      enddo
      close(5)
      open (5,file='bck/h10RHFE')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), erhf(i)
      enddo
      close(5)
      open(6,file='h10RHFE')
      do i=1,nline
       write(6,*) alpha(i), y(i), erhf(i)
      enddo
      close(6)

      end
