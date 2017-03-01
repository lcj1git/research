c23456------------------------------------------------------------------
      program comboutput
      implicit real*8 (a-h,o-z)
      parameter (nlinefwd = 41)
      parameter (nlinebck = 40)
c alpha, y describing triangle geometry
      real*8 alpha(nlinefwd+nlinebck), y(nlinefwd+nlinebck)
c rrhf and crhf hessian eigenvalues
      real*8 rrhf(nlinefwd+nlinebck), crhf(nlinefwd+nlinebck)
c ruhf hessian eigen value, e of 2nd highest occ mo
      real*8 ruhf(nlinefwd+nlinebck), eshomo(nlinefwd+nlinebck)
c e of homo, e of lumo
      real*8 ehomo(nlinefwd+nlinebck), elumo(nlinefwd+nlinebck)
c e of 2nd lowest virt mo, homo lumo e gap
      real*8 eslumo(nlinefwd+nlinebck), egap(nlinefwd+nlinebck)
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
       read(5,*) alpha(i), y(i), eshomo(i), ehomo(i), elumo(i),
     & eslumo(i), egap(i)
      enddo
      close(5)
      open (5,file='bck/h10MOE')
      do i=nlinefwd+1,nline
       read(5,*) alpha(i), y(i), eshomo(i), ehomo(i), elumo(i),
     & eslumo(i), egap(i)
      enddo
      close(5)
      open(6,file='h10MOE')
      do i=1,nline
       write(6,*) alpha(i), y(i), eshomo(i), ehomo(i), elumo(i),
     & eslumo(i), egap(i)
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
