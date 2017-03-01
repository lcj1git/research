c23456-----------------------------------------------------------------
      program combineEs
      implicit real*8 (a-h,o-z)
      parameter (nlinesf1 = 100 )
      parameter (nlinesf2 = 100 )
      integer nstates=3
      integer nlines=nlinesf1+nlinesf2
      real*8 e(nstates,nlines)


      open(5,file='bck/h10singlet')
      do i=nlinesf1,1,-1
       read(5,*) e(1,i), e(2,i), e(3,i)
      enddo
      close(5)
      open(5,file='fwd/h10singlet')
      do i=nlinesf1+1,nlines
       read(5,*) e(1,i), e(2,i), e(3,i)
      enddo
      close(5)
      open(6,file='h10singlet')
      do i=1,nlines
       write(6,*) (e(j,i),j=1,3)
      enddo


      end
