c23456------------------------------------------------------------------
      program reorderEs
      implicit real*8 (a-h,o-z)
      parameter(npt=224)
      parameter(nrt=5)
      parameter(pi=3.1415926535898)
      real*8 y(npt), etmp(nrt), e(npt,nrt), alpha(npt)

      open(5,file='Es')
      do i=1,npt
       read(5,*) y(i), (etmp(j),j=1,nrt)
       call reorder(nrt,etmp)
       do j=1,nrt
        e(i,j)=etmp(j)
       enddo
       a=y(i)
       b=1.d0/(a*dsqrt(3.d0))
       alpha(i)=datan(a/b)*180.d0/pi
      enddo
      close(5)

      open(6,file='reorderedEs')
      do i=1,npt
       write(6,*) alpha(i), (e(i,j),j=1,nrt)
      enddo
      close(6)

      end
c23456------------------------------------------------------------------
      subroutine reorder(npt,ar)
      implicit real*8 (a-h,o-z)
      real*8 ar(npt)

      do i=1,npt
       do j=i+1,npt
        if (ar(j) .lt. ar(i)) then
         tmp=ar(i)
         ar(i)=ar(j)
         ar(j)=tmp
        endif
       enddo
      enddo

      return
      end
