c23456-----------------------------------------------------------------
      program writez
      implicit real*8 (a-h, o-z)
      parameter(pi=3.1415926535898)

      area = 1.d0/dsqrt(3.d0)
      read(5,*) alpha
c      alpha = 45.d0
      beta = 90.d0 - alpha
      alpharad=alpha*pi/180.d0
      a = dsqrt(area*dtan(alpharad))
      b = dsqrt(area/dtan(alpharad))
      c = dsqrt(a*a + b*b)
      write(8,*) a

      write(6,*) '%chk=h6.chk'
      write(6,*) '#rhf/gen scf(tight,novaracc,noincfock,fulldiag)
     & stable(crhf) guess=mix geom=nocrowd test'
      write(6,*) ''
      write(6,*) 'hf stability analysis of h6 triangular lattice, y=',
     & a,'alpha=',alpha
      write(6,*) ''
      write(6,*) '0 1'
      write(6,*) 'H'                                    !1
      write(6,*) 'H ','1',c                             !2
      write(6,*) 'H ','1',c,' 2 ',2.d0*beta             !3
      write(6,*) 'H ','2',c,' 1 ',2.d0*alpha,' 3 ',0.d0 !4
      write(6,*) 'H ','3',c,' 1 ',180.d0,' 2 ',0.d0     !5
      write(6,*) 'H ','4',c,' 2 ',180.d0,' 1 ',0.d0     !6
      write(6,*) ''
      write(6,'(a)') '@H.gbs'
      write(6,*) ''

      end
