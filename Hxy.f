      PROGRAM Cnivel
      IMPLICIT REAL *8(a-h,o-z)
      
c      OPEN(1,FILE='nivel2.grd',STATUS='UNKNOWN')
      
      OPEN(2,FILE='1077.dat',STATUS='UNKNOWN')
            OPEN(3,FILE='xy1077.dat',STATUS='UNKNOWN')
            
      epsilon=1.077d0
      
c	write(1,'(4hDSAA)')
c	write(1,*)'101 101'
c	write(1,*)'-0.001 0.001'
c	write(1,*)'-0.001 +0.001'
c   	write(1,*)'1.E+014 +9.E+020'

   	
      pi=4.0d0*datan(1.0d0)
      conv=pi/180.0d0


c     DIA
      Periodo=1.01445d0
      dn=2.0d0*pi/Periodo
      dn2=dn**2.0d0
      exc=0.00139d0

c     kg, km
      dm=0.73d20
      R=252.1d0
      dmoi=0.335d0
c     dmoi=C/(dm*R**2.0d0)
      C=dmoi*(dm*R**2.0d0)
c c          write(*,*)C
           
c                      write(*,*)
       aa2=256.3d0**2.0d0
       bb2=247.3d0**2.0d0
      C=dm*(aa2+bb2)/5.0d0
c           write(*,*)C
           
c           pause
c       epsilon1=(1.0d0/3.0d0)/(1.0d0+(exc**2.0d0)*(9.0d0/16.0d0)**2.0d0)
c      epsilon= epsilon1

c                write(*,*) epsilon

c      alfa=epsilon*dn*C
c      gama=(3.0d0/16.0d0)*exc*dsqrt(epsilon)
      
c      write(*,*)'gama',gama
      
       delta=epsilon-1.0d0
       write(*,*) epsilon   ,delta

      do y=-0.7d0*dsqrt(epsilon),+0.7d0*dsqrt(epsilon),5.d-4
      do x=-.75d0/dsqrt(epsilon),.65d0/dsqrt(epsilon),5.d-4

      H=       delta*(x**2.0d0+y**2.0d0)/2.0d0
     .         -( (x**2.0d0+y**2.0d0)**2.0d0 )/16.0d0
     .         + 3.0d0*exc*dsqrt(epsilon)*(x**2.0d0+y**2.0d0)*x/8.0d0
     .         - exc*dsqrt(epsilon**3.0d0)*x

      write(2,*)y/dsqrt(epsilon),x*dsqrt(epsilon)+1.0d0,H
             write(3,*)x,y,H
      write(1,*)H

      enddo
      enddo
      
c      pause
      end
