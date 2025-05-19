      PROGRAM MODELO
c      	USE MSFLIB
      IMPLICIT REAL *8(a-h,o-z)
            COMPLEX*16 roots(4),a(5)
c      real*8 array(4)

      DIMENSION X(5),V(5),TT(1500000),TET(1500000),
     .                 	FA(1500000),WWZ(1500000)


c	constantes
      pi=4.0d0*datan(1.0d0)      
      pi2=pi+pi
      conv=pi/180.0d0

c	arquivo dados
      OPEN(1,FILE='r1.dat',STATUS='UNKNOWN')
      OPEN(2,FILE='r2.dat',STATUS='UNKNOWN')
      OPEN(3,FILE='r3.dat',STATUS='UNKNOWN')

      eps0=0.96d0
      epsf=1.08d0
      epsi=eps0
      do while(epsi.le.epsf)

      exc=0.0014d0
c      epsilon=1.077d0
      epsilon=epsi
      gama1=3.0d0*exc*dsqrt(epsilon)/8.0d0
      gama2=exc*dsqrt(epsilon**3.0d0)
      delta=epsilon-1.0d0
      bb=0.25d0
      cc=-3.0d0*gama1
      dd=-delta
      ee=gama2


        a(1)=dcmplx(ee)
	a(2)=dcmplx(dd)
	a(3)=dcmplx(cc)
	a(4)=dcmplx(bb)
c	a(5)=dcmplx(bb)

	CALL zroots(a,3,roots,.true.)

        dlim=1.1874d-8
ccccccc ORDENANDO AS RAÍZES

	if( dabs(dimag(roots(1))).gt.dlim ) then
	real1=0.0d0
c        goto 23
	else
	real1=dreal(roots(1))
	endif

c23    continue

	if( dabs(dimag(roots(2))).gt.dlim ) then
c	real2=0.0d0
        goto 24
	else
	real2=dreal(roots(2))
	endif

24    continue

       if( dabs(dimag(roots(3))).gt.dlim ) then
c      real3=0.0d0
        goto 25
	else
	real3=dreal(roots(3))
	endif

25    continue
      write(1,*)epsilon,1.0d0+dsqrt(epsilon)*real1
      write(2,*)epsilon,1.0d0+dsqrt(epsilon)*real2
      write(3,*)epsilon,1.0d0+dsqrt(epsilon)*real3

      deps=(epsf-eps0)/2000.0d0
      epsi=epsi+deps
      
      enddo
                 pause

      END

c--------------------------------------------------------------------
      SUBROUTINE zroots(a,m,roots,polish)

      INTEGER m,MAXM

      REAL EPS

      COMPLEX*16 a(m+1),roots(m)

      LOGICAL polish

      PARAMETER (EPS=1.e-6,MAXM=101)

CU    USES laguer

      INTEGER i,j,jj,its

      COMPLEX*16 ad(MAXM),x,b,c

      do 11 j=1,m+1

        ad(j)=a(j)

11    continue

      do 13 j=m,1,-1

        x=cmplx(0.0,0.0)

        call laguer(ad,j,x,its)

        if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)

        roots(j)=x

        b=ad(j+1)

        do 12 jj=j,1,-1

          c=ad(jj)

          ad(jj)=b

          b=x*b+c

12      continue

13    continue

      if (polish) then

        do 14 j=1,m

          call laguer(a,m,roots(j),its)

14      continue

      endif

      do 16 j=2,m

        x=roots(j)

        do 15 i=j-1,1,-1

          if(real(roots(i)).le.real(x))goto 10

          roots(i+1)=roots(i)

15      continue

        i=0

10      roots(i+1)=x

16    continue

      return

      END

c----------------------------------------------------------------------
      SUBROUTINE laguer(a,m,x,its)

      INTEGER m,its,MAXIT,MR,MT

      REAL EPSS

      COMPLEX*16 a(m+1),x

      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)

      INTEGER iter,j

      REAL abx,abp,abm,err,frac(MR)

      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2

      SAVE frac

      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./

      do 12 iter=1,MAXIT

        its=iter

        b=a(m+1)

        err=abs(b)

        d=cmplx(0.,0.)

        f=cmplx(0.,0.)

        abx=abs(x)

        do 11 j=m,1,-1

          f=x*f+d

          d=x*d+b

          b=x*b+a(j)

          err=abs(b)+abx*err

11      continue

        err=EPSS*err

        if(abs(b).le.err) then

          return

        else

          g=d/b

          g2=g*g

          h=g2-2.*f/b

          sq=sqrt((m-1)*(m*h-g2))

          gp=g+sq

          gm=g-sq

          abp=abs(gp)

          abm=abs(gm)

          if(abp.lt.abm) gp=gm

          if (max(abp,abm).gt.0.) then

            dx=m/gp

          else

            dx=exp(cmplx(log(1.+abx),float(iter)))

          endif

        endif

        x1=x-dx

        if(x.eq.x1)return

        if (mod(iter,MT).ne.0) then

          x=x1

        else

          x=x-dx*frac(iter/MT)

        endif

12    continue

      pause 'too many iterations in laguer'

      return

      END


