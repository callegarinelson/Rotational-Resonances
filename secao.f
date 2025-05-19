      PROGRAM secaorot
      IMPLICIT REAL *8(a-h,o-z)
      DIMENSION X(3),V(3),TT(1500000),TET(1500000),
     .                 	FA(1500000),WWZ(1500000)
      COMMON/IN/pi,pi2,conv
      COMMON/const/c1,c2,eH,dK,dn,ttp

	is=2
c	write(*,*)'escolha o satelite: 1: Mimas, 2: Enceladus, 3: Hyperion'
c	read(*,*)is
      CALL ENTRE(is)
c	write(*,*)ttp
c	pause

	if(is.eq.1) then
c  	tf=0.94d0/20.0d0
  	tf=ttp/20.0d0	
      tint=140.0d0*ttp
c	write(*,*)tint
	endif

	if(is.eq.2) then
c  	tf=1.37d0/20.0d0
  	tf=ttp/100.0d0
      tint=500.0d0*ttp
	endif

	if(is.eq.3) then
c  	tf=21.3d0/20.0d0
  	tf=ttp/20.0d0	
      tint=100.0d0*ttp
	endif


c      OPEN(2,FILE='f.dat',STATUS='UNKNOWN')
c      OPEN(22,FILE='fwrite.dat',STATUS='UNKNOWN')
   	     
c**** DADOS DE ENTRADA                                           
c	do exc=ei,ef,de
	do tetaloop=-pi,pi,pi/40.0d0
c	write(*,*)tetaloop/conv
	f=0.0d0
	wz=dn
	teta=tetaloop
	
	X(1)=f
	X(2)=wz
	X(3)=teta

      f=f-aint(f/pi2)*pi2
c      if (f.lt.0.0d0) f=f+pi2
c	f=f/conv

      teta=teta-aint(teta/pi2)*pi2
c      if (teta.lt.0.0d0) teta=teta+pi2

c	write(2,*)t,teta,wz,f
c	write(22,*)t,f/conv,teta

	TT(1)=t
	TET(1)=teta
	FA(1)=f
	WWZ(1)=wz


c**** INTEGRACAO NUMERICA E SAIDA DOS DADOS
	jj=2
      do while (t.lt.tint)      
      CALL RA15(TM,X,V,TF,12,XL,3,1)      
C      write(*,*)'call ra15'                 
      t=t+tf   

	f=X(1)
	wz=X(2)
	teta=X(3)

      f=f-aint(f/pi2)*pi2
c      if (f.lt.0.0d0) f=f+pi2
c	f=f/conv
      teta=teta-aint(teta/pi2)*pi2
c      if (teta.lt.0.0d0) teta=teta+pi2


cccccccccc	if(f.gt.270.0d0*conv) tf=1.37d0/1500.0d0	
 
	if(f.ge.340.0d0*conv.AND.f.le.360.0d0*conv) then
	f=-f
	tf=ttp/3600.d0
	endif

	if(f.gt.0.0d0*conv) then
  	tf=ttp/30.0d0	
	endif

	TT(jj)=t
	TET(jj)=teta
	FA(jj)=f
	WWZ(jj)=wz
	jj=jj+1

c	write(2,*)t,teta,wz,f
c	write(22,*)t,f/conv,teta
c	write(*,*)jj
      end do   


c	pause   
	CALL TROCA1(TT,TET,FA,WWZ,jj)
	
	T=0.0D0
	enddo
c	enddo
 
      
      END                               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TROCA1(TT,TETA,F,WZ,jj)
      IMPLICIT REAL *8(a-h,o-z)
	DIMENSION tt(1500000),teta(1500000),wz(1500000),f(1050000)
      COMMON/const/c1,c2,eH,dK,dn,ttp
      COMMON/IN/pi,pi2,conv
c*************************************
C     ABRINDO ARQUIVOS DAS INTEGRAÇÕES
c	write(*,*)jj,'oi'
CCCCCCCCCCCC	OPEN(2,FILE='f.dat',STATUS='UNKNOWN')   
c	write(1,*)t,teta,wz,f

c	atualizando leitor do arquivos da integração
CCCCCCCCCCCCC      rewind(2)

c**********************************
C	ARQUIVO PARA SEÇÃO
      OPEN(20,FILE='secao.dat',STATUS='UNKNOWN')   
c      OPEN(200,FILE='secaoVERIFICA.dat',STATUS='UNKNOWN')   
c**********************************
c	Lendo arquivos para fazer seção: condições iniciais

c	i=1
CCCCCCCCCCCCC	read(2,*,err=2)tt(1),teta(1),wz(1),f(1)

c	Lendo arquivos para fazer seção: loop

CCCCCCCCCCCC1	continue
CCCCCCCCCCCCC	i=i+1

CCCCCCCCCCCCC	read(2,*,err=2)tt(i),teta(i),wz(i),f(i)

cccccccccccccccccccccccccccccccccccccccccccc

	i=2
	do while(i.le.jj-1)

c	seção:  procuro f=0

c	verificando se cruzamento do zero
	prod=f(i)*f(i-1)

c	inicio if produto (prod)	
	if(prod.lt.0.0d0.AND.f(i-1).lt.0.0d0) then

	alfa=( f(i)-f(i-1) )/( tt(i)-tt(i-1) )
	bb1=f(i)-alfa*tt(i)
	ts=-bb1/alfa

c	calculando secao 

	ateta=( teta(i)-teta(i-1) )/( tt(i)-tt(i-1) )
	bteta=teta(i)-ateta*tt(i)
	tetas=ateta*ts+bteta

	awz=( wz(i)-wz(i-1) )/( tt(i)-tt(i-1) )
	bwz=wz(i)-awz*tt(i)
	wzs=awz*ts+bwz

c      tetas=tetas-aint(tetas/pi2)*pi2
c      if (tetas.lt.0.0d0) tetas=tetas+pi2
 
ccccccccc	if (tetas.gt.pi) goto 444
c	tetas=tetas-pi2

c      write(20,*)tetas,wzs/dn,ts
      
c      write(200,*)ts,f(i),f(i-1)
ccccccccccc444	continue


	if (tetas.gt.pi) tetas=tetas-pi2

      write(20,*)tetas/conv,wzs/dn,ts

c      write(200,*)ts,f(i),f(i-1)
ccccccccccc444	continue


c	fim endif do produto
	endif

ccccccccccccccc	go to 1

ccccccccccccccccc2	continue

cccccccccc	fechando arquivos das integrações

cccccccccc	CLOSE(2)

	i=i+1
	enddo


c	pulando uma linha pra distinguir das demais secoes (no caso de varias condicoes iniciais)
	write(20,*)

	RETURN
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                 
      SUBROUTINE FORCE (TM,X,V,F1)
      IMPLICIT REAL *8(a-h,o-z)
      DIMENSION X(3),V(3),F1(3)
      COMMON/const/c1,c2,eH,dK,dn,ttp
c	aalfa=(bb-aa)/cc
c	aalfa=0.26d0
c	dK=-1.5d0*G*dm0*aalfa
c	c1=aH*(1.0d0-eH**2.0d0)
c	c2=dsqrt( amiH*aH*(1.0d0-eH**2.0d0) )
c	write(*,*)c1,c2,dk

	f=X(1)
	wz=X(2)
	teta=X(3)
	r=c1/(1.0d0+eH*dcos(f))

c	df/dt
	F1(1)=c2/(r*r)
c	dwz/dt
	F1(2)=dK*dsin( 2.0d0*(f-teta) )/(r**3.0d0)
c	dteta/dt
	F1(3)=wz
      
      return
      end      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      SUBROUTINE ENTRE(is)
      implicit real *8(a-h,o-z)
      COMMON/const/c1,c2,eH,dK,dn,ttp
      COMMON/IN/pi,pi2,conv
      pi=4.0d0*datan(1.0d0)      
      pi2=pi+pi
      conv=pi/180.0d0
 
      OPEN(2588,FILE='Pan.dat',STATUS='UNKNOWN')

	dm0=1.0d0  
	grav=4.0d0*pi*pi
ccccc      G=39.476926421373020d0 

c	HORIZONS DATA SYSTEM PARA DATA 2007-JAN 01- 00:00
      ua=1.49597870691d11                                                   metros

      requat=6.0268d7
      esc=ua/requat

c	em	quilos
c	dmsatt=1.9891d30/5.68462313752d26
c      write(*,*)dmsat
c	=3499.08859722190

c	JACOBSON 2006
	dmsun=3497.9018d0*37940585.2d0/6.6742d-23
c	write(*,*)dmsun

c	Esta quantidade e' 1/M
	dmsatt=dmsun/(37931207.7d0/6.6742d-23)
c      write(*,*)dmsatt,i
c	pause

      grav=grav*esc*esc*esc
	G=grav/(dmsatt*365.25d0*365.25d0)


ccccccccccccccc Enceladus
	if(is.eq.2) then
c	aalfa=(bb-aa)/cc
c	aalfa=0.26d0
	aalfa=0.03853d0
	
c      aalfa= ( (aaa**2.0d0+ccc**2.0d0)-(bbb^2+ccc^2) )/(aaa^2+bbb^2)
      
      aaa=1.94d0
      bbb=1.29d0
      aalfa= (aaa**2.0d0-bbb**2.0d0)/(aaa**2.0d0+bbb**2.0d0)


	dK=1.5d0*G*dm0*aalfa

      dmH=1.901d-7
      dmH=0.0d0
      
c10.805d22/5.68462313752d29
      amiH=G*(dm0+dmH)
      dmimiH=(dm0+dmH)/(dm0*dmH)
      aH=0.001593665160605842d0*ua/requat

      aH=133583.0d0/60268.0d0
      aH=1.342636274754209d+05/60268.0d0
      
         aH= 1.947075740548770E+05 /60268.0d0

      
	ttp=2.0d0*pi*aH**(3.0d0/2.0d0)/dsqrt(amiH)
c	write(*,*)aH,ttp
c	pause
	dn=dsqrt( amiH/(aH**3.0d0) )

      eH=1.391889851375553d-03

        domega2=3.0d0*dn*dn*aalfa

        dkk11=2.0d0*domega2*eH/(dn*dn-domega2)

        shift=1.0d0-dkk11

          write(2588,*)'(B-A)/C, 1-k1:1',aalfa,shift
          write(2588,*)'k1:1',dkk11
          write(2588,*)'w02, denominador',domega2,dn*dn-domega2
          write(2588,*)'T e Tw0',ttp,2.0d0*pi/(dsqrt(domega2))

	c1=aH*(1.0d0-eH**2.0d0)
	c2=dsqrt( amiH*aH*(1.0d0-eH**2.0d0) )

c	write(*,*)c1,c2,dk
	endif




cccccccccccccccccccccccccccccccc

c	JACOBSON 2005
c      dj2= 0.0162906d0
c	dj4=-0.000936d0
C	J2, J4
C	JACOBSON 2006
      dj2= 0.01629071d0
	dj4=-0.0009358d0

c     grav=G      
c     ua=1.4959787d11  metros        
 
c:    GRAV acima em UA**3/(Massa Solar*Ano**2)
c:    UA acima em metros  e Raio equatorial=requat sera em METROS (aqui)
 
      
      return
	end       

C----------------------------------------------------------------------

      SUBROUTINE RA15(TM,X,V,TF,LL,XL,NV,NCLASS)
C  Integrador RADAU de Everhart - Version de orden 15
C  Integra ecuaciones de segundo orden:
C    y' = F(y,t)    es NCLASS = 1
C    y" = F(y',y,t) es NCLASS = 2  
C    y" = F(y,t)    es NCLASS = -2
C  TF es t(final) - t(initial). Debe ser negativo para 
C    integracion hacia atras
C  NV es el numero de ecuaciones diferenciales simultaneas
C  LL controla el tamanio de paso, monitoreando el error en los
C    terminos en base a la tolerancia SS = 10**(-LL). Un valor
C    tipico para LL esta entre 6 y 12
C  XL es el tamanio de paso constante si LL < 0. Si no, es el
C    tamanio de paso inicial
C  TM, X y V son tiempo, coordenadas y velocidades iniciales en
C    la entrada. En la salida son devueltas como tiempo, 
C    coordenadas y velocidades finales.
      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(60),V(60),F1(60),FJ(60),C(21),D(21),R(21),
     >     Y(60),Z(60),B(7,60),G(7,60),E(7,60),BD(7,60),H(8),
     >     W(7),U(7),NW(8)
      LOGICAL NSF,NPER,NPQ,NCL,NES
      DATA NW/0,0,1,3,6,10,15,21/
      DATA ZERO, HALF, ONE, SR/0.0D0, 0.5D0, 1.0D0, 1.4D0/
      DATA H/         0.D0, .05626256053692215D0, .18024069173689236D0,
     >.35262471711316964D0, .54715362633055538D0, .73421017721541053D0,
     >.88532094683909577D0, .97752061356128750D0/
C  Estos valores de H son los espaciadores de Gauss - Radau,
C  escalados en el intervalo [0,1]

      NPER=.FALSE.
      NSF=.FALSE.
      NCL=NCLASS.EQ.1
      NPQ=NCLASS.LT.2
      NES=LL.LT.0
C  NCLASS =  1   ==>   NCL = .TRUE.    NPQ = .TRUE.
C  NCLASS = -2   ==>   NCL = .FALSE.   NPQ = .TRUE.
C  NCLASS =  2   ==>   NCL = .FALSE.   NPQ = .FALSE.
C  NPER es .TRUE. solo en la ultima secuencia de integracion.
C  NSF es .FALSE. solo en la secuencia inicial
C  NES es .TRUE. solo si LL es negativo. Entonces el tamanio de
C  paso es XL.
      
      IF (TF.LT.ZERO) THEN
       DIR=-ONE
      ELSE
       DIR=ONE
      END IF
      
      XL=DIR*DABS(XL)
      PW=1./9.
      
      DO N=2,8
       WW=N+N*N
       IF (NCL) WW=N
       W(N-1)=ONE/WW
       WW=N
       U(N-1)=ONE/WW
      END DO
      
      DO K=1,NV
       IF (NCL) V(K)=ZERO
       DO L=1,7
	BD(L,K)=ZERO
	B(L,K)=ZERO
       END DO
      END DO
      
      W1=HALF
      IF (NCL) W1=ONE
      C(1)=-H(2)
      D(1)=H(2)
      R(1)=ONE/(H(3)-H(2))
      LA=1
      LC=1
      
      DO K=3,7
       LB=LA
       LA=LC+1
       LC=NW(K+1)
       C(LA)=-H(K)*C(LB)
       C(LC)=C(LA-1)-H(K)
       D(LA)=H(2)*D(LB)
       D(LC)=-C(LC)
       R(LA)=ONE/(H(K+1)-H(2))
       R(LC)=ONE/(H(K+1)-H(K))
       
       IF (K.NE.3) THEN
	DO L=4,K
	 LD=LA+L-3
	 LE=LB+L-4
	 C(LD)=C(LE)-H(K)*C(LE+1)
	 D(LD)=D(LE)+H(L-1)*D(LE+1)
	 R(LD)=ONE/(H(K+1)-H(L-1))
	END DO
       END IF
      
      END DO
      
      SS=10.**(-LL)
C  Las instrucciones anteriores son calculadas solo una vez
C  en una integracion para inicializar las constantes

C  A continuacion se inicializa el tamanio de paso inicial TP
      IF (NES) THEN
       TP=XL
      ELSE IF (XL.NE.ZERO) THEN
       TP=XL
      ELSE
       TP=0.1D0*DIR
      END IF
      
      IF (TP/TF.GT.HALF) TP=HALF*TF
      NCOUNT=0

C  La linea 1000 es el comienzo de la primera secuencia.
C  NS es el numero de secuencia. 
C  NF es el numero de llamados a la subrrutina FORCE. 
C  NI es el numero de iteraciones en cada secuencia
1000  NS=0
      NF=0
      NI=6
      TM=ZERO
      CALL FORCE (TM, X, V, F1)
      NF=NF+1

C  La linea 2000 es el comienzo de cada secuencia despues de
C  la primera
2000  DO K=1,NV
       G(1,K)=B(1,K)+D( 1)*B(2,K)+D(2)*B(3,K)+
     >   D(4)*B(4,K)+D( 7)*B(5,K)+D(11)*B(6,K)+D(16)*B(7,K)
       G(2,K)=             B(2,K)+D(3)*B(3,K)+
     >   D(5)*B(4,K)+D( 8)*B(5,K)+D(12)*B(6,K)+D(17)*B(7,K)
       G(3,K)=                         B(3,K)+
     >   D(6)*B(4,K)+D( 9)*B(5,K)+D(13)*B(6,K)+D(18)*B(7,K)
       G(4,K)=B(4,K)+D(10)*B(5,K)+D(14)*B(6,K)+D(19)*B(7,K)
       G(5,K)=             B(5,K)+D(15)*B(6,K)+D(20)*B(7,K)
       G(6,K)=                          B(6,K)+D(21)*B(7,K)
       G(7,K)=                                       B(7,K)
      END DO
      
      T=TP
      T2=T*T
      IF (NCL) T2=T
      TVAL=DABS(T)

C  Comienzo de las iteraciones para cada secuencia. Realiza 6
C  iteraciones en la primera secuencia (NI=6)
      DO M=1,NI
       DO J=2,8
	JD=J-1
	JDM=J-2
	S=H(J)
	Q=S
	IF (NCL) Q=ONE
	
C  Aqui calcula los predictores de posicion y velocidad para
C  cada subsecuencia        
	DO K=1,NV
	 A=W(3)*B(3,K)+S*(W(4)*B(4,K)+S*(W(5)*B(5,K)+S*(W(6)*
     >     B(6,K)+S*W(7)*B(7,K))))
	 Y(K)=X(K)+Q*(T*V(K)+T2*S*(F1(K)*W1+S*(W(1)*B(1,K)+
     >        S*(W(2)*B(2,K)+S*A))))
	 
	 IF (.NOT.NPQ) THEN
C         Esta parte solo es calculada cuando NCLASS = 2          
	  A=U(3)*B(3,K)+S*(U(4)*B(4,K)+S*(U(5)*B(5,K)+
     >      S*(U(6)*B(6,K)+S*U(7)*B(7,K))))
	  Z(K)=V(K)+S*T*(F1(K)+S*(U(1)*B(1,K)+S*(U(2)*B(2,K)+
     >      S*A)))
	 END IF
	
	END DO
	
C  Calcula la fuerza FJ al final de cada subsecuencia
	CALL FORCE(TM+S*T, Y, Z, FJ)
	NF=NF+1
	
C  Calcula los nuevos valores de G y de B para la fuerza FJ        
	DO K=1,NV
	 TEMP=G(JD,K)
	 GK=(FJ(K)-F1(K))/S
	 
	 GO TO (102,102,103,104,105,106,107,108),J
 102        G(1,K)=GK
	    GO TO 100
 103        G(2,K)=(GK-G(1,K))*R(1)
	    GO TO 100
 104        G(3,K)=((GK-G(1,K))*R(2)-G(2,K))*R(3)
	    GO TO 100
 105        G(4,K)=(((GK-G(1,K))*R(4)-G(2,K))*R(5)-G(3,K))*
     >             R(6)
	    GO TO 100
 106        G(5,K)=((((GK-G(1,K))*R(7)-G(2,K))*R(8)-G(3,K))*
     >             R(9)-G(4,K))*R(10)
	    GO TO 100
 107        G(6,K)=(((((GK-G(1,K))*R(11)-G(2,K))*R(12)-
     >             G(3,K))*R(13)-G(4,K))*R(14)-G(5,K))*R(15)
	    GO TO 100
 108        G(7,K)=((((((GK-G(1,K))*R(16)-G(2,K))*R(17)-
     >             G(3,K))*R(18)-G(4,K))*R(19)-G(5,K))*R(20)-
     >             G(6,K))*R(21)
 100        CONTINUE
	 
	 TEMP=G(JD,K)-TEMP
	 B(JD,K)=B(JD,K)+TEMP
	 
	 GO TO (200,200,203,204,205,206,207,208),J
 203        B(1,K)=B(1,K)+C(1)*TEMP
	    GO TO 200
 204        B(1,K)=B(1,K)+C(2)*TEMP
	    B(2,K)=B(2,K)+C(3)*TEMP
	    GO TO 200
 205        B(1,K)=B(1,K)+C(4)*TEMP
	    B(2,K)=B(2,K)+C(5)*TEMP
	    B(3,K)=B(3,K)+C(6)*TEMP
	    GO TO 200
 206        B(1,K)=B(1,K)+C(7)*TEMP
	    B(2,K)=B(2,K)+C(8)*TEMP
	    B(3,K)=B(3,K)+C(9)*TEMP
	    B(4,K)=B(4,K)+C(10)*TEMP
	    GO TO 200
 207        B(1,K)=B(1,K)+C(11)*TEMP
	    B(2,K)=B(2,K)+C(12)*TEMP
	    B(3,K)=B(3,K)+C(13)*TEMP
	    B(4,K)=B(4,K)+C(14)*TEMP
	    B(5,K)=B(5,K)+C(15)*TEMP
	    GO TO 200
 208        B(1,K)=B(1,K)+C(16)*TEMP
	    B(2,K)=B(2,K)+C(17)*TEMP
	    B(3,K)=B(3,K)+C(18)*TEMP
	    B(4,K)=B(4,K)+C(19)*TEMP
	    B(5,K)=B(5,K)+C(20)*TEMP
	    B(6,K)=B(6,K)+C(21)*TEMP
 200        CONTINUE
	
	END DO
       END DO
       
C  Final de la secuencia. Calculo del control (HV) del tamanio
C  de paso       
       IF (.NOT.NES.OR.M.GE.NI) THEN
	HV=ZERO
	DO K=1,NV
	 HV=DMAX1(HV,DABS(B(7,K)))
	END DO
	TVAL1=TVAL*TVAL*TVAL*TVAL*TVAL*TVAL*TVAL
	HV=HV*W(7)/TVAL1
       END IF
      
      END DO
      
C  Calcula el nuevo tamanio de paso (TP) y lo compara con HV.
C  Si es menor que lo previsto, reinicia la secuencia anterior
C  usando un paso 0.8 veces menor
      IF (.NOT.NSF) THEN
       
       IF (NES) THEN
	TP=XL
       ELSE
	TP=(SS**PW)/(HV**PW)*DIR
	IF (TP/T.LE.ONE) THEN
	 TP=.8D0*TP
	 NCOUNT=NCOUNT+1
	 IF (NCOUNT.GT.10) RETURN
	 GO TO 1000
	END IF
       END IF
       
       NSF=.TRUE.
      END IF

C  Coordenadas y velocidades y tiempo al final de la secuencia
      DO K=1,NV
       X(K)=X(K)+V(K)*T+T2*(F1(K)*W1+B(1,K)*W(1)+B(2,K)*W(2)+
     >      B(3,K)*W(3)+B(4,K)*W(4)+B(5,K)*W(5)+B(6,K)*W(6)+
     >      B(7,K)*W(7))
       IF (.NOT.NCL) THEN
	V(K)=V(K)+T*(F1(K)+B(1,K)*U(1)+B(2,K)*U(2)+B(3,K)*U(3)
     >       +B(4,K)*U(4)+B(5,K)*U(5)+B(6,K)*U(6)+B(7,K)*U(7))
       END IF
      END DO
      
      TM=TM+T
      NS=NS+1
      
C  Si termino retorna. Caso contrario controla el tamanio de
C  la siguiente secuencia y la ajusta para cubrir exactamente
C  el intervalo de total de integracion (TF)
      IF (NPER) RETURN
      CALL FORCE (TM, X, V, F1)
      NF=NF+1
      
      IF (NES) THEN
       TP=XL
      ELSE
       TP=DIR*(SS**PW)/(HV**PW)
       IF (TP/T.GT.SR) TP=T*SR
      END IF
      
      IF (DIR*(TM+TP).GE.DIR*TF-1.D-8) THEN
       TP=TF-TM
       NPER=.TRUE.
      END IF
  
C  Ahora predice los nuevos valores de B para utilizarlos en
C  la proxima secuencia      
      Q=TP/T
      
      DO K=1,NV
       
       IF (NS.NE.1) THEN
	DO J=1,7
	 BD(J,K)=B(J,K)-E(J,K)  
	END DO
       END IF
       
       E(1,K)=     Q*(B(1,K)+ 2.D0*B(2,K)+ 3.D0*B(3,K)+
     >           4.D0*B(4,K)+ 5.D0*B(5,K)+ 6.D0*B(6,K)+ 
     >           7.D0*B(7,K))
       E(2,K)=               Q**2*(B(2,K)+ 3.D0*B(3,K)+
     >           6.D0*B(4,K)+10.D0*B(5,K)+15.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(3,K)=                            Q**3*(B(3,K)+
     >           4.D0*B(4,K)+10.D0*B(5,K)+20.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(4,K)=  Q**4*(B(4,K)+ 5.D0*B(5,K)+15.D0*B(6,K)+
     >          35.D0*B(7,K))
       E(5,K)=               Q**5*(B(5,K)+ 6.D0*B(6,K)+
     >          21.D0*B(7,K))
       E(6,K)=                            Q**6*(B(6,K)+
     >           7.D0*B(7,K))
       E(7,K)=   Q**7*B(7,K)
       
       DO L=1,7
	B(L,K)=E(L,K)+BD(L,K)
       END DO
      
      END DO

C  A partir de la segunda secuencia, solo realiza dos 
C  iteraciones por cada secuencia
      NI=2
      GO TO 2000

      END
	


c13    FORMAT(1F8.1,1X,1F9.6,1X,1F10.6,1X,1F8.6,1X,1F12.6,
c     /1X,1F12.6,1X,1F12.6)
C14    FORMAT(1F8.1,1X,1F9.6,1X,1F10.6,1X,1F8.6,1X,1F12.6,
C     /1X,1F12.6,1X,1F10.6,1X,i1)
c     ,1X,1F10.6,1X,i1)
c      write(*,13)t,au,am0u/conv,eu,wu/conv,wu/conv,diu/conv
c      ,omu/conv,lei       
c      write(1,13)t,au,am0u/conv,eu,wu/conv,wu/conv,diu/conv
c      ,omu/conv,lei 
c      write(1,*)
c      write(*,12)t,au,am0u/conv,eu,wu/conv,diu/conv,omu/conv,lei
c      write(1,12)t,au,am0u/conv,eu,wu/conv,diu/conv,omu/conv,lei 
c      write(*,*)     
c      write(*,12)t,an,am0n/conv,en,wn/conv,din/conv,omn/conv,lei        
c      write(2,12)t,an,am0n/conv,en,wn/conv,din/conv,omn/conv,lei  
C      write(*,13)t,an,am0n/conv,en,wn/conv,wn/conv,din/conv
c      ,omn/conv,lei        
C      write(2,13)t,an,am0n/conv,en,wn/conv,wn/conv,din/conv
c      ,omn/conv,lei  
c      write(2,*)           
c      write(*,*)   

C      write(*,13)t,au,amu/conv,eu,wu/conv,wwu/conv,diu/conv
C      ,omu/conv,lei       
C      write(1,13)t,au,amu/conv,eu,wu/conv,wwu/conv,diu/conv
C      ,omu/conv,lei
c      write(*,12)t,au,amu/conv,eu,wu/conv,diu/conv,omu/conv,lei       
c      write(1,12)t,au,amu/conv,eu,wu/conv,diu/conv,omu/conv,lei 
c      write(1,*)t,wtifu/conv,fu/conv,xura(4),rcosu,alatu/conv      
c      write(1,*)t,wtifu/conv,xura(5),xura(4),xura(5)/xura(4),
c     ?fu/conv,rcosfu,rsinfu 
c      write(1,*)             
c      write(*,*)     
C      write(*,13)t,an,amn/conv,en,wn/conv,wwn/conv,din/conv
c      ,omn/conv,lei        
C      write(2,13)t,an,amn/conv,en,wn/conv,wwn/conv,din/conv
c      ,omn/conv,lei
c      write(*,12)t,an,amn/conv,en,wn/conv,din/conv,omn/conv,lei        
c      write(2,12)t,an,amn/conv,en,wn/conv,din/conv,omn/conv,lei
c      write(2,*)t,wtifn/conv,fn/conv,xnet(4),rcosn,alatn/conv       
c      write(2,*)t,wtifn/conv,xnet(5),xnet(4),xnet(5)/xnet(4),
c     ?fn/conv,rcosfn,rsinfn 
c      write(2,*)          
c      write(*,*)




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

c      SUBROUTINE ENERGIA(H,X)
c      IMPLICIT REAL *8(a-h,o-z)
c      DIMENSION X(12)
c      COMMON/MASSAS/dm0,G,dj2,dj4


c      xr=x(1)-x(4)
c      yr=x(2)-x(5)
c      zr=x(3)-x(6) 
c      d=sqrt(xr**2+yr**2+zr**2)    
c      ru=sqrt(x(1)**2+x(2)**2+x(3)**2)
c      rn=sqrt(x(4)**2+x(5)**2+x(6)**2)                                         

      
c      U0=-cmu/ru-cmn/rn

c      pura2=x(7)**2+x(8)**2+x(9)**2
c      pnet2=x(10)**2+x(11)**2+x(12)**2            
c      T0=0.5d0*dmiu*pura2+0.5d0*dmin*pnet2

c      H0=U0+T0

c      U1=-G*dkm/d        
      
c      T1=x(7)*x(10)+x(8)*x(11)+x(9)*x(12)
      
c      H1=U1+T1      


c      H=H0+H1 
      
c      return
c      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C	SUBROUTINE SISTEMA(b,ami,N)
C	implicit real *8(a-h,o-z)
C      DIMENSION a(24,24),b(24,1),ami(8)

c	write(*,*)'ordem da matriz, n='
c	read(*,*)n
	
C	m=1

C	k=1
C	do 2 i=1,n
C	do 1 j=1,n

C	if(i.eq.j) then
C	a(i,j)=ami(k)
C	k=k+1
C	else
C	a(i,j)=1.0d0
C	endif

C1     continue
C2     continue

C	call  gaussj(a,n,b,m)

c     Sistema qualquer
     
c     do 4 i=1,n
c	do 3 j=1,n
c	write(*,*)'vetor a(i,j)',i,j,a(i,j)
c	pause
c	read(*,*)a(i,j)
c3	continue
c	write(*,*)'vetor b(i,1)',i,1
c	read(*,*)b(i,1)
c4	continue
c     conferir valores de a(i,j)
c      do 4 l=1,n
c	do 3 ll=1,n
c	write(*,*)'vetor a(l,ll)',l,ll,a(l,ll)
c	pause
c3	continue
c4     continue
	
C	SAIDA
	     
c      do 5 i=1,n

c     matriz inversa
      
c	do 4 j=1,n

c	write(*,*)'i,j,a(i,j)',i,j,a(i,j)	
c	pause
c4	continue
c	write(*,*)'i,b(i,1)',i,b(i,1)
c	pause

c5	continue
C      return
C	END

cccccccccccccccccccccccccccccccccccccccccccccccc

c     reducao de wu e wn para 0-2pi
c      if (lei.eq.1) go to 3
c     se lei=1 => singularidade em eu, en (ver xyzorb), 
c     e wu, wn recebem 9.d10 ficticios, e portanto nao existe reducao
c        inwu=wu/pi2
c        wu=wu-inwu*pi2
c        if(wu.lt.0.0d0) wu=wu+pi2                   
c        inwn=wn/pi2
c        wn=wn-inwn*pi2
c        if(wn.lt.0.0d0) wn=wn+pi2                 
c3     continue          

cccccccccccccccccccccccccccccccccccccccccccccccccc

c      CALL ENERGIA(H,X)       
c      H=H*100000.0D0
c7     format(1f8.1,1x,1f20.15)      
c      write(3,*)t,H          
c      write(*,*)t,H
c      write(*,*)

C      dhu=eu*dsin(wu)
C      dku=eu*dcos(wu)
C      write(*,*)dku,dhu
C      write(4,*)dku,dhu      
C      dhn=en*dsin(wn)
C      dkn=en*dcos(wn)
C      write(*,*)dkn,dhn
C      write(5,*)dkn,dhn      
C      write(*,*) 



C      dhu=eu*dsin(wu)
C      dku=eu*dcos(wu)
C      write(*,*)dku,dhu
C      write(4,*)dku,dhu      
C      dhn=en*dsin(wn)
C      dkn=en*dcos(wn)
C      write(*,*)dkn,dhn
C      write(5,*)dkn,dhn      
C      write(*,*)


cccccccccccccccccccccccccccccccccccccccccccccccc	

c     Conferir os valores de xyz     
c      do 77 linha=1,N
c      write(*,*)'f(i)',f(linha),linha
c      do 66 lcol=1,6
c      write(*,*)'xpla(i,j),linha,coluna',xpla(linha,lcol),linha,lcol
c66    continue
c      write(*,*)                            
c      pause
c77    continue
c      pause


ccccccccccccccccccccccccccccccccccccccccccccccccc
C      CHARACTER NOME1*20  
C      CHARACTER NOME2*20      
C      CHARACTER NOME3*20      
c      CHARACTER NOME4*20  
c      CHARACTER NOME5*20      
C      WRITE(*,*) 'Arq. Plan.:'
C      READ(*,48)NOME1
C      WRITE(*,*) 'Arq. U:'
C      READ(*,48)NOME2      
C      WRITE(*,*) 'Arq. N:'
C      READ(*,48)NOME3     
C48     FORMAT(A)   
cccccccccc

c      do 7 i=1,5	
c      write(*,12)t,a(i),am(i)/conv,e(i),w(i)/conv,di(i)/conv,
c     |om(i)/conv,lei
c      write(1,12)t,a(i),am(i)/conv,e(i),w(i)/conv,di(i)/conv,
c     |om(i)/conv,lei 
c7     continue   
c      write(*,*)
c      write(1,*)
cccccccccccccccccccccccccccc


  
c: Mercurio
c      requat=2.4397d6
c      dj2=0.027d-3
C      write(*,*)'entre dj2 de Mercurio,sabendo que o de Venus=0.027d-3'
C      read(*,*) dj2
c      dm(1)=1.d0/6.023600d6  
c      ami(1)=G*(dm0+dm(1))
c      dmimi(1)=(dm0+dm(1))/(dm0*dm(1))

c      a(1)=0.387098d0
c      e(1)=0.205633d0      
c      w(1)=8.d0*conv
c      om(1)=0.d0*conv 
c      am(1)=1.d0*conv
c      di(1)=0.0d0*conv
c: atual    epsil=0.0d0 (pag. E87-Nautical Almanac 96)        
C        write(*,*)'entre epsil de Mercurio , sendo o atual=0.d0'
C        read(*,*)epsil

  
c: Venus
c      requat=6.0518d6
c:c:     dj2=0.027d-3     c: c: c:
C      write(*,*)' entre J2 de Venus c/ rotacao (hoje 0.027d-3) '
c      , sendo atualmente=0.027d-3'
C      read(*,*)dj2
c      dm(2)=1.d0/408523.5d0
c      ami(2)=G*(dm0+dm(2))
c      dmimi(2)=(dm0+dm(2))/(dm0*dm(2))      

c      a(2)=0.7233d0
c      e(2)=0.0067d0      
c      w(2)=7.d0*conv
c      om(3)=0.d0*conv 
c      am(3)=2.d0*conv
c      di(2)=0.0d0*conv
c: atual    epsil=177.36d0         
C      write(*,*)'entre epsil de Venus , sendo o atual=177.36'
C      read(*,*)epsil


c: Terra                                                                
c      requat=0.637814d7                                                 Naut
c      dj2=1.08263d-3                                                    Naut
                                                    
c: P/Terra:tomarei: do Connaissaince:  semi-eixo, massa,J2,Requat,epsil
c: para aj vamos tomar aj=1.000001d0 (Connaissance 92) em UA
c: esol: tirado Nautical Almanac 96 (~=6 casas) 
c      dm(3)=1.d0/332946.d0
c      ami(3)=G*(dm0+dm(3))
c      dmimi(3)=(dm0+dm(3))/(dm0*dm(3))                                                                         Conn
C      epsil=23.43929115d0                                                Conn
C      write(*,*)'entre epsil da Terra, sendo o atual 23.44 graus'
C      read(*,*) epsil
c      a(3)=1.000001d0     
c      e(3)=0.01672d0                                                      Naut
c      w(3)=6.d0*conv
c      om(3)=0.d0*conv        
c      am(3)=3.d0*conv
c      di(3)=0.0d0*conv


c: Marte
c      requat=3397.d3                                                    Naut
c      dj2=0.001964d0                                                    Naut
c      dm(4)=1.d0/3098710.d0
c      ami(4)=G*(dm0+dm(4))                                               Conn
c      dmimi(4)=(dm0+dm(4))/(dm0*dm(4))

                                                                      
c      a(4)=1.5236d0                                                       Naut
c      e(4)=0.0933d0                                                       Naut
c      w(4)=5.d0*conv
c      om(4)=0.d0*conv 
c      am(4)=4.d0*conv
c      di(4)=0.0d0*conv
c: atual        epsil=25.19d0                                           Naut
C       write(*,*)'entre epsil de Marte, sendo  o atual= 25.19'
C       read (*,*) epsil  
c: se tomo epsil=Incl=29.2 e hnod=w=0 , ha um grande salto em excent.    
c: abaixo dados usados por Adrian
c:    requat=3.3972d6
c:    dj2=0.001964d0
c:    ej=0.0934183
c:    dmj=1.d0/3098710.d0
c:    aj=1.52236636d0
c:    esc=ua/requat
c:    grav=grav*esc*esc*esc
c:    write(*,*)' entre epsil de Marte'
c:    read(*,*) epsil
