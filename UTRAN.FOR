	SUBROUTINE UTRAN(ZT,NS)        
C-----------------------------------------------------------------------
C Subroutine to generate the up-going waves and down-going waves.
C   UP: U+(s+1)-    
C   UD: U_(s)+        
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	REAL*8 ZT,HU,HD
	INTEGER NS
C
	INTEGER I,J,J1,J2
	COMPLEX*16 X1(N,N),X2(N,N),X3(N,N),X4(N,N)
	COMPLEX*16 E1(N),E2(N),V1(N),V2(N),COEF
C
	DO 2 I=1,N
	  V1(I)=(1.0D0,0.0D0)/V(NS,I)
	  DO 2 J=1,N
	  X1(I,J)=0.0D0
	  X2(I,J)=0.0D0
	  IF(NS.NE.NBOU) X1(I,J)=QU(NS+1,I,J)
	  IF(NS.NE.0) X2(I,J)=QD(NS,I,J)
	  X3(I,J)=C(NS,J,I)
2       CONTINUE          
C	
	IF(NS.NE.NBOU) THEN
	HU=HOR(NS+1)-ZT
	DO 5 I=1,N
5     E1(I)=CDEXP(-V(NS,I)*HU)
	CALL AMULT(E1,X1,X4,0,N)
	CALL AMULT(E1,X4,X1,1,N)
	END IF
C----- Where X1 is Q+(s+1)-
C
	IF(NS.NE.0) THEN
	HD=ZT-HOR(NS)
	DO 7 I=1,N
7     E2(I)=CDEXP(-V(NS,I)*HD)
	CALL AMULT(E2,X2,X4,0,N)
	CALL AMULT(E2,X4,X2,1,N)
	END IF
C----- Where X2 is Q-(s)+
C	
	CALL GCAL(RADIUST,G)
	DO 10 I=1,N
10    V2(I)=CMPLX(G(I))
C
	CALL AMULT(V1,X3,X4,0,N)
	CALL BMULT(X4,V2,V1,N)
C	
	COEF=-(0.0D0,0.5D0)*OMEGA*U0*CURENT
	DO 12 I=1,N
12    V1(I)=COEF*V1(I)
C----- Where V1 is bs.
C
	IF(NS.EQ.0) THEN
	  DO 15 I=1,N
	  UP(I)=V1(I)
15      X1(I,I)=(1.0D0,0.0D0)+X1(I,I)        
	  CALL BMULT(X1,V1,UD,N)
	ELSEIF(NS.EQ.NBOU) THEN
	  DO 17 I=1,N
	  UD(I)=V1(I)
17      X2(I,I)=(1.0D0,0.0D0)+X2(I,I)        
	  CALL BMULT(X2,V1,UP,N)
	ELSE
	  CALL MULT(X1,X2,X3,N)
	  CALL MULT(X2,X1,X4,N)
C	  
	  DO 20 I=1,N
	  DO 20 J=1,N
	  IF(I.EQ.J) THEN
	  X1(I,J)=(1.0D0,0.0D0)+X1(I,J)
	  X2(I,J)=(1.0D0,0.0D0)+X2(I,J)
	  X3(I,J)=(1.0D0,0.0D0)-X3(I,J)
	  X4(I,J)=(1.0D0,0.0D0)-X4(I,J)
	  ELSE
	  X3(I,J)=-X3(I,J)
	  X4(I,J)=-X4(I,J)
	  END IF
20      CONTINUE
C
	  CALL INV(X3,N)
	  CALL INV(X4,N)
	  CALL BMULT(X1,V1,V2,N)
	  CALL BMULT(X3,V2,UD,N)
	  CALL BMULT(X2,V1,V2,N)
	  CALL BMULT(X4,V2,UP,N)
C
	END IF
C
C----- For up-going wave
C
	IF(NS.NE.NBOU) THEN
	DO 22 I=1,N
22    V1(I)=E1(I)*UP(I)
	DO 25 I=NS+1,NBOU
	  DO 27 J1=1,N
	  DO 27 J2=1,N
27      X1(J1,J2)=SU(I,J1,J2)
	  CALL BMULT(X1,V1,V2,N)
	  DO 30 J=1,N
	  IF(I.NE.NBOU) V1(J)=E(I,J)*V2(J)
30      U(I,J)=V2(J)
25      CONTINUE
	END IF
C
C----- For down-going wave
C
	IF(NS.NE.0) THEN
	DO 32 I=1,N
32    V1(I)=E2(I)*UD(I)
	DO 35 I=NS,1,-1
	  DO 37 J1=1,N
	  DO 37 J2=1,N
37      X1(J1,J2)=SD(I,J1,J2)
	  CALL BMULT(X1,V1,V2,N)
	  DO 40 J=1,N
	  IF(I.NE.1) V1(J)=E(I-1,J)*V2(J)
40      U(I,J)=V2(J)
35      CONTINUE
	END IF
C
	END