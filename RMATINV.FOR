	SUBROUTINE INVER(A,N)
	IMPLICIT NONE
C
	REAL*8 EP
	INTEGER N
	PARAMETER (EP=1.E-6)
C
	INTEGER ME(N),I,J,K,IW
	REAL*8 A(N,N),B(N),C(N),DE,Y,W
C
	DE=1.0
	DO 10 J=1,N
10      ME(J)=J
	DO 20 I=1,N
	Y=0.
	DO 30 J=I,N
	IF(DABS(A(I,J)).LE.DABS(Y)) GOTO 30
	K=J
	Y=A(I,J)
30      CONTINUE
	DE=DE*Y
	IF(DABS(Y).LT.EP) GOTO 32
	Y=1./Y
	DO 40 J=1,N
	C(J)=A(J,K)
	A(J,K)=A(J,I)
	A(J,I)=-C(J)*Y
	B(J)=A(I,J)*Y
40      A(I,J)=A(I,J)*Y
	A(I,I)=Y
	J=ME(I)
	ME(I)=ME(K)
	ME(K)=J
	DO 11 K=1,N
	IF(K.EQ.I) GOTO 11
	DO 12 J=1,N
	IF(J.EQ.I) GOTO 12
	A(K,J)=A(K,J)-B(J)*C(K)
12      CONTINUE
11      CONTINUE
20      CONTINUE
	DO 33 I=1,N
	DO 44 K=1,N
	IF(ME(K).EQ.I) GOTO 55
44      CONTINUE
55      IF(K.EQ.I) GOTO 33
	DO 66 J=1,N
	W=A(I,J)
	A(I,J)=A(K,J)
66      A(K,J)=W
	IW=ME(I)
	ME(I)=ME(K)
	ME(K)=IW            
	DE=-DE
33      CONTINUE
	RETURN
32      DE=0.
	END