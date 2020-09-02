	SUBROUTINE QSTRAN       
C-----------------------------------------------------------------------
C Subroutine to evaluate the generalized reflection matrices and 
C generlized transmission matrices at the formation boundries.
C   P       : Pm,m+1
C   E       : CEXP(-Vm*hm)
C   Q+_(m)+_: generalized reflection matrix
C   S(m)+_  : generalized transmission matrix
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	COMPLEX*16 E1(N),V1(N),V2(N)
	COMPLEX*16 X1(N,N),X2(N,N),X3(N,N),P0(N,N)
C
	INTEGER NBED,I,J,J1,J2,IFLAG
C
	NBED=NBOU-1
	DO 2 I=1,NBED
	DO 2 J=1,N
2     E(I,J)=CDEXP(-V(I,J)*(HOR(I+1)-HOR(I)))
	
	DO 5 I=1,NBOU
	  DO 7 J1=1,N
	  DO 7 J2=1,N
	  X1(J1,J2)=C(I-1,J2,J1)
7       X2(J1,J2)=C(I,J1,J2)
	  CALL CMULT(X1,B0,X3,N)
	  CALL MULT(X3,X2,X1,N)
	  DO 10 J1=1,N
	  DO 10 J2=1,N
10      P(I,J1,J2)=X1(J1,J2)
5       CONTINUE        
C
C----- For up-going wave : Q+(m)+_,S(m)+        
	DO 15 I=NBOU,1,-1
	  IFLAG=1
	  IF(I.EQ.NBOU) IFLAG=0
	  DO 17 J1=1,N
	  IF(I.NE.1) E1(J1)=E(I-1,J1)
	  V1(J1)=V(I-1,J1)
	  V2(J1)=V(I,J1)
	  DO 17 J2=1,N
17      P0(J1,J2)=P(I,J1,J2)
	  CALL TRAN(P0,V1,V2,X1,X2,X3,IFLAG,N)
	  DO 20 J1=1,N
	  DO 20 J2=1,N
	  QU(I,J1,J2)=X2(J1,J2)
20      SU(I,J1,J2)=X3(J1,J2)
	  IF(I.EQ.1) GOTO 15
	  CALL AMULT(E1,X2,X3,0,N)
	  CALL AMULT(E1,X3,X1,1,N)
15      CONTINUE
C
C----- For down-going wave : Q-(m)+_,S(m)-        
	DO 27 I=1,NBOU
	  IFLAG=1
	  IF(I.EQ.1) IFLAG=0
	  DO 30 J1=1,N
	  IF(I.NE.NBOU) E1(J1)=E(I,J1)
	  V1(J1)=V(I,J1)
	  V2(J1)=V(I-1,J1)
	  DO 30 J2=1,N
30      P0(J1,J2)=P(I,J2,J1)
	  CALL TRAN(P0,V1,V2,X1,X2,X3,IFLAG,N)
	  DO 32 J1=1,N
	  DO 32 J2=1,N
	  QD(I,J1,J2)=X2(J1,J2)
32      SD(I,J1,J2)=X3(J1,J2)
	  IF(I.EQ.NBOU) GOTO 27
	  CALL AMULT(E1,X2,X3,0,N)
	  CALL AMULT(E1,X3,X1,1,N)
27      CONTINUE
C	
	END