	SUBROUTINE TRAN(P,V1,V2,Q1,Q2,S,IFLAG,N)   
C-----------------------------------------------------------------------
C Subroutine to calculate the generalized reflection matrix and 
C generalized transmission matrix in the inner layer using the recursive 
C formula when the generalized reflection matrix in the outer layer is
C known.
C   P	   :	 Pm,m+_1
C   V1   :	 Vm
C   V2   :	 Vm+_1
C   Q1   : Q+_(m+_1)_+
C   Q2   : Q+_(m)_+
C   S    : S(m)+_
C   IFLAG: indicator, =0, means m=0 or m=NBOU
C-----------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER N,IFLAG,I,J
	COMPLEX*16 P(N,N),Q1(N,N),V1(N),V2(N)
	COMPLEX*16 Q2(N,N),S(N,N),X(N,N)
C
	DO 2 I=1,N
	V1(I)=(1.0D0,0.0D0)/V1(I)
	IF(IFLAG.NE.0) THEN
	  DO 3 J=1,N
	  IF(I.EQ.J) THEN
	  S(I,J)=(1.0D0,0.0D0)-Q1(I,J)
	  Q1(I,J)=(1.0D0,0.0D0)+Q1(I,J)
	  ELSE
	  S(I,J)=-Q1(I,J)
	  END IF
3       CONTINUE
	END IF
2     CONTINUE
C	
	IF(IFLAG.NE.0) THEN
	  CALL MULT(P,Q1,X,N)
	  CALL AMULT(V2,S,Q2,0,N)
	  CALL MULT(P,Q2,S,N)
	  CALL AMULT(V1,S,Q2,0,N)
	ELSE
	  CALL AMULT(V1,P,S,0,N)
	  CALL AMULT(V2,S,Q2,1,N)
	  DO 4 I=1,N
	  DO 4 J=1,N
4       X(I,J)=P(I,J)
	END IF
C
	DO 5 I=1,N
	DO 5 J=1,N
	S(I,J)=X(I,J)+Q2(I,J)
5     X(I,J)=X(I,J)-Q2(I,J)
C
	CALL INV(S,N)
	CALL MULT(X,S,Q2,N)
C	
	DO 7 I=1,N
	DO 7 J=1,N
7     S(I,J)=(2.0D0,0.0D0)*S(I,J)
C
	END