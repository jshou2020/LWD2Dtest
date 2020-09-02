	SUBROUTINE GCAL(X,GF)
C-----------------------------------------------------------------------
C Subroutine to evaluate the vector of basis functions at a given point.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	REAL*8 GF(N),L1,L2
	INTEGER I
	REAL*8 R1,R2,X,X1,X2
C	 
	L1(X,X1,X2)=(X2-X)/(X2-X1)
	L2(X,X1,X2)=(X-X1)/(X2-X1)
C	
	DO 2 I=1,N
2     GF(I)=0.0D0
C	
	I=1
	DO 5 WHILE(R(I).LT.X) 
5     I=I+1  
C	
	R1=R(I-1)
	R2=R(I)
	IF(I.EQ.1) THEN
	  GF(I)=L2(X,R1,R2)*L2(X,R1,R2)
	ELSE
	  GF(I-1)=L1(X,R1,R2)
	  IF(I.LE.N) GF(I)=L2(X,R1,R2)
	END IF
C	
	END