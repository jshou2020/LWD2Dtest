	SUBROUTINE ASEMB1
C-----------------------------------------------------------------------
C Subroutine to assemble matrix B and first part of matrix A. 
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	REAL*8 XL(0:MAXORDER),XR(0:MAXORDER),R1,R2,RX,DRL,DR2L,DRR,DR2R
	INTEGER I,J
C
	DO 5 I=1,N
	DO 5 J=1,N
	A(I,J)=0.0D0
5     B(I,J)=0.0D0
C
	DO 10 I=1,N
	  R1=R(I-1)
	  R2=R(I)
	  DRL=R2-R1
	  DR2L=DRL*DRL
	  CALL XCAL(R1,RX,R2,XL,0)
	  R1=R(I)
	  R2=R(I+1)
	  DRR=R2-R1
	  DR2R=DRR*DRR
	  CALL XCAL(R1,RX,R2,XR,0)
	  IF(I.EQ.1) THEN
	  A(I,I)=4.0D0/DR2L*XL(2)+XR(0)/DR2R
	  B(I,I)=XL(4)+XR(2)-2.0D0*XR(1)+XR(0)
	  ELSE
	  A(I,I)=XL(0)/DR2L+XR(0)/DR2R
	  B(I,I)=XL(2)+XR(2)-2.0D0*XR(1)+XR(0)
	  END IF
	  IF(I.EQ.N) GOTO 15
	  A(I,I+1)=-XR(0)/DR2R
	  B(I,I+1)=-XR(2)+XR(1)
	  B(I+1,I)=B(I,I+1)
	  A(I+1,I)=A(I,I+1)
15      CONTINUE
10    CONTINUE        
	
	END 
