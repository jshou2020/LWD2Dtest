	SUBROUTINE ASEMB2(XA,NBED)
C-----------------------------------------------------------------------
C Subroutine to assemble the second part of matrix A. 
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	COMPLEX*16 XA(N,N)
	INTEGER NBED
C
	REAL*8 XL(0:MAXI,0:MAXORDER),XR(0:MAXI,0:MAXORDER)
	REAL*8 X1(0:MAXORDER),X2(0:MAXORDER)
C
	INTEGER I,J,K
	REAL*8 R1,R2,RX,DRL,DR2L,DRR,DR2R
C
	DO 5 I=1,N
	DO 5 J=1,N
5     XA(I,J)=(0.0D0,0.0D0)
C
	DO 10 I=1,N
	  R1=R(I-1)
	  R2=R(I)
	  CALL XCAL(R1,R2,R2,X2,0)
	  IF(IN(NBED,I-1).EQ.0) THEN
	    DO 15 K=0,MAXORDER
15        XL(0,K)=X2(K)
	  ELSE
	    DO 20 J=IN(NBED,I-1),0,-1
	      IF(J.EQ.0) THEN
	        DO 25 K=0,MAXORDER
25            XL(J,K)=X2(K)
	      ELSE
	        RX=RM(NBED,I-1,J)
	        CALL XCAL(R1,RX,R2,X1,1)
	        DO 30 K=0,MAXORDER
	        XL(J,K)=X2(K)-X1(K)
30            X2(K)=X1(K)      
	      END IF
20        CONTINUE
	  END IF
	  DRL=R(I)-R(I-1)
	  DR2L=DRL*DRL
C	  
	  R1=R(I)
	  R2=R(I+1)
	  CALL XCAL(R1,R2,R2,X2,0)
	  IF(IN(NBED,I).EQ.0) THEN
	    DO 35 K=0,MAXORDER
35        XR(0,K)=X2(K)
	  ELSE
	    DO 40 J=IN(NBED,I),0,-1
	      IF(J.EQ.0) THEN
	        DO 45 K=0,MAXORDER
45            XR(J,K)=X2(K)
	      ELSE
	        RX=RM(NBED,I,J)
	        CALL XCAL(R1,RX,R2,X1,1)
	        DO 50 K=0,MAXORDER
	        XR(J,K)=X2(K)-X1(K)
50            X2(K)=X1(K)      
	      END IF
40        CONTINUE
	  END IF
	  DRR=R(I+1)-R(I)
	  DR2R=DRR*DRR
C	  
	  DO 55 J=0,IN(NBED,I-1)
	  IF(I.EQ.1) THEN
	    XA(I,I)=XA(I,I)+K2(NBED,I-1,J)*XL(J,4)
	  ELSE
	    XA(I,I)=XA(I,I)+K2(NBED,I-1,J)*XL(J,2)
	  END IF
55        CONTINUE
	  DO 60 J=0,IN(NBED,I)
	    XA(I,I)=XA(I,I)+K2(NBED,I,J)*(XR(J,2)-2.0D0*XR(J,1)+XR(J,0))
60        IF(I.NE.N) XA(I,I+1)=XA(I,I+1)+K2(NBED,I,J)*(-XR(J,2)+XR(J,1))
65        IF(I.NE.N) XA(I+1,I)=XA(I,I+1)

10      CONTINUE        
C
	END 