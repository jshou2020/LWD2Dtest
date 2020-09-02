	SUBROUTINE XCAL(R1,RX,R2,X,IFLAG)
	USE GlobalVarbls
	IMPLICIT NONE
C
	INTEGER IFLAG,I
	REAL*8 X(0:MAXORDER),R1,R2,RX,DR
C
	IF(IFLAG.EQ.0) THEN
	  IF(R1.EQ.0.0D0) THEN
	    X(0)=0.0D0
	    DO 5 I=1,MAXORDER
5           X(I)=1.0D0/DBLE(I)
	  ELSE
	    DR=R2-R1
	    X(0)=DLOG(R2/R1)
	    DO 10 I=1,MAXORDER
10          X(I)=1.0D0/DBLE(I)-X(I-1)*R1/DR
	  END IF
	ELSE
	  IF(R1.EQ.0.0D0) THEN
	    DR=R2-R1
	    X(0)=0.0D0
	    DO 15 I=1,MAXORDER
15          X(I)=(RX/DR)**I/DBLE(I)
	  ELSE
	    DR=R2-R1
	    X(0)=DLOG(RX/R1)
	    DO 20 I=1,MAXORDER
20        X(I)=((RX-R1)/DR)**I/DBLE(I)-X(I-1)*R1/DR
	  END IF
	END IF
C
	END