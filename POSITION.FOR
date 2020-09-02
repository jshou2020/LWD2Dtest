	SUBROUTINE POSITION(Z0,N0)
	USE GlobalVarbls
	IMPLICIT NONE
C
	REAL*8 Z0
	INTEGER N0,K
C
	IF(NBOU == 0) THEN
	   N0 = 0
      ELSE
	   IF(Z0.GE.HOR(NBOU)) THEN
	      N0=NBOU
	   ELSE
	      K=1
	      DO 20 WHILE(HOR(K).LE.Z0)
20          K=K+1
	      N0=K-1
         END IF
	END IF
C	    
	END