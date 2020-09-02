	SUBROUTINE PREDT
C-----------------------------------------------------------------------
C Subroutine to prepare electrical properties for the evaluation of A 
C and B matrices.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	COMPLEX*16, DIMENSION(:, :), ALLOCATABLE :: K2Temp
	INTEGER, DIMENSION(:), ALLOCATABLE :: NSECTR
C
	INTEGER :: I, J, K, IX, NSEC
	REAL*8 :: K02
C
	ALLOCATE(K2Temp(0 : NBOU,0 : MAXS), NSECTR(0 : NBOU))
C	
	K02=OMEGA*OMEGA*U0
	DO 2 I=0,NBOU
	DO 2 J=0,MAXS
	  K2Temp(I,J)=K02*DCMPLX(EPS(I,J)*E0,-SIG(I,J)/OMEGA)
2	CONTINUE
C
	DO 5 I=0,NBOU 
	  J=1
        IF(J.GT.MAXS) GOTO 10
	  DO 7 WHILE(SIG(I,J).GT.0.0D0)
	    J=J+1
          IF(J.GT.MAXS) GOTO 10
7	  CONTINUE
10      NSECTR(I)=J-1 ! Number of electrical interfaces in a layer.
5	CONTINUE
C
	DO 12 I=0,NBOU
	  DO 12 J=0,N
	    IX=0
	    DO 15 K=1,NSECTR(I)
	      IF(VECOD(I,K).GE.R(J).AND.VECOD(I,K).LE.R(J+1)) THEN
	        IX=IX+1
	        IF(IX.GT.MAXI) 
     #		STOP 'THE NUMBER OF DEFAULT SECTIONS IN AN ELEMENT IS SMALL!'
	        RM(I,J,IX)=VECOD(I,K)
	        K2(I,J,IX-1)=K2Temp(I,K-1)
	        K2(I,J,IX)=K2Temp(I,K)
	      END IF
15        CONTINUE  
C
	    IN(I,J)=IX
	    IF(IX.EQ.0) THEN
	      NSEC=NSECTR(I)
	      IF(R(J+1).GE.VECOD(I,NSEC)) THEN
		    K2(I,J,0)=K2Temp(I,NSEC)
	      ELSE
		    K=0
		    DO 17 WHILE(R(J+1).GT.VECOD(I,K))
17              K=K+1
		    K2(I,J,0)=K2Temp(I,K-1)
	      END IF
	    END IF  
12      CONTINUE           
C
	DEALLOCATE(K2Temp, NSECTR)
C
	END