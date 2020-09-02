	SUBROUTINE GRID
C-----------------------------------------------------------------------
C Subroutine to generate the grid for the one-dimensional eigenvalue 
C problem. Here, N is the number of inner grid points.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	INTEGER I,K
	!REAL TEMP
	REAL(8):: TEMP
!===========================================================
C	
	R(0)=0.0D0; TEMP=RMAX*(1-Q)/(1-Q**N)
	DO 2 I=0,N-1
      R(I+1)=R(I)+TEMP*Q**I
      !WRITE(*,*)I,R(I)
2     CONTINUE;  !PAUSE
C
      !WRITE(*,*)RADIUST,RMAX;PAUSE
	I=1
	DO 5 WHILE(R(I).LT.RADIUST) 
5     I=I+1  
      !WRITE(*,*)I;PAUSE
C	
	DO 7 K=N,I,-1
7       R(K+1)=R(K)
C	
	R(I)=RADIUST
      RETURN
      OPEN(300, FILE = 'GRID.OUT')
      WRITE(300,*)'N = ', N;  WRITE(300,*)'Q = ', Q
      WRITE(300,*)'RMAX = ', RMAX;  WRITE(300,*)'TEMP = ', TEMP
	do i=1, N+1; write(300,*)i,R(i); ENDDO; CLOSE(100)
    	END