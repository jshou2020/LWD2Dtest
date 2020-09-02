	SUBROUTINE MULT(A,B,C,N)
	IMPLICIT NONE
	INTEGER N,I,J,K
	COMPLEX*16 A(N,N),B(N,N),C(N,N)
C	
	DO 5 I=1,N
	  DO 5 J=1,N
	  C(I,J)=(0.0D0,0.0D0)
	  DO 5 K=1,N
5       C(I,J)=C(I,J)+A(I,K)*B(K,J)
C	
	END
C	
C	
	SUBROUTINE CMULT(A,B,C,N)
	IMPLICIT NONE
	INTEGER N,I,J,K
	COMPLEX*16 A(N,N),C(N,N)
	REAL*8 B(N,N)
C
	DO 5 I=1,N
	  DO 5 J=1,N
	  C(I,J)=(0.0D0,0.0D0)
	  DO 5 K=1,N
5       C(I,J)=C(I,J)+A(I,K)*B(K,J)
C	
	END
C
C		
	SUBROUTINE DMULT(A,B,C,N)
	IMPLICIT NONE
	INTEGER N,I,K
	COMPLEX*16 B(N,N),C(N)
	REAL*8 A(N)
C
	DO 5 I=1,N
	  C(I)=(0.0D0,0.0D0)
	  DO 5 K=1,N
5       C(I)=C(I)+A(K)*B(K,I)
C	
	END
C
C		
	SUBROUTINE RMULT(A,B,C,N)
	IMPLICIT NONE
	INTEGER N,I,J,K
	REAL*8 A(N,N),B(N,N),C(N,N)
C	
	DO 5 I=1,N
	  DO 5 J=1,N
	  C(I,J)=0.0D0
	  DO 5 K=1,N
5       C(I,J)=C(I,J)+A(I,K)*B(K,J)
C	
	END
C
C		
	SUBROUTINE AMULT(A,B,C,K,N)
C-----------------------------------------------------------------------
C IF K=0,C=diag(a1,a2,...,aN)*B
C IF K=1,C=B*diag(a1,a2,...,aN)
C-----------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER N,I,J,K
	COMPLEX*16 A(N),B(N,N),C(N,N)
C	
	DO 5 I=1,N
	  DO 5 J=1,N
	  IF(K.EQ.0) THEN
	  C(I,J)=A(I)*B(I,J)
	  ELSE
	  C(I,J)=A(J)*B(I,J)
	  END IF
5       CONTINUE        

	END 
C
C
	SUBROUTINE BMULT(A,B,C,N)
C-----------------------------------------------------------------------
C c=A*b
C-----------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER N,I,J
	COMPLEX*16 A(N,N),B(N),C(N)
	
	DO 5 I=1,N
	  C(I)=(0.0D0,0.0D0)
	  DO 5 J=1,N
5       C(I)=C(I)+A(I,J)*B(J)
	
	END 