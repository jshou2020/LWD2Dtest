	SUBROUTINE COPYRIGHT
C-----------------------------------------------------------------------
C Subroutine to generate the copyright information.
C-----------------------------------------------------------------------
	CHARACTER*1, X1(100),X2(100)
	CHARACTER*18, TITLE
	CHARACTER*22, CP
	CHARACTER*14, ONWER
	CHARACTER*14, TIME
	DATA TITLE/'  INDMod2DMed 1.0 '/
	DATA CP/'Copyright @ 1997-2008 '/
	DATA ONWER/'Gong Li Wang. '/
	DATA TIME/'Dec 02, 2008  '/
C
	DO 2 I=1,100
	X1(I)=CHAR(205)
	X2(I)=CHAR(32)
2	CONTINUE
C
	PRINT*
	PRINT*
	WRITE(*,*) (X1(I),I=1,76)
	WRITE(*,*) (X1(I),I=1,4), TITLE,CP,ONWER,TIME,(X1(I),I=1,4)
	WRITE(*,*) (X1(I),I=1,76)
	PRINT*
	PRINT*
C
	END