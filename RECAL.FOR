	SUBROUTINE RECAL(ZT,ZR,NT,NR,RE)        
C-----------------------------------------------------------------------
C Subroutine to evaluate r*E at the observation.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	REAL*8 ZT,ZR
	INTEGER NT,NR
	COMPLEX*16 RE
C
	COMPLEX*16 V1(N),V2(N),V3(N),X1(N,N),X2(N,N)
C
	INTEGER I,J
	REAL*8 H
C
	CALL GCAL(RADIUSR,G)
	DO 2 I=1,N
	DO 2 J=1,N
2       X1(I,J)=C(NR,I,J)
	CALL DMULT(G,X1,V2,N)
C
	IF(NT.EQ.NR) THEN
	  H=ZT-ZR
	  IF(H.LE.0) THEN
	    IF(NR.NE.NBOU) THEN
	      DO 5 I=1,N
	      V1(I)=CDEXP(-V(NR,I)*(HOR(NR+1)-ZR))
	      V3(I)=CDEXP(-V(NR,I)*(ZR-ZT))*UP(I)
	      DO 5 J=1,N
5             X1(I,J)=QU(NR+1,I,J)
	      CALL AMULT(V1,X1,X2,0,N)
	      CALL AMULT(V1,X2,X1,1,N)
	      DO 7 I=1,N
7             X1(I,I)=(1.0D0,0.0D0)+X1(I,I)   
	      CALL BMULT(X1,V3,V1,N)
	      RE=(0.0D0,0.0D0)
	      DO 10 I=1,N
10            RE=RE+V1(I)*V2(I)
	    ELSE
	      RE=(0.0D0,0.0D0)
	      DO 12 I=1,N
12            RE=RE+CDEXP(-V(NR,I)*(ZR-ZT))*UP(I)*V2(I)
	    END IF
	  ELSE
	    IF(NR.NE.0) THEN
	      DO 15 I=1,N
	      V1(I)=CDEXP(-V(NR,I)*(ZR-HOR(NR)))
	      V3(I)=CDEXP(-V(NR,I)*(ZT-ZR))*UD(I)
	      DO 15 J=1,N
15            X1(I,J)=QD(NR,I,J)
	      CALL AMULT(V1,X1,X2,0,N)
	      CALL AMULT(V1,X2,X1,1,N)
	      DO 17 I=1,N
17            X1(I,I)=(1.0D0,0.0D0)+X1(I,I)   
	      CALL BMULT(X1,V3,V1,N)
	      RE=(0.0D0,0.0D0)
	      DO 20 I=1,N
20            RE=RE+V1(I)*V2(I)
	    ELSE
	      RE=(0.0D0,0.0D0)
	      DO 22 I=1,N
22            RE=RE+CDEXP(-V(NR,I)*(ZT-ZR))*UD(I)*V2(I)
	    END IF
	  END IF
	ELSEIF(NT.LT.NR) THEN
	  IF(NR.NE.NBOU) THEN
	    DO 25 I=1,N
	    V1(I)=CDEXP(-V(NR,I)*(HOR(NR+1)-ZR))
	    V3(I)=CDEXP(-V(NR,I)*(ZR-HOR(NR)))*U(NR,I)
	    DO 25 J=1,N
25          X1(I,J)=QU(NR+1,I,J)
	    CALL AMULT(V1,X1,X2,0,N)
	    CALL AMULT(V1,X2,X1,1,N)
	    DO 27 I=1,N
27          X1(I,I)=(1.0D0,0.0D0)+X1(I,I)   
	    CALL BMULT(X1,V3,V1,N)
	    RE=(0.0D0,0.0D0)
	    DO 30 I=1,N
30          RE=RE+V1(I)*V2(I)
	  ELSE
	    RE=(0.0D0,0.0D0)
	    DO 32 I=1,N
32          RE=RE+CDEXP(-V(NR,I)*(ZR-HOR(NR)))*U(NR,I)*V2(I)
	  END IF
	ELSE
	  IF(NR.NE.0) THEN
	    DO 35 I=1,N
	    V1(I)=CDEXP(-V(NR,I)*(ZR-HOR(NR)))
	    V3(I)=CDEXP(-V(NR,I)*(HOR(NR+1)-ZR))*U(NR+1,I)
	    DO 35 J=1,N
35          X1(I,J)=QD(NR,I,J)
	    CALL AMULT(V1,X1,X2,0,N)
	    CALL AMULT(V1,X2,X1,1,N)
	    DO 37 I=1,N
37          X1(I,I)=(1.0D0,0.0D0)+X1(I,I)   
	    CALL BMULT(X1,V3,V1,N)
	    RE=(0.0D0,0.0D0)
	    DO 40 I=1,N
40          RE=RE+V1(I)*V2(I)
	  ELSE
	    RE=(0.0D0,0.0D0)
	    DO 42 I=1,N
42          RE=RE+CDEXP(-V(NR,I)*(HOR(NR+1)-ZR))*U(NR+1,I)*V2(I)
	  END IF
	END IF
C
	END