!
!**********************************************************************************
!     File name    :  SOLVE.for   
!     DATE         :  Wed, 09/02/2020
!     LAST UPDATE  :  Wed, 09/02/2020
!**********************************************************************************
	SUBROUTINE SOLVE_EWR(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep)
C-----------------------------------------------------------------------
	USE GlobalVarbls; use Tool5_par,only: NARRAYzz5;
      use Tool52_par,only: NARRAYzz52;  IMPLICIT NONE
C
	COMPLEX*16 :: RE, VRTEMP, Rearray(12,12), Vtemp(2), VR12
	INTEGER :: I, JT, JR, NT, NR, IFLAG, NSpacing, NARRAYzz, J
	REAL*8 :: ZS, ZT, ZR
	REAL(8):: Phase(6),Att(6); !REAL(8):: AllPhase(300,6),AllAtt(300,6);
      REAL(8):: Pha2(2,6),Att2(2,6);
      REAL(8):: AllPhase(IMAX,6),AllAtt(IMAX,6);
	REAL(8):: RaPhase(6),RaAtt(6);	!REAL(8):: RaAllPhase(300,6),RaAllAtt(300,6);
	REAL(8):: RaAllPhase(IMAX,6),RaAllAtt(IMAX,6), logDep(IMAX);
	INTEGER:: IFREQ, JdownUp
	REAL(8):: Freq
C===========================================================================
!     
      Freq=OMEGA/(2*PI)/1000.0; 
      IF(TRIM(toolname)=='EWR-M5')THEN
        NARRAYzz=NARRAYzz5; JdownUp=2;
      ELSE IF(TRIM(toolname)=='EWR-P4')THEN
        NARRAYzz=NARRAYzz52; JdownUp=1;   
      ENDIF
      !
	IFLAG = 0; NSpacing = 20
	IF(IMAX > 100) THEN
         NSpacing = 50
	ELSEIF(IMAX > 1000) THEN
         NSpacing = 100
	END IF
C
	DO 2 I = 1, IMAX
	   ZS = SAPOINT + (I - 1) * DELTA; VRTEMP = CZERO
         logDep(I)= ZS / FFOOT
         !--------------------------------------------
         DO J=1,JdownUp !2
	   DO  JT = 1, NTT
	      IF(J==1)ZT = ZS + PTT(JT)  !T-DOWN
	      IF(J==2)ZT = ZS - PTT(JT)  !T-UP
	      CALL POSITION(ZT, NT);  CALL UTRAN(ZT, NT)
	      DO  JR = 1, NRR
	         ZR = ZS + PRR(JR)
	         CALL POSITION(ZR, NR); CALL RECAL(ZT, ZR, NT, NR, RE)
	         VRTEMP = VRTEMP + MMT(JT) * MMR(JR) * RE; !WRITE(*,*)RE
	         Rearray(JT,JR)=RE
!5              CONTINUE
               Vtemp(JR)=RE
            ENDDO !JR LOOP
             IF(J==1)VR12=Vtemp(2)/Vtemp(1);  !T-DOWN - 2/1
             IF(J==2)VR12=Vtemp(1)/Vtemp(2);  !T-UP   - 1/2
             VR12=CONJG(VR12);
             CALL ComplexV_2PhaseAtt (VR12, pha2(J,JT),Att2(J,JT))
         ENDDO !JT LOOP
         ENDDO !JLOOP
         IF(JdownUp==2)THEN
          Phase=0.5*(pha2(1,:)+pha2(2,:)); Att=0.5*(Att2(1,:)+Att2(2,:))
         ELSE IF(JdownUp==1)THEN
          Phase=pha2(1,:); Att=Att2(1,:);
         ENDIF
         !-----------------------------------------------------------------------
         !write(*,*)NTT,NRR
C-----------------------------------------------------------------------
C Unit for apparent conductivity is mS/m, and unit for depth is foot.
C-----------------------------------------------------------------------
	   !
	   VRTEMP = DTWO * PI * VRTEMP / TTK * 1.0D3;  VR(I) = CONJG(VRTEMP); !cycle
	   !
	   !Rearray=CONJG(Rearray);
	   !V -- > Phase, att
	   !CALL Compute_3coilPhaseAtt(Rearray(1,:),Phase,Att);!WRITE(*, 10)Phase(1:3),Att(1:3)
	   AllPhase(I,1:NARRAYzz)=Phase(1:NARRAYzz); 
         AllAtt(I,1:NARRAYzz)=Att(1:NARRAYzz)
	   !Phase, att  --> Ra
         IF(TRIM(toolname)=='EWR-M5')THEN
	    IF(ABS(Freq-250.0)<1.0D-8)IFREQ=1;
          IF(ABS(Freq-500.0)<1.0D-8)IFREQ=2;
          IF(ABS(Freq-2000.0)<1.0D-8)IFREQ=3;
         ELSE IF(TRIM(toolname)=='EWR-P4')THEN
	    IF(ABS(Freq-1000.0)<1.0D-8)IFREQ=1; 
	    IF(ABS(Freq-2000.0)<1.0D-8)IFREQ=2;
         ENDIF
	   CALL PhaseAtt_2Ra(IFREQ, Phase,Att, RaPhase,RaAtt); !WRITE(*, 10)RaPhase(1:3),RaAtt(1:3)
	   RaAllPhase(I,1:NARRAYzz)=RaPhase(1:NARRAYzz); 
         RaAllAtt(I,1:NARRAYzz)=RaAtt(1:NARRAYzz)
         !
	   IF(MOD(I, NSpacing) == 0) THEN
	      IF(IFLAG == 0) THEN; !WRITE(*,'(1X, 20HLOGGING POINTS LEFT:)')
	         PRINT*;           !WRITE(*,'(12X, I8)') IMAX - I
		     IFLAG = 1
            ELSE
              ! WRITE(*,'(12X, I8)') IMAX - I
            END IF 
	   END IF
2     CONTINUE 
      RETURN
      !==============================================================================
      !
	OPEN(UNIT,FILE = 'IND.OUT')
	   DO 7 I = 1, IMAX
	      ZS = SAPOINT + (I - 1) * DELTA;   ZS = ZS / FFOOT
            WRITE(*, 10) ZS, VR(I);   WRITE(UNIT, 10) ZS, VR(I)
10          FORMAT(1X, F16.7, 20F16.7)
7	  CONTINUE
	  CLOSE(UNIT)
C
      Freq=OMEGA/(2*PI)/1000
      
	IF(ABS(Freq-2000.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF2000.OUT')
	IF(ABS(Freq-500.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF500.OUT')
	IF(ABS(Freq-250.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF250.OUT')
	   DO 77 I = 1, IMAX
	      ZS = SAPOINT + (I - 1) * DELTA
	      ZS = ZS / FFOOT
            !WRITE(*, 10) ZS, VR(I); 
            !WRITE(UNIT, 10) ZS, VR(I)
            WRITE(*, 10) ZS, AllPhase(I,1:NARRAYzz),
     1                       AllAtt(I,1:NARRAYzz) ; 
            WRITE(UNIT, 110)ZS, AllPhase(I,1:NARRAYzz),
     1                          AllAtt(I,1:NARRAYzz) ;
110          FORMAT(1X, F16.7, 20F16.7)
77	  CONTINUE
	  CLOSE(UNIT)
!
	IF(ABS(Freq-2000.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF2000.OUT')
	IF(ABS(Freq-500.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF500.OUT')
	IF(ABS(Freq-250.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF250.OUT')
	   DO I = 1, IMAX
	      ZS = SAPOINT + (I - 1) * DELTA
	      ZS = ZS / FFOOT
            WRITE(*, 10) ZS, RaAllPhase(I,1:NARRAYzz),
     1                       RaAllAtt(I,1:NARRAYzz) ; 
            WRITE(UNIT, 110)ZS, RaAllPhase(I,1:NARRAYzz),
     1                          RaAllAtt(I,1:NARRAYzz) ; 
	  ENDDO
	  CLOSE(UNIT)
	END SUBROUTINE SOLVE_EWR
!	
!**********************************************************************************
	SUBROUTINE SOLVE_ARC(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep)
	USE GlobalVarbls; !use Tool5_par,only: NARRAYzz5
      use Tool6_par,only: NARRAYzz6
	IMPLICIT NONE
C
	COMPLEX*16 :: RE, VRTEMP, Rearray(12,12), Vtemp(2), VR12
	INTEGER :: I, JT, JR, NT, NR, IFLAG, NSpacing, NARRAYzz, J
	REAL*8 :: ZS, ZT, ZR
	REAL(8):: Phase(6),Att(6); REAL(8):: Pha2(2,6),Att2(2,6);
      REAL(8):: AllPhase(IMAX,6),AllAtt(IMAX,6); 
	REAL(8):: RaPhase(6),RaAtt(6);	;
	REAL(8):: RaAllPhase(IMAX,6),RaAllAtt(IMAX,6), logDep(IMAX);
	INTEGER:: IFREQ; REAL(8):: Freq
C===========================================================================
!     
      Freq=OMEGA/(2*PI)/1000.0; NARRAYzz=NARRAYzz6;  !NARRAYzz5
C
	DO 2 I = 1, IMAX
	   ZS = SAPOINT + (I - 1) * DELTA; VRTEMP = CZERO
         logDep(I)= ZS / FFOOT
         !--------------------------------------------
         DO J=1,1 !2
	   DO  JT = 1, NTT
	      IF(J==1)ZT = ZS + PTT(JT)  !T-DOWN
	      IF(J==2)ZT = ZS - PTT(JT)  !T-UP
	      CALL POSITION(ZT, NT);  CALL UTRAN(ZT, NT)
	      DO  JR = 1, NRR
	         ZR = ZS + PRR(JR)
	         CALL POSITION(ZR, NR); CALL RECAL(ZT, ZR, NT, NR, RE)
	         VRTEMP = VRTEMP + MMT(JT) * MMR(JR) * RE; !WRITE(*,*)RE
	         Rearray(JT,JR)=RE
!5              CONTINUE
               Vtemp(JR)=RE
            ENDDO !JR LOOP
             !IF(J==1)VR12=Vtemp(2)/Vtemp(1);  !T-DOWN - 2/1
             !IF(J==2)VR12=Vtemp(1)/Vtemp(2);  !T-UP   - 1/2
             IF(PTT(JT)>0.D0)VR12=Vtemp(2)/Vtemp(1);
             IF(PTT(JT)<0.D0)VR12=Vtemp(1)/Vtemp(2);
             VR12=CONJG(VR12);
             CALL ComplexV_2PhaseAtt (VR12, pha2(J,JT),Att2(J,JT))
         ENDDO !JT LOOP
         ENDDO !JLOOP
         !Phase=0.5*(pha2(1,:)+pha2(2,:)); Att=0.5*(Att2(1,:)+Att2(2,:))
         Phase=pha2(1,:); Att=Att2(1,:)
         !-----------------------------------------------------------------------
         !write(*,*)NTT,NRR
	   !V -- > Phase, att
	   !CALL Compute_3coilPhaseAtt(Rearray(1,:),Phase,Att);!WRITE(*, 10)Phase(1:3),Att(1:3)
	   AllPhase(I,1:NARRAYzz)=Phase(1:NARRAYzz); 
         AllAtt(I,1:NARRAYzz)=Att(1:NARRAYzz)
	   !Phase, att  --> Ra
	   IF(ABS(Freq-400.)<1.0D-8)IFREQ=1; IF(ABS(Freq-2000.)<1.0D-8)IFREQ=2;
	   CALL PhaseAtt_2Ra(IFREQ, Phase,Att, RaPhase,RaAtt); !WRITE(*, 10)RaPhase(1:3),RaAtt(1:3)
	   RaAllPhase(I,1:NARRAYzz)=RaPhase(1:NARRAYzz); 
         RaAllAtt(I,1:NARRAYzz)=RaAtt(1:NARRAYzz)
         !
2     CONTINUE 
      !RETURN
	END SUBROUTINE SOLVE_ARC
!	
!**********************************************************************************
	SUBROUTINE SOLVE_MPR(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep)
	USE GlobalVarbls; !use Tool5_par,only: NARRAYzz5
      use Tool7_par,only: NARRAYzz7; !use Tool6_par,only: NARRAYzz6
	IMPLICIT NONE
C
	COMPLEX*16 :: RE, VRTEMP, Rearray(12,12), Vtemp(2), VR12
	INTEGER :: I, JT, JR, NT, NR, IFLAG, NSpacing, NARRAYzz, J
	REAL*8 :: ZS, ZT, ZR
	REAL(8):: Phase(6),Att(6); REAL(8):: Pha2(2,6),Att2(2,6);
      REAL(8):: AllPhase(IMAX,6),AllAtt(IMAX,6); 
	REAL(8):: RaPhase(6),RaAtt(6);	;
	REAL(8):: RaAllPhase(IMAX,6),RaAllAtt(IMAX,6), logDep(IMAX);
	INTEGER:: IFREQ; REAL(8):: Freq
C===========================================================================
!     
      Freq=OMEGA/(2*PI)/1000.0; NARRAYzz=NARRAYzz7;  !NARRAYzz5, NARRAYzz6
C
	DO 2 I = 1, IMAX
	   ZS = SAPOINT + (I - 1) * DELTA; VRTEMP = CZERO
         logDep(I)= ZS / FFOOT
         !--------------------------------------------
         DO J=1,2
	   DO  JT = 1, NTT
	      IF(J==1)ZT = ZS + PTT(JT)  !T-DOWN
	      IF(J==2)ZT = ZS - PTT(JT)  !T-UP
	      CALL POSITION(ZT, NT);  CALL UTRAN(ZT, NT)
	      DO  JR = 1, NRR
	         ZR = ZS + PRR(JR)
	         CALL POSITION(ZR, NR); CALL RECAL(ZT, ZR, NT, NR, RE)
	         VRTEMP = VRTEMP + MMT(JT) * MMR(JR) * RE; !WRITE(*,*)RE
	         Rearray(JT,JR)=RE
!5              CONTINUE
               Vtemp(JR)=RE
            ENDDO !JR LOOP
             IF(J==1)VR12=Vtemp(2)/Vtemp(1);  !T-DOWN - 2/1
             IF(J==2)VR12=Vtemp(1)/Vtemp(2);  !T-UP   - 1/2
             VR12=CONJG(VR12);
             CALL ComplexV_2PhaseAtt (VR12, pha2(J,JT),Att2(J,JT))
         ENDDO !JT LOOP
         ENDDO !JLOOP
         Phase=0.5*(pha2(1,:)+pha2(2,:)); Att=0.5*(Att2(1,:)+Att2(2,:))
         !-----------------------------------------------------------------------
         !write(*,*)NTT,NRR
	   !V -- > Phase, att
	   !CALL Compute_3coilPhaseAtt(Rearray(1,:),Phase,Att);!WRITE(*, 10)Phase(1:3),Att(1:3)
	   AllPhase(I,1:NARRAYzz)=Phase(1:NARRAYzz); 
         AllAtt(I,1:NARRAYzz)=Att(1:NARRAYzz)
	   !Phase, att  --> Ra
	   IF(ABS(Freq-400.)<1.0D-8)IFREQ=1; IF(ABS(Freq-2000.)<1.0D-8)IFREQ=2;
	   CALL PhaseAtt_2Ra(IFREQ, Phase,Att, RaPhase,RaAtt); !WRITE(*, 10)RaPhase(1:3),RaAtt(1:3)
	   RaAllPhase(I,1:NARRAYzz)=RaPhase(1:NARRAYzz); 
         RaAllAtt(I,1:NARRAYzz)=RaAtt(1:NARRAYzz)
         !
2     CONTINUE 
      !RETURN
	END SUBROUTINE SOLVE_MPR
!	
!**********************************************************************************
!
	SUBROUTINE SOLVE_glw
C-----------------------------------------------------------------------
C Subroutine to compute apparent conductivity.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	COMPLEX*16 :: RE, VRTEMP
	INTEGER :: I, JT, JR, NT, NR, IFLAG, NSpacing
	REAL*8 :: ZS, ZT, ZR
C
	IFLAG = 0
	NSpacing = 20
	IF(IMAX > 100) THEN
         NSpacing = 50
	ELSEIF(IMAX > 1000) THEN
         NSpacing = 100
	END IF
C
	DO 2 I = 1, IMAX
	   ZS = SAPOINT + (I - 1) * DELTA
	   VRTEMP = CZERO
	   DO 5 JT = 1, NTT
	      ZT = ZS + PTT(JT)
	      CALL POSITION(ZT, NT) 
	      CALL UTRAN(ZT, NT)
	      DO 5 JR = 1, NRR
	         ZR = ZS + PRR(JR)
	         CALL POSITION(ZR, NR)
	         CALL RECAL(ZT, ZR, NT, NR, RE)
	         VRTEMP = VRTEMP + MMT(JT) * MMR(JR) * RE
5        CONTINUE
C-----------------------------------------------------------------------
C Unit for apparent conductivity is mS/m, and unit for depth is foot.
C-----------------------------------------------------------------------
	   VRTEMP = DTWO * PI * VRTEMP / TTK * 1.0D3
	   VR(I) = CONJG(VRTEMP)
C
	   IF(MOD(I, NSpacing) == 0) THEN
	      IF(IFLAG == 0) THEN
               WRITE(*,'(1X, 20HLOGGING POINTS LEFT:)')
	         PRINT*
            WRITE(*,'(12X, I8)') IMAX - I
		     IFLAG = 1
            ELSE
               WRITE(*,'(12X, I8)') IMAX - I
            END IF 
	   END IF
2     CONTINUE        
C
	OPEN(UNIT,FILE = 'IND.OUT')
	   DO 7 I = 1, IMAX
	      ZS = SAPOINT + (I - 1) * DELTA
	      ZS = ZS / FFOOT
            WRITE(*, 10) ZS, VR(I)
            WRITE(UNIT, 10) ZS, VR(I)
10          FORMAT(1X, F16.7, 2F16.7)
7	  CONTINUE
	  CLOSE(UNIT)
C
      END	
!
!***********************************************************************************
c   The End 
!***********************************************************************************
!