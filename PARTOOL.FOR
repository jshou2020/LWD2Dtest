!
!*******************************************************************************
!     Filename: PARTOOL.FOR
!     DATE         : Wed, 09/02/2020
!     LAST UPDATE  : Wed, 09/02/2020
!*******************************************************************************
!
      SUBROUTINE PARTOOL_02(FREQ)    !SUBROUTINE PARTOOL_hal(FREQ)
C-----------------------------------------------------------------------
C Subroutine to read in tool configration.
C-----------------------------------------------------------------------
	USE GlobalVarbls; USE Tool5_par; USE Tool52_par; 
      use Tool6_par; use Tool7_par;
	IMPLICIT NONE
	REAL*8 :: PP(MAXT + MAXR), ITR(MAXT + MAXR), NN(MAXT + MAXR)
C
      !CHARACTER*10:: ToolName
	REAL*8 FREQ, RADIUS , T1, T2; 	INTEGER I, IT, IR
C
	!OPEN(UNIT, FILE = 'TOOL_HAL.INP', STATUS = 'OLD')
	   !READ(UNIT, *) FREQ     !tool frequency input in kHz
	   FREQ=1000.0*FREQ      !kKz - > Hz
	   !IF(FREQ < 1.0D3) STOP 'Frequency too small!'!IF(FREQ > 1.0D5) STOP 'Frequency too large!'
	   !IF(FREQ < 1.0D3) write(*,*)'Frequency too small!!IF(FREQ > 1.0D5) write(*,*)'Frequency too large!'
         
	   !READ(UNIT, *) RADIUS 
         !RADIUS = 0.025 !INCH
         RADIUS = 0.0125
	   !I = 1
!2        !CONTINUE
	   !READ(UNIT, *, END = 5) PP(I), ITR(I), NN(I)
	   !I = I + 1
	   !IF(I > MAXT + MAXR) STOP 'Too many coils!'
	   !GOTO 2
!5        !CLOSE(UNIT); !WRITE(*,*)I;PAUSE
C
      I=7   
      PP(1:I)=(/0,12,20,28,36,44,52/); ITR(1:I)=(/1,0,0,0,0,0,0/); NN=5;
      !I=7; !NOCL = I - 1
      NOCL = I;	NTT = 0; NRR = 0
	DO 7 I = 1, NOCL
	   IF(ITR(I) == 1) THEN
	      NTT = NTT + 1;  IF(NTT > MAXT) STOP 'Too many transmitters!'
	      PTT(NTT) = PP(I) * FINCH; MMT(NTT) = NN(I);  GOTO 7
	   ELSEIF(ITR(I) == 0) THEN
	      NRR = NRR + 1;  IF(NRR > MAXR) STOP 'Too many receivers!'
	      PRR(NRR) = PP(I) * FINCH; MMR(NRR) = NN(I);   GOTO 7
	   ELSE
	      STOP 'Coil flag has to be 0 or 1!'		  
	   END IF
7     CONTINUE
!
      IF(Itool==-4)THEN  !EWR-M5
        IF(TRIM(toolname)=='EWR-M5')THEN
         NTT = 5; NRR = 2;
         PTT(1:NTT)=(/16,24,32,40,48/); PTT = PTT * FINCH; MMT = 5;
         PRR(1:NRR)=(/-4,4/); PRR=PRR*FINCH; MMR=5;
         Rspacing5(1:NReceiver5)=(/12,20,28,36,44,52/)  !NReceiver5=6; =NTT+1
        ELSE IF(TRIM(toolname)=='EWR-P4')THEN
         NTT = 4; NRR = 2;
         PTT(1:NTT)=(/9,15,27,39/); PTT = PTT * FINCH; MMT = 5;
         PRR(1:NRR)=(/-3,3/); PRR=PRR*FINCH; MMR=5;
         Rspacing5(1:NReceiver52)=(/6,12,12,18,24,30,36,42/)  !NReceiver52=8; =2*NTT
        ENDIF
        !
      ELSE IF(Itool==-5)THEN  !ARC675/825/900; ARC312/475; IMPulse
        NTT = 5; 	NRR = 2;
        IF(TRIM(toolname)=='ARC675'.OR.TRIM(toolname)=='ARC825'
     1                             .OR.TRIM(toolname)=='ARC900')THEN
          PTT(1:NTT)=(/16,-22,28,-34,40/); 
          Rspacing6(1:NReceiver6)=(/13,19,25,31,37,43/)  !NReceiver6=6; =NTT+1
        ELSE IF(TRIM(toolname)=='ARC312'.OR.TRIM(toolname)=='ARC475'
     1                             .OR.TRIM(toolname)=='IMPulse')THEN
          PTT(1:NTT)=(/10,-16,22,-28,34/); 
          Rspacing6(1:NReceiver6)=(/7,13,19,25,31,37/)
        ENDIF
        PTT = PTT * FINCH; MMT = 5; 
        PRR(1:NRR)=(/-3,3/); PRR=PRR*FINCH; MMR=5;
        !
      ELSE IF(Itool==-6)THEN  !MPR; MPR95
        NTT = 2; 	NRR = 2;
        IF(TRIM(toolname)=='MPR')THEN
         PTT(1:NTT)=(/22.375,35.625/); 
         Rspacing7(1:NReceiver7)=(/18.375,26.375,31.625,39.625/) !NReceiver7=4; =2*NTT
        ELSE IF(TRIM(toolname)=='MPR95')THEN
         PTT(1:NTT)=(/31.375,44.625/); 
         Rspacing7(1:NReceiver7)=(/27.375,35.375,40.625,48.625/)
        ENDIF
        PTT = PTT * FINCH; MMT = 5;
        PRR(1:NRR)=(/-4,4/); PRR=PRR*FINCH; MMR=5;
      ENDIF
C
	OMEGA = 2.0D0 * PI * FREQ
	RADIUST = RADIUS * FINCH
	RADIUSR = RADIUS * FINCH
	S0T = PI * RADIUST * RADIUST
	S0R = PI * RADIUSR * RADIUSR
C
C-----------------------------------------------------------------------
C Compute tool constants.
C-----------------------------------------------------------------------
	TTK = DZERO
	T1 = - OMEGA ** 2 * U0 ** 2 * S0T * S0R * CURENT /(4.0D0 * PI)
	DO 10 IT = 1, NTT
	   DO 10 IR = 1, NRR
	      T2 = DABS(PTT(IT) - PRR(IR))
	      TLK(IT, IR) = T1 * MMT(IT) * MMR(IR) / T2
            TTK = TTK + TLK(IT, IR)
10    CONTINUE
C	
      END SUBROUTINE PARTOOL_02  !END SUBROUTINE PARTOOL_hal
!
!*******************************************************************************
!
      SUBROUTINE PARTOOL(NumFrequency,FreqArray); !SUBROUTINE PARTOOL(FREQ)
C-----------------------------------------------------------------------
C Subroutine to read in tool configration.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
	REAL*8 :: PP(MAXT + MAXR), ITR(MAXT + MAXR), NN(MAXT + MAXR)
C
	REAL*8 FREQ, RADIUS , T1, T2
	INTEGER I, IT, IR
      !CHARACTER*10:: ToolName
      INTEGER :: NumFrequency; REAL(8):: FreqArray(10)
C====================================================================
!
	OPEN(UNIT, FILE = 'TOOL.INP', STATUS = 'OLD')
         READ(UNIT, *) ToolName; !READ(UNIT, *) !FREQ     !tool frequency input in kHz
         READ(UNIT, *) NumFrequency
         READ(UNIT, *)FreqArray(1:NumFrequency) !tool frequency input in kHz
         CLOSE(UNIT)
      IF(TRIM(ToolName)=='EWR-M5'.OR.TRIM(ToolName)=='EWR-P4')then
          Itool=-4;
      Endif
      IF(TRIM(ToolName)=='ARC675'.OR.TRIM(ToolName)=='ARC825'
     1                           .OR.TRIM(ToolName)=='ARC900'.or.
     2   TRIM(ToolName)=='ARC312'.OR.TRIM(ToolName)=='ARC475'
     3                           .OR.TRIM(ToolName)=='IMPulse')THEN
          Itool=-5;
      ENDIF    
      IF(TRIM(ToolName)=='MPR'.or.TRIM(ToolName)=='MPR95')THEN
          Itool=-6; 
      ENDIF!   
         SELECT CASE (ITOOL)
         CASE(-4)
             !CHECK fREQUECNY NUMBER
             IF(TRIM(ToolName)=='EWR-M5')then
                IF(NumFrequency>3)THEN !EWR-M5
                 WRITE(*,*); WRITE(*,*)'Logging tool name: EWR-M5';
                 WRITE(*,*)'NumFrequency input error; '
                 WRITE(*,*)'Num OF Frequency sholud be <=3';
                 WRITE(*,*); STOP
                ENDIF
             ELSE IF(TRIM(ToolName)=='EWR-P4')then   
                IF(NumFrequency>2)THEN !EWR-P4
                 WRITE(*,*); WRITE(*,*)'Logging tool name: EWR-P4';   
                 WRITE(*,*); WRITE(*,*)'NumFrequency input error; '
                 WRITE(*,*)'Num OF Frequency sholud be <=2';
                 WRITE(*,*); STOP
                ENDIF
             ENDIF
             !CHECK fREQUECNY VALUES
             DO I=1,NumFrequency
              IF(TRIM(ToolName)=='EWR-M5')then   
                IF(DABS(FreqArray(I)-250.0)>1.0D-6.AND.
     1            DABS(FreqArray(I)-500.0)>1.0D-6.AND.
     2            DABS(FreqArray(I)-2000.0)>1.0D-6)THEN
                  WRITE(*,*); WRITE(*,*)'Logging tool name: EWR-M5';
                  WRITE(*,*); WRITE(*,*)'Frequency-VALUE input error; '
                  WRITE(*,*)'Frequency-VALUE sholud be 250, 500K, 2000K'
                  WRITE(*,*); STOP
                ENDIF
              ELSE IF(TRIM(ToolName)=='EWR-P4')then   
                IF(DABS(FreqArray(I)-1000.0)>1.0D-6.AND.
     2            DABS(FreqArray(I)-2000.0)>1.0D-6)THEN
                  WRITE(*,*); WRITE(*,*)'Logging tool name: EWR-P4';
                  WRITE(*,*); WRITE(*,*)'Frequency-VALUE input error; '
                  WRITE(*,*)'Frequency-VALUE sholud be 1000K, 2000K'
                  WRITE(*,*); STOP
                 ENDIF
              ENDIF
             ENDDO
         CASE(-5,-6) !ARC; MPR
             IF(NumFrequency>2)THEN
              WRITE(*,*); WRITE(*,*)'Logging tool name: ',TRIM(ToolName);
              WRITE(*,*); WRITE(*,*)'NumFrequency input error; '
              WRITE(*,*)'Num OF Frequency sholud be <=2'
              WRITE(*,*); STOP
             ENDIF
             !CHECK fREQUECNY VALUES
             DO I=1,NumFrequency
              IF(DABS(FreqArray(I)-400.0)>1.0D-6.AND.
     2           DABS(FreqArray(I)-2000.0)>1.0D-6)THEN
               WRITE(*,*);WRITE(*,*)'Logging tool name: ',TRIM(ToolName)
               WRITE(*,*);WRITE(*,*)'Frequency-VALUE input error; '
               WRITE(*,*)'Frequency-VALUE sholud be 400, 2000K'
               WRITE(*,*); STOP
              ENDIF
             ENDDO
         !CASE(-6) 
         END SELECT
         !CHECK NumFrequency, FreqArray(1:NumFrequency)
         CALL Output_PARTOOLinput(NumFrequency,FreqArray)
         RETURN
         !============================================================================
         !
	   FREQ=1000.0*FREQ      !kKz - > Hz
	   !IF(FREQ < 1.0D3) STOP 'Frequency too small!'
	   !IF(FREQ > 1.0D5) STOP 'Frequency too large!'
	   IF(FREQ < 1.0D3) write(*,*)'Frequency too small!'
	   IF(FREQ > 1.0D5) write(*,*)'Frequency too large!'
	   READ(UNIT, *) RADIUS 
	   I = 1
2        CONTINUE
	   READ(UNIT, *, END = 5) PP(I), ITR(I), NN(I)
	   I = I + 1
	   IF(I > MAXT + MAXR) STOP 'Too many coils!'
	   GOTO 2
5     CLOSE(UNIT)
C
	NOCL = I - 1
	NTT = 0
	NRR = 0
	DO 7 I = 1, NOCL
	   IF(ITR(I) == 1) THEN
	      NTT = NTT + 1
	      IF(NTT > MAXT) STOP 'Too many transmitters!'
	      PTT(NTT) = PP(I) * FINCH;  MMT(NTT) = NN(I)
	      GOTO 7
	   ELSEIF(ITR(I) == 0) THEN
	      NRR = NRR + 1
	      IF(NRR > MAXR) STOP 'Too many receivers!'
	      PRR(NRR) = PP(I) * FINCH; MMR(NRR) = NN(I)
	      GOTO 7
	   ELSE
	      STOP 'Coil flag has to be 0 or 1!'		  
	   END IF
7     CONTINUE
C
	OMEGA = 2.0D0 * PI * FREQ
	RADIUST = RADIUS * FINCH
	RADIUSR = RADIUS * FINCH
	S0T = PI * RADIUST * RADIUST
	S0R = PI * RADIUSR * RADIUSR
C
C-----------------------------------------------------------------------
C Compute tool constants.
C-----------------------------------------------------------------------
	TTK = DZERO
	T1 = - OMEGA ** 2 * U0 ** 2 * S0T * S0R * CURENT /(4.0D0 * PI)
	DO 10 IT = 1, NTT
	   DO 10 IR = 1, NRR
	      T2 = DABS(PTT(IT) - PRR(IR))
	      TLK(IT, IR) = T1 * MMT(IT) * MMR(IR) / T2
            TTK = TTK + TLK(IT, IR)
10    CONTINUE
C	
      END SUBROUTINE PARTOOL
!
!*******************************************************************************
!
      SUBROUTINE Output_PARTOOLinput(NumFrequency,FreqArray)
	USE GlobalVarbls,ONLY: UNIT0,ToolName
	IMPLICIT NONE
      !
      INTEGER :: NumFrequency; REAL(8):: FreqArray(10)
C====================================================================
!
	   WRITE(UNIT0,*)'INPUT FILE NAME: TOOL.INP'
         WRITE(UNIT0,*) TRIM(ToolName);
         WRITE(UNIT0,*) NumFrequency
         WRITE(UNIT0,*)FreqArray(1:NumFrequency) !tool frequency input in kHz
C	
      END SUBROUTINE Output_PARTOOLinput
!
!*******************************************************************************
!
      SUBROUTINE PARTOOL_glw   !org version
C-----------------------------------------------------------------------
C Subroutine to read in tool configration.
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
	REAL*8 :: PP(MAXT + MAXR), ITR(MAXT + MAXR), NN(MAXT + MAXR)
C
	REAL*8 FREQ, RADIUS , T1, T2
	INTEGER I, IT, IR
C
	OPEN(UNIT, FILE = 'TOOL.INP', STATUS = 'OLD')
	   READ(UNIT, *) FREQ     !tool frequency input in Hz
	   IF(FREQ < 1.0D3) STOP 'Frequency too small!'  !f>1kHz
	   IF(FREQ > 1.0D5) STOP 'Frequency too large!'  !f<100kHz
	   READ(UNIT, *) RADIUS 
	   I = 1
2        CONTINUE
	   READ(UNIT, *, END = 5) PP(I), ITR(I), NN(I)
	   I = I + 1
	   IF(I > MAXT + MAXR) STOP 'Too many coils!'
	   GOTO 2
5     CLOSE(UNIT)
C
	NOCL = I - 1
	NTT = 0
	NRR = 0
	DO 7 I = 1, NOCL
	   IF(ITR(I) == 1) THEN
	      NTT = NTT + 1
	      IF(NTT > MAXT) STOP 'Too many transmitters!'
	      PTT(NTT) = PP(I) * FINCH;  MMT(NTT) = NN(I)
	      GOTO 7
	   ELSEIF(ITR(I) == 0) THEN
	      NRR = NRR + 1
	      IF(NRR > MAXR) STOP 'Too many receivers!'
	      PRR(NRR) = PP(I) * FINCH; MMR(NRR) = NN(I)
	      GOTO 7
	   ELSE
	      STOP 'Coil flag has to be 0 or 1!'		  
	   END IF
7     CONTINUE
C
	OMEGA = 2.0D0 * PI * FREQ
	RADIUST = RADIUS * FINCH
	RADIUSR = RADIUS * FINCH
	S0T = PI * RADIUST * RADIUST
	S0R = PI * RADIUSR * RADIUSR
C
C-----------------------------------------------------------------------
C Compute tool constants.
C-----------------------------------------------------------------------
	TTK = DZERO
	T1 = - OMEGA ** 2 * U0 ** 2 * S0T * S0R * CURENT /(4.0D0 * PI)
	DO 10 IT = 1, NTT
	   DO 10 IR = 1, NRR
	      T2 = DABS(PTT(IT) - PRR(IR))
	      TLK(IT, IR) = T1 * MMT(IT) * MMR(IR) / T2
            TTK = TTK + TLK(IT, IR)
10    CONTINUE
C	
      END
!
!********************************************************************************************
!    The End
!********************************************************************************************
!