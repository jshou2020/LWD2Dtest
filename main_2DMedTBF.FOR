!
!********************************************************************************************
!     Filename: main_2DMedTBF.for  
!     DATE         : Wed, 09/02/2020
!     LAST UPDATE  : Wed, 09/02/2020
!********************************************************************************************
!
!     PROGRAM main_2DMedTBF
!********************************************************************************************
!
	PROGRAM main_2DMedTBF
C-----------------------------------------------------------------------
	USE GlobalVarbls, ONLY: UNIT0, UNIT, Itool, ToolName
	IMPLICIT NONE
	!integer:: Itool;  
      REAL*8 T1, T2; REAL(8):: FREQ
	!REAL(8):: AllPhase(300,6),AllAtt(300,6);!REAL(8):: RaAllPhase(300,6),RaAllAtt(300,6)
      REAL*8, DIMENSION(:, :), ALLOCATABLE :: AllPhase,AllAtt;
      REAL*8, DIMENSION(:, :), ALLOCATABLE :: RaAllPhase,RaAllAtt
C----------------------------------------------------------------------------------
	REAL*8, DIMENSION(:, :), ALLOCATABLE :: EPS, SIG, VECOD
	REAL*8, DIMENSION(:), ALLOCATABLE :: HOR, logDep
      INTEGER :: NBOU, MAXS, MAXI
      INTEGER :: NumLyrs, NumCldLyrs
      integer:: IMAX
	REAL*8 :: SAPOINT, DELTA;  !CHARACTER*10:: ToolName
      INTEGER :: NumFrequency; REAL(8):: FreqArray(10); INTEGER :: IIF
!================================================================================
      !
	!CALL COPYRIGHT; 
      !
      OPEN(UNIT0, FILE='LWD2DRTsimulation.out', STATUS='replace') 
      CALL CPU_TIME(T1);	PRINT*,'MODELING BEGINS!';	PRINT*
C
	OPEN(UNIT, FILE = 'IND.INP', STATUS='OLD') !model pars input file
	READ(UNIT, *) NumLyrs, NumCldLyrs  !!!
	NBOU = NumLyrs - 1
	MAXS = NumCldLyrs - 1;  MAXI = MAXS
	ALLOCATE(EPS(0 : NBOU, 0 : MAXS), SIG(0 : NBOU, 0 : MAXS))
	ALLOCATE(VECOD(0 : NBOU, 0 : MAXS));  ALLOCATE(HOR(NBOU));
	!CALL INDAT 
	CALL INDAT0(NBOU, MAXS, SIG,VECOD,HOR, SAPOINT,DELTA,IMAX)
	!CALL AllocGlbVecsMats
C
	WRITE(*,'(1X,25HNUMBER OF LOGGING POINTS:,I8)') IMAX;	PRINT*;
      !WRITE(*,*)NBOU, MAXS; !WRITE(*,*) SAPOINT/0.3048,DELTA/0.0254D0
      !PAUSE
      ALLOCATE(AllPhase(IMAX,6),AllAtt(IMAX,6)); 
      ALLOCATE(RaAllPhase(IMAX,6),RaAllAtt(IMAX,6), logDep(IMAX))
C
	CALL PARTOOL(NumFrequency,FreqArray); !CALL PARTOOL(FREQ); 
    ! IF(TRIM(ToolName)=='EWR-M5')then
    !     Itool=-4;
    ! Endif
    ! IF(TRIM(ToolName)=='ARC675'.OR.TRIM(ToolName)=='ARC825'
    !1                           .OR.TRIM(ToolName)=='ARC900'.or.
    !2   TRIM(ToolName)=='ARC312'.OR.TRIM(ToolName)=='ARC475'
    !3                           .OR.TRIM(ToolName)=='IMPulse')THEN
    !     Itool=-5;
    ! ENDIF    
    ! IF(TRIM(ToolName)=='MPR'.or.TRIM(ToolName)=='MPR95')THEN
    !     Itool=-6; 
    ! ENDIF!   
      DO IIF =1, NumFrequency
         FREQ= FreqArray(IIF); !write(*,*)trim(ToolName),FREQ,Itool;pause
       !CALL INDMod2DMedTBF(FREQ, AllPhase,AllAtt, RaAllPhase,RaAllAtt)
       !CALL INDMod2DMedTBF(NBOU, MAXS, SIG,VECOD,HOR, SAPOINT,DELTA,IMAX,
       CALL INDMod2DMedTBF(NumLyrs, NumCldLyrs, SIG,VECOD,HOR, 
     1                    SAPOINT,DELTA,IMAX, FREQ, 
     !                    AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDEP)
	 !CALL GRID; CALL PREDT
C
       !LWD-OD look-up table for mapping from ph att to apparent Ra
       !Itool =-4; CALL Build_LWDODtable(Itool); !STOP
c	
	 !CALL EIGENCAL;	CALL QSTRAN; !CALL SOLVE;
       CALL Output_calcPA_RaPA(AllPhase,AllAtt, RaAllPhase,RaAllAtt,   
     1                        logDep)
      ENDDO
      PRINT*; PRINT*,'MODELING ENDS!';	PRINT*
	CALL CPU_TIME(T2);	WRITE(*,'(1X,18HTOTAL CPU TIME(S):,F10.2)') T2-T1; !CALL COPYRIGHT; 
      CALL CleanUp
      CLOSE(UNIT0)
      END PROGRAM main_2DMedTBF
!
!********************************************************************************************
!     SUBROUTINE Output_calcPA_RaPA(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep)
!********************************************************************************************
	SUBROUTINE Output_calcPA_RaPA(AllPhase,AllAtt, RaAllPhase,RaAllAtt,
     1                              logDep)
C-----------------------------------------------------------------------
	USE GlobalVarbls; use Tool5_par,only: NARRAYzz5; 
      use Tool52_par,only: NARRAYzz52;
      use Tool6_par,only: NARRAYzz6;  use Tool7_par,only: NARRAYzz7
	IMPLICIT NONE
C
	!COMPLEX*16 :: RE, VRTEMP, Rearray(12,12)
	INTEGER :: I, NARRAYzz ! JT, JR, NT, NR, IFLAG, NSpacing, 
	REAL*8 :: ZS, ZT, ZR
      REAL(8):: AllPhase(IMAX,6),AllAtt(IMAX,6);
	REAL(8):: RaAllPhase(IMAX,6),RaAllAtt(IMAX,6), logDep(IMAX);
	!INTEGER:: IFREQ
	REAL(8):: Freq
C===========================================================================
!     
      !Freq=OMEGA/(2*PI)/1000.0; 
      IF(Itool==-4.AND.TRIM(toolname)=='EWR-M5')THEN
          NARRAYzz=NARRAYzz5; 
      ELSE IF(Itool==-4.AND.TRIM(toolname)=='EWR-P4')THEN
          NARRAYzz=NARRAYzz52; 
      ENDIF
      !
      IF(Itool==-5)NARRAYzz=NARRAYzz6; IF(Itool==-6)NARRAYzz=NARRAYzz7;
      !RETURN
!	OPEN(UNIT,FILE = 'IND.OUT')
!	   DO 7 I = 1, IMAX
!	      !ZS = SAPOINT + (I - 1) * DELTA;   ZS = ZS / FFOOT
!            ZS= logDep(I);
!            WRITE(*, 10) ZS, VR(I);   WRITE(UNIT, 10) ZS, VR(I)
10          FORMAT(1X, F16.7, 20F16.7)
!7	  CONTINUE
!	  CLOSE(UNIT)
C
      Freq=OMEGA/(2*PI)/1000
!      
	IF(ABS(Freq-2000.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF2000.OUT')
	IF(ABS(Freq-500.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF500.OUT')
	IF(ABS(Freq-250.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF250.OUT')
	IF(ABS(Freq-400.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF400.OUT')
	IF(ABS(Freq-1000.0)<1.0D-8)OPEN(UNIT,FILE ='PhasAttF1000.OUT')
!      
	   DO 77 I = 1, IMAX
	      !ZS = SAPOINT + (I - 1) * DELTA;   ZS = ZS / FFOOT
            ZS= logDep(I);
            !WRITE(*, 10) ZS, VR(I);  !WRITE(UNIT, 10) ZS, VR(I)
            !WRITE(*, 10) ZS, AllPhase(I,1:NARRAYzz),
     1      !                 AllAtt(I,1:NARRAYzz) ; 
            WRITE(UNIT, 110)ZS, AllPhase(I,1:NARRAYzz),
     1                          AllAtt(I,1:NARRAYzz) ;
110          FORMAT(1X, F16.7, 20F16.7)
77	  CONTINUE
	  CLOSE(UNIT)
!
	IF(ABS(Freq-2000.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF2000.OUT')
	IF(ABS(Freq-500.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF500.OUT')
	IF(ABS(Freq-250.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF250.OUT')
	IF(ABS(Freq-400.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF400.OUT')
	IF(ABS(Freq-1000.0)<1.0D-8)OPEN(UNIT,FILE = 'RaPhasAttF1000.OUT')
!      
	   DO I = 1, IMAX
	      !ZS = SAPOINT + (I - 1) * DELTA;  ZS = ZS / FFOOT
            ZS= logDep(I);
            !WRITE(*, 10) ZS, RaAllPhase(I,1:NARRAYzz), RaAllAtt(I,1:NARRAYzz) ; 
    !!       IF(TRIM(toolname)=='EWR-P4'.AND.ABS(Freq-2000.0)<1.0D-8)THEN
    !!        WRITE(UNIT, 110)ZS, RaAllPhase(I,1:NARRAYzz-1),
    !!1                          RaAllAtt(I,1:NARRAYzz-1) ; 
    !!       ELSE IF(TRIM(toolname)=='EWR-P4'.AND.
    !!1                                  ABS(Freq-1000.)<1.0D-8)THEN
    !!        WRITE(UNIT, 110)ZS, RaAllPhase(I,NARRAYzz),
    !!1                          RaAllAtt(I,NARRAYzz) ; 
!            ELSE
             WRITE(UNIT, 110)ZS, RaAllPhase(I,1:NARRAYzz),
     1                           RaAllAtt(I,1:NARRAYzz) ; 
            !ENDIF
	  ENDDO
	  CLOSE(UNIT)
        END SUBROUTINE Output_calcPA_RaPA
!
!********************************************************************************************
!    The End
!********************************************************************************************
!