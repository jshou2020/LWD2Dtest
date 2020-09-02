!
!********************************************************************************************
!    File name: INDMod2DMedTBF.for
!    DATE         : Wed, 09/02/2020
!    LAST UPDATE  : Wed, 09/02/2020
!********************************************************************************************
!
!********************************************************************************************
!    SUBROUTINE INDMod2DMedTBF(NumLyrs, NumCldLyrs, SIG0,VECOD0,HOR0,...,  logDEP)
!********************************************************************************************
!
	SUBROUTINE INDMod2DMedTBF(NumLyrs, NumCldLyrs, SIG0,VECOD0,HOR0, 
     1                          SAPOINT0,DELTA0,IMAX0, FREQ, 
     2                    AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDEP)
C-----------------------------------------------------------------------
!
!Input parameters (formation model, logging setting, & frequency) of SUBROUTINE INDMod2DMedTBF()
!    NumLyrs,   integer,  the number of formation beds;
!    NumCldLyrs, integer,  the number of radial layers in each bed 
!    SIG0(i,j),  REAL*8, DIMENSION(0 : NumLyrs-1, 0 : NumCldLyrs-1) 
!                the resistivity of radial section j in i-th bed
!    VECOD0(i,j),  REAL*8, DIMENSION(0 : NumLyrs-1, 0 : NumCldLyrs-1)
!                the radiu of radial section j in i-th bed (unit=inch)
!    HORO(i), REAL*8, DIMENSION(NumLyrs-1),
!             the position of i-th formation bounary, also the 
C               position of lower boundary of layer i. Here, by lower 
C			  we mean the boundary whose coordinate is smaller. In 
C               the present context, the z-axis is pointed downward, 
C               so the lower boundary is physically above the upper 
C               boundary 
!  SAPOINT0,   REAL(8),   the position (unit: feet) of first logging point
!  DELTA0,  REAL(8),    the spacing (unit: inches) between two consecutive
!                      logging points 
!  IMAX0:   total logging point number
!  FREQ: REAL(8), tool's the operating frequency (unit = kHz)
!      
!Output parameters (phase, attenuation & their apparent resistivities) of SUBROUTINE INDMod2DMedTBF()   
!   AllPhase(i,j), REAL(8), DIMENSION(IMAX0,6), (j=1-5 from short to long arrays)
!                  phase of j-th array at i-th log point (unit =deg);
!   AllAtt(i,j), REAL(8), DIMENSION(IMAX0,6), (j=1-5 from short to long arrays)
!                  attenuation of j-th array at i-th log point (unit=db);
!   RaAllPhase(i,j), REAL(8), DIMENSION(IMAX0,6), (j=1-5 from short to long arrays)
!                  phase resistivity of j-th array at i-th log point (unit = ohm-m);
!   RaAllAt(i,j), REAL(8), DIMENSION(IMAX0,6), (j=1-5 from short to long arrays)
!                  attnuation resistivity of j-th array at i-th log point (unit = ohm-m);
!    logDEP(i):   REAL(8), DIMENSION(IMAX0), log depth at i-th log point (unit = ft).
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
      INTEGER :: NumLyrs, NumCldLyrs;   INTEGER:: NBOU0, MAXS0
	REAL*8, DIMENSION(0 : NumLyrs-1, 0 : NumCldLyrs-1) :: SIG0, VECOD0
	REAL*8, DIMENSION(NumLyrs-1) :: HOR0
      REAL(8):: SAPOINT0,DELTA0; INTEGER:: IMAX0
      !
	!integer:: Itool;  
      REAL*8 T1, T2
	REAL(8):: AllPhase(IMAX0,6),AllAtt(IMAX0,6);
	REAL(8):: RaAllPhase(IMAX0,6),RaAllAtt(IMAX0,6), logDep(IMAX0);
      REAL(8):: FREQ
C---------------------------------------------------------------------------
!
	!CALL COPYRIGHT; CALL CPU_TIME(T1);	PRINT*,'MODELING BEGINS!';	PRINT*
C
	!CALL INDAT  
      NBOU0=NumLyrs-1; MAXS0=NumCldLyrs-1
      CALL INDAT(NBOU0, MAXS0, SIG0,VECOD0,HOR0, SAPOINT0,DELTA0,IMAX0)
	CALL AllocGlbVecsMats
C
	!WRITE(*,'(1X,25HNUMBER OF LOGGING POINTS:,I8)') IMAX;	PRINT*
C
      CALL PARTOOL_02(FREQ); !CALL PARTOOL_HAL(FREQ); 
	!IF(Itool==-4)CALL PARTOOL_HAL(FREQ); !IF(Itool==-5)CALL PARTOOL_ARC(FREQ)
!      
	CALL GRID; 	CALL PREDT
C
      !LWD-OD look-up table for mapping from ph att to apparent Ra;!Itool =-4; 
      CALL Build_LWDODtable(Itool); !STOP
c	
	CALL EIGENCAL;	CALL QSTRAN; 
      !
      IF(Itool==-4)THEN
       CALL SOLVE_EWR(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep);  !EWR-M5; EWR-P4;
      ELSEIF(Itool==-5)THEN
       CALL SOLVE_ARC(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep);  !ARC
       !                                                              'ARC675,ARC825,ARC900',
       !                                                              'ARC312,ARC475','IMPulse';   
      ELSEIF(Itool==-6)THEN
       CALL SOLVE_MPR(AllPhase,AllAtt, RaAllPhase,RaAllAtt, logDep);  !MPR; MPR95
      ENDIF
      !
      !PRINT*;	PRINT*,'MODELING ENDS!';	PRINT*
	!CALL CPU_TIME(T2);	WRITE(*,'(1X,18HTOTAL CPU TIME(S):,F10.2)') T2-T1
	!CALL COPYRIGHT; CALL CleanUp
C
      END SUBROUTINE INDMod2DMedTBF    
!
!********************************************************************************************
!     The End
!********************************************************************************************
!