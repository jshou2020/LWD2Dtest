!
!****************************************************************************************
!   File name: Build_LWDODtable.f90
!   Routine  : MODULE LWDODtable, Build_LWDODtable(Itool) 
!   Purpose  : LWDODtable
!   DATE         : Wed, 07/22/2020
!   LAST UPDATE  : Wed, 09/02/2020
!   ONLY CALLED BY 
!****************************************************************************************
!
!   MODULE LWDODtable
!****************************************************************************************
!
    MODULE LWDODtable
    IMPLICIT NONE
    INTEGER:: NRHOD
    REAL(8),ALLOCATABLE:: RhOD(:); !e.g., RhOD(NRHOD);
    REAL(8),ALLOCATABLE:: PhaseOD(:,:,:); !e.g., PhaseOD(NF,NSUB,NRHOD); 
    REAL(8),ALLOCATABLE:: AttOD(:,:,:);   !e.g., AttOD(NF,NSUB,NRHOD);
    !
    END MODULE LWDODtable
!
!****************************************************************************************
!   SUBROUTINE Build_LWDODtable(Itool)
!****************************************************************************************
!    
    SUBROUTINE Build_LWDODtable(Itool)
    USE GlobalVarbls, ONLY: toolname
    USE LWDODtable; USE Tool5_par ! EWR-M5-(Itool==-4)
    USE Tool52_par ! EWR-P4-(Itool==-4)
    USE Tool6_par;  !ARC675-(Itool==-5)
    USE Tool7_par;  !MPR   -(Itool==-6)
    IMPLICIT NONE
    INTEGER:: Itool
    INTEGER:: IRT, IFRE, IReceiver, NF,NSUB , NReceiver
    REAL(8):: DIP,FREQ,Z0,RH,RV; COMPLEX(8):: HZZ(12); REAL(8):: PhaseTemp(6),AttTemp(6);
    REAL(8),allocatable:: Rspacing(:),FreqArray(:)
!==============================================================
    !
    NRHOD=68; !58; !45,42; !34; !32
    IF (ALLOCATED(RhOD)) DEALLOCATE(RhOD); ALLOCATE(RhOD(NRHOD));
    !
    RhOD=(/0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.,1.5,2.,3.,4.,5.,6.,7.,8.,9.,10., &
           10.25,10.5,11.,12.,13.,14.,15.,16.,17.,17.5,18.,18.5,19.,19.5,&
           20.,21.,22.,23.,24.,25.0,&
           26.,27.,28.,29.,30.,31.,32., 34., 36.,38., 40.,42., 44., 46., 48., &
           50.,52.,54.,58.,60.,70.,80.,90., &
           100.,150.,200./)
    IF(Itool==-4)THEN  !EWR-M5 OR EWR-P4
      IF(TRIM(toolname)=='EWR-M5')THEN
       NF=NFREQ35; NSUB=NARRAYzz5; NReceiver=NReceiver5;
       allocate(Rspacing(NReceiver),FreqArray(NF)); Rspacing=Rspacing5; FreqArray=FREQ35;
      ELSE IF(TRIM(toolname)=='EWR-P4')THEN
       NF=NFREQ352; NSUB=NARRAYzz52; NReceiver=NReceiver52; 
       allocate(Rspacing(NReceiver),FreqArray(NF)); Rspacing=Rspacing52; FreqArray=FREQ352;
      ENDIF
      !allocate(Rspacing(NReceiver),FreqArray(NF)); Rspacing=Rspacing5; FreqArray=FREQ35;
    ELSE IF(Itool==-5)THEN !'ARC675,ARC825,ARC900', 'ARC312,ARC475','IMPulse';  
      NF=NFREQ36; NSUB=NARRAYzz6; NReceiver=NReceiver6;
      allocate(Rspacing(NReceiver),FreqArray(NF)); Rspacing=Rspacing6; FreqArray=FREQ36;
    ELSE IF(Itool==-6)THEN !MPR, MPR95
      NF=NFREQ37; NSUB=NARRAYzz7; NReceiver=NReceiver7;
      allocate(Rspacing(NReceiver),FreqArray(NF)); Rspacing=Rspacing7; FreqArray=FREQ37;
    ENDIF
    !
    IF(ALLOCATED(PhaseOD)) DEALLOCATE(PhaseOD, AttOD); 
    ALLOCATE(PhaseOD(NF,NSUB,NRHOD), AttOD(NF,NSUB,NRHOD));!WRITE(*,*)SIZE(PhaseOD,3);PAUSE
    !
    DO IFRE=1, NF
    DO IRT=1, NRHOD
     DO IReceiver=1,NReceiver
      FREQ=FreqArray(IFRE); !FREQ=FREQ35(IFRE); 
      Z0=Rspacing(IReceiver); !Z0=Rspacing5(IReceiver)
      dip = 0.0; RH=RhOD(IRT); RV=RH
      CALL HOMO_2COIL02(DIP,FREQ,Z0,RH,RV,HZZ(IReceiver))
     ENDDO !Receiver loop
     CALL Compute_3coilPhaseAtt(HZZ,PhaseTemp,AttTemp);
     !
     PhaseOD(IFRE,1:NSUB,IRT)=PhaseTemp(1:NSUB); AttOD(IFRE,1:NSUB,IRT)=AttTemp(1:NSUB)
    ENDDO !Rt loop
    ENDDO !Frequency loop
    deallocate(Rspacing,FreqArray); 
    RETURN
    !=============================================================================
    !
    DO IFRE=1, NF
         IF(ITOOL==-4.and.IFRE==1)OPEN(100, FILE='ODF250_CZZ.out')
         IF(ITOOL==-4.and.IFRE==2)OPEN(100, FILE='ODF500_CZZ.out')
         IF(ITOOL==-4.and.IFRE==3)OPEN(100, FILE='ODF2000_CZZ.out')
         !
         write(100,*)'% frequency = ',FreqArray(IFRE)
         write(100,*)'%  Rh, PhaseTemp(1:NSUB), AttTemp(1:NSUB)'
         DO IRT=1, NRHOD
          write(100,500)RhOD(IRT),PhaseOD(IFRE,1:NSUB,IRT),AttOD(IFRE,1:NSUB,IRT)
         enddo
    ENDDO
500 format(G14.4,1X, <2*NSUB>G17.7)
!500 format(G14.4,1X,6G17.7)
    !stop !
    END SUBROUTINE Build_LWDODtable
!
!****************************************************************************************
!   THE END
!****************************************************************************************
!    