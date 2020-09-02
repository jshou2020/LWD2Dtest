!
!****************************************************************************************
!   File name: CalcPhaAttFromVoltage.f90
!   Routine  : ComplexV_2PhaseAtt (C, phase,Att)
!   Purpose  : Calculation of  Phase & Att 
!   DATE         : Wed, 07/15/2020
!   LAST UPDATE  : Mon, 08/10/2020
!   ONLY CALLED BY 
!****************************************************************************************
!   SUBROUTINE ComplexV_2PhaseAtt (C, phase,Att)
    ! Purpose:
    ! Accept a complex number C = RE + i*IM, return its Attenuation & phase "phase" 
    ! This sub returns the Phase in deg.
    ! Record of revisions:
    ! Date Programmer Description of change
    ! ==== ========== =====================
    ! 01/10/97 S. J. Chapman Original code
!****************************************************************************************
!    
    SUBROUTINE ComplexV_2PhaseAtt (C, phase,Att)
    IMPLICIT NONE
    REAL(8), PARAMETER:: PI  = 3.141592653589793238  !Pi constant
    REAL(8), PARAMETER:: radian2deg = 180.0/PI       !radian to deg 
    ! List of dummy arguments:
    INTEGER :: N,j
    complex(8),INTENT(IN) :: C      ! Input complex number
    REAL(8),INTENT(OUT)   :: Phase  ! Phase in deg
    REAL(8),INTENT(OUT)   :: Att    ! Att in db
!    
    ! Get att and phase for a complex number.
    !---------------------------------------------------------------
    !! z= x+iy
    !do j=1,N
     !Amp(j) = max(ABS (C(j)), 1.0D-16)
     !amp(j) = Dmax1(ABS ( c(j) ),1.0D-16)
     !phase(j) = ATAN2 (AIMAG(c(j)), REAL(c(j)))
     !phase(j) = ATAN2 (IMAG(C(j)), real(C(j)))
     !Phase(j) = ATAN2 (IMAG(C(j)), real(C(j)))      !tg-1(imag/real)
     !phase(j) = DATAN2D (IMAG(C(j)), real(C(j)))
    !enddo
     Att = max(ABS (C), 1.0D-16); Att=-20.D0*log10(Att);
     Phase = ATAN2(IMAG(C), real(C)) 
     Phase = Phase*radian2deg; !radian2deg=180.D0/pi;    
    !RETURN 
    END SUBROUTINE ComplexV_2PhaseAtt
!
!*****************************************************************************
!     SUBROUTINE Compute_3coilPhaseAtt(Vecc12,Phase,Att)
!     Purpose: complex voltage -> phase and att (T--- R1 - R2)
!     only called by: 
!*****************************************************************************
!
      SUBROUTINE Compute_3coilPhaseAtt(Vecc12,Phase,Att)
      USE GlobalVarbls,only: Itool, toolname; !USE Tool_pars,only: NARRAYzz
      USE Tool5_par,ONLY: NARRAYzz5   ! EWR-M5, Itool=-4
      USE Tool52_par,ONLY: NARRAYzz52 ! EWR-P4, Itool=-4
      USE Tool6_par,ONLY: NARRAYzz6   ! ARC,    Itool=-5
      USE Tool7_par,ONLY: NARRAYzz7   ! MPR,    Itool=-6
      !
      IMPLICIT NONE
	  complex*16,INTENT(IN):: Vecc12(12)        !in - complex v 
	  REAL(8),INTENT(OUT)::   Phase(6),Att(6);  !out 
! Locals:
      INTEGER:: Iary, I1,I2;  complex*16:: V; INTEGER:: NARRAYzz;
!=======================
!
      IF(Itool==-4.and.TRIM(toolname)=='EWR-M5')NARRAYzz=NARRAYzz5; 
      IF(Itool==-4.and.TRIM(toolname)=='EWR-P4')NARRAYzz=NARRAYzz52; 
      !
      IF(Itool==-5)NARRAYzz=NARRAYzz6
      IF(Itool==-6)NARRAYzz=NARRAYzz7; 
      DO Iary=1, NARRAYzz;    
        IF(Itool==-6.OR.(Itool==-4.and.TRIM(toolname)=='EWR-P4'))then
          I1=2*Iary-1; I2=2*Iary  
        !else IF(Itool==-4.or.Itool==-5)then
        else IF((Itool==-4.and.TRIM(toolname)=='EWR-M5').or.Itool==-5)then
          I1=Iary; I2=I1+1;      
        endif
        V=Vecc12(I1)/Vecc12(I2); !V=Vnear/Vfar
        CALL ComplexV_2PhaseAtt(V, phase(Iary),Att(Iary))
      ENDDO
      !RETURN     
      END SUBROUTINE Compute_3coilPhaseAtt
!
!****************************************************************************************
!   THE END
!****************************************************************************************
!    
!****************************************************************************************
!
!   SUBROUTINE Complex_AmpPhase3D (N, c, amp, phase)
!****************************************************************************************
!    
    SUBROUTINE Complex_AmpPhase3D (N, C, amp, phase)
    ! Purpose:
    ! Accept a complex number C = RE + i*IM, return its amplitude "amp" & phase "phase" 
    ! This sub returns the phase in radians.
    !
    ! Record of revisions:
    ! Date Programmer Description of change
    ! ==== ========== =====================
    ! 01/10/97 S. J. Chapman Original code
    !
    IMPLICIT NONE
    ! List of dummy arguments:
    INTEGER :: N,j
    complex(8),INTENT(IN) :: C(N)      ! Input complex number
    REAL(8),INTENT(OUT)   :: Amp(N)    ! Amplitude
    REAL(8),INTENT(OUT)   :: Phase(N)  ! Phase in radians
    ! Get amplitude and phase for a complex number.
    !---------------------------------------------------------------
    !! z= x+iy
    do j=1,N
     Amp(j) = max(ABS (C(j)), 1.0D-16)
     !amp(j) = Dmax1(ABS ( c(j) ),1.0D-16)
     !phase(j) = ATAN2 (AIMAG(c(j)), REAL(c(j)))
     !phase(j) = ATAN2 (IMAG(C(j)), real(C(j)))
     Phase(j) = ATAN2 (IMAG(C(j)), real(C(j)))      !tg-1(imag/real)
     !phase(j) = DATAN2D (IMAG(C(j)), real(C(j)))
    enddo
    !RETURN 
    END SUBROUTINE Complex_AmpPhase3D