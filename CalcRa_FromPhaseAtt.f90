!
!****************************************************************************************
!   File name: Calc Ra from phase and att.f90
!   Routine  : SUBROUTINE PhaseAtt_2Ra(Ifre, Phase,Att, RaPhase,RaAtt)
!   Purpose  : Calc Ra from phase and att
!   DATE         : Thu, 07/23/2020
!   LAST UPDATE  : Mon, 08/10/2020
!   ONLY CALLED BY 
!****************************************************************************************
!
!   SUBROUTINE PhaseAtt_2Ra(Ifre, Phase,Att, RaPhase,RaAtt)
!****************************************************************************************
!
    SUBROUTINE PhaseAtt_2Ra(Ifre, Phase,Att, RaPhase,RaAtt)
    USE LWDODtable;  IMPLICIT NONE
    INTEGER:: Ifre
    REAL(8):: Phase(6),Att(6), RaPhase(6),RaAtt(6);
    !locals:
    INTEGER:: Isub, I1,I2, IRh, Nsub; REAL(8):: W1,W2,W12
    !=====================================================
    !
    Nsub= size(PhaseOD,2)
    !RaPhase(1:Nsub)
    DO Isub = 1, Nsub  !SUB loop
      I1=0; I2=0
      DO IRh=1, NRHOD-1
       IF(Phase(Isub)>=PhaseOD(Ifre,Isub,IRh).AND.Phase(Isub)<=PhaseOD(Ifre,Isub,IRh+1))THEN
        I1=IRh; I2=I1+1; 
        W1=Phase(Isub)-PhaseOD(Ifre,Isub,IRh); W2=PhaseOD(Ifre,Isub,IRh+1)-Phase(Isub)
        W12=PhaseOD(Ifre,Isub,IRh+1)-PhaseOD(Ifre,Isub,IRh)
        W1=W1/W12; W2=W2/W12
        RaPhase(Isub)=W2*RhOD(I1)+W1*RhOD(I2); !RaPhase(Isub)=0.5*(RhOD(I1)+RhOD(I2))
       ENDIF
      ENDDO !Rh loop
      IF(I1==0.AND.I2==0)THEN
       IF(Phase(Isub)<=PhaseOD(Ifre,Isub,1))THEN; RaPhase(Isub)=RhOD(1); ENDIF
       IF(Phase(Isub)>=PhaseOD(Ifre,Isub,NRHOD))THEN; RaPhase(Isub)=RhOD(NRHOD); ENDIF
      ENDIF
    ENDDO !SUB loop
    !
    !RaAtt(1:Nsub)
    DO Isub = 1, Nsub  !SUB loop
      I1=0; I2=0
      DO IRh=1, NRHOD-1
       IF(Att(Isub)>=AttOD(Ifre,Isub,IRh).AND.Att(Isub)<=AttOD(Ifre,Isub,IRh+1))THEN
        I1=IRh; I2=I1+1; 
        W1=Att(Isub)-AttOD(Ifre,Isub,IRh); W2=AttOD(Ifre,Isub,IRh+1)-Att(Isub)
        W12=AttOD(Ifre,Isub,IRh+1)-AttOD(Ifre,Isub,IRh)
        W1=W1/W12; W2=W2/W12
        RaAtt(Isub)=W2*RhOD(I1)+W1*RhOD(I2); !RaAtt(Isub)=0.5*(RhOD(I1)+RhOD(I2))
       ENDIF
      ENDDO !Rh loop
      IF(I1==0.AND.I2==0)THEN
       IF(Att(Isub)<=AttOD(Ifre,Isub,1))THEN; RaAtt(Isub)=RhOD(1); ENDIF
       IF(Att(Isub)>=AttOD(Ifre,Isub,NRHOD))THEN; RaAtt(Isub)=RhOD(NRHOD); ENDIF
      ENDIF
    ENDDO !SUB loop
      !
    END SUBROUTINE PhaseAtt_2Ra
!
!****************************************************************************************
!   THE END
!****************************************************************************************
!  