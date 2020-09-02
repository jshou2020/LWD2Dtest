	SUBROUTINE AllocGlbVecsMats
C-----------------------------------------------------------------------
C Subroutine to allocate global vectors and matrices
C Called by: INDMod2DMedTBF
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
      IF (ALLOCATED(RM)) THEN; DEALLOCATE(RM,K2,IN, C,V); 
      ENDIF
	ALLOCATE(RM(0 : NBOU, 0 : N, MAXI))
	ALLOCATE(K2(0 : NBOU, 0 : N, 0 : MAXI))
	ALLOCATE(IN(0 : NBOU, 0 : N))
C
	ALLOCATE(C(0 : NBOU, N, N), V(0 : NBOU, N))
C
	IF(NBOU > 0) THEN
         IF (ALLOCATED(QU)) THEN; DEALLOCATE(QU,QD,SU,SD,P); ENDIF 
	   ALLOCATE(QU(NBOU, N, N), QD(NBOU, N, N))
	   ALLOCATE(SU(NBOU, N, N), SD(NBOU, N, N))
	   ALLOCATE(P(NBOU, N, N))
	   IF(NBOU > 1) THEN
             IF (ALLOCATED(E)) DEALLOCATE(E); ALLOCATE(E(NBOU - 1, N))
         ENDIF
      END IF
C
	IF(NBOU > 0) THEN
          IF (ALLOCATED(U)) DEALLOCATE(U); ALLOCATE(U(NBOU, N))
      ENDIF
C
      IF (ALLOCATED(VR)) DEALLOCATE(VR); ALLOCATE(VR(IMAX))
C
	END SUBROUTINE AllocGlbVecsMats
