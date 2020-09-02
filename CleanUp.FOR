	SUBROUTINE CleanUp
C-----------------------------------------------------------------------
C Subroutine to clean up the memory used by the program.
C Called by: INDMod2DMedTBF
C-----------------------------------------------------------------------
	USE GlobalVarbls
      IMPLICIT NONE
C
	DEALLOCATE(EPS, SIG, VECOD)
	IF(NBOU > 0) DEALLOCATE(HOR)
	DEALLOCATE(RM, K2, IN, C, V)
	IF(NBOU > 0) DEALLOCATE(QU, QD, SU, SD, P)
	IF(NBOU > 1) DEALLOCATE(E)
	IF(NBOU > 0) DEALLOCATE(U)
C
	END SUBROUTINE CleanUp
