!
!********************************************************************************************
!     Filename: INDAT.for  
!     DATE         : Wed, 09/02/2020
!     LAST UPDATE  : Wed, 09/02/2020
!********************************************************************************************
      !SUBROUTINE INDAT0
C-----------------------------------------------------------------------
C Subroutine to read in the model data.
C   NBOU      : the actual number of formation layers
C   EPS(i,j)  : the dielectric of section j in i-th layer
C   SIG(i,j)  : the condivity of section j in i-th layer
C   HOR(i)    : the position of i-th formation bounary, also the 
C               position of lower boundary of layer i. Here, by lower 
C			  we mean the boundary whose coordinate is smaller. In 
C               the present context, the z-axis is pointed downward, 
C               so the lower boundary is physically above the upper 
C               boundary
C   VECOD(i,j): the leftend of section j in i-th layer
C-----------------------------------------------------------------------
C								   |	
C						layer 0    |
C					---------------|---------------  z1
C						layer 1	   |
C					---------------|---------------  z2
C								   .
C								   .
C								   .
C					---------------|---------------  zNBOU-1
C					  layer NBOU-1 |
C					---------------|---------------  zNBOU
C								   |
C					   layer NBOU  |
C								   V										
C-----------------------------------------------------------------------
      !
      SUBROUTINE INDAT0(NBOU,MAXS, SIG,VECOD,HOR, SAPOINT,DELTA,IMAX)
	USE GlobalVarbls, only: unit0,unit, FFOOT, FINCH, DONE, DZERO
	IMPLICIT NONE
	REAL*8, DIMENSION(0 : NBOU, 0 : MAXS) :: EPS, SIG, VECOD
	REAL*8, DIMENSION(NBOU) :: HOR
	INTEGER :: NBOU, MAXS, MAXI
      !
	INTEGER :: IMAX
	REAL*8  :: SAPOINT, DELTA
C
      !INTEGER :: unit = 2
	INTEGER :: I, J, NumLyrs, NumCldLyrs; REAL*8 :: ZSTART, ZEND 
!============================================================================
C
	!OPEN(UNIT, FILE = 'IND.INP', STATUS='OLD')
	!   READ(UNIT, *) NumLyrs, NumCldLyrs  !!!
	!   NBOU = NumLyrs - 1
	!   MAXS = NumCldLyrs - 1
	!   MAXI = MAXS
C
	!   ALLOCATE(EPS(0 : NBOU, 0 : MAXS), SIG(0 : NBOU, 0 : MAXS))
	!   ALLOCATE(VECOD(0 : NBOU, 0 : MAXS))
	   EPS = DZERO; SIG = DZERO; VECOD = DZERO
C
	   IF(NBOU > 0) THEN
	      !ALLOCATE(HOR(NBOU));  
            HOR = DZERO
	      READ(UNIT, *) (HOR(I), I = 1, NBOU)  !!HOR = HOR * FFOOT
         END IF
C
	   DO I = 0, NBOU; READ(UNIT, *) (VECOD(I, J), J = 1, MAXS); END DO
	   !VECOD = VECOD * FINCH
C
	   DO I = 0, NBOU; READ(UNIT, *) (SIG(I, J), J = 0, MAXS); END DO
	   !SIG = DONE / SIG;  EPS = DONE
         !EPS = 10.0
C
	   READ(UNIT, *) ZSTART, ZEND, DELTA
	   ZSTART = ZSTART * FFOOT; ZEND = ZEND * FFOOT; DELTA = DELTA * FINCH
	   IMAX = INT(ABS((ZSTART - ZEND) / DELTA)) + 1
	   SAPOINT = ZSTART
         !
         SAPOINT = ZSTART/FFOOT; DELTA = DELTA/FINCH !SAPOINT in ft, DELTA in inch
      CLOSE(UNIT)
      !
      call Output_INDAT0input(NBOU,MAXS, SIG,VECOD,HOR, SAPOINT,
     1                        DELTA,IMAX)
C
      END SUBROUTINE INDAT0
! 
!***********************************************************************************
!
      SUBROUTINE Output_INDAT0input(NBOU,MAXS, SIG,VECOD,HOR, SAPOINT,
     1                              DELTA,IMAX)
	USE GlobalVarbls, only: unit0, FFOOT, FINCH, DONE, DZERO, INCH2FOOT
	IMPLICIT NONE
	REAL*8, DIMENSION(0 : NBOU, 0 : MAXS) :: EPS, SIG, VECOD
	REAL*8, DIMENSION(NBOU) :: HOR
	INTEGER :: NBOU, MAXS, MAXI
      !
	INTEGER :: IMAX
	REAL*8  :: SAPOINT, DELTA
C
      !INTEGER :: unit = 2
	INTEGER :: I, J, NumLyrs, NumCldLyrs; REAL*8 :: ZSTART, ZEND 
!============================================================================
C
      write(unit0,*)'INPUT FILE : IND.INP'
      write(unit0,*)NBOU+1,MAXS+1
	   IF(NBOU > 0) THEN
	      write(unit0,*) (HOR(I), I = 1, NBOU)  !!HOR = HOR * FFOOT
         END IF
C
	   DO I = 0, NBOU;  
             write(unit0,*)  (VECOD(I, J), J = 1, MAXS); END DO
C
	   DO I = 0, NBOU;  write(unit0,*) (SIG(I, J), J = 0, MAXS); END DO
C
         ZEND=SAPOINT+(IMAX-1)*DELTA*INCH2FOOT
	   write(unit0,*) SAPOINT, ZEND, DELTA, ' !ZSTART, ZEND, DELTA '
         write(unit0,*) IMAX, '! TOTAL LOG POINT NUMBER'
         !
      END SUBROUTINE Output_INDAT0input
! 
!***********************************************************************************
!
      SUBROUTINE INDAT(NBOU0, MAXS0, SIG0,VECOD0,HOR0, 
     1                                     SAPOINT0,DELTA0,IMAX0)
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
      INTEGER:: NBOU0, MAXS0
	REAL*8, DIMENSION(0 : NBOU0, 0 : MAXS0) :: SIG0, VECOD0
	REAL*8, DIMENSION(NBOU0) :: HOR0
      REAL(8):: SAPOINT0,DELTA0; INTEGER:: IMAX0
      !
C
	INTEGER :: I, J, NumLyrs, NumCldLyrs; REAL*8 :: ZSTART, ZEND 
C
	!OPEN(UNIT, FILE = 'IND.INP', STATUS='OLD')
	!   READ(UNIT, *) NumLyrs, NumCldLyrs  !!!
	!   NBOU = NumLyrs - 1
	!   MAXS = NumCldLyrs - 1
	!   MAXI = MAXS
      NBOU = NBOU0; MAXS = MAXS0;  MAXI = MAXS
C
      IF (ALLOCATED(EPS)) DEALLOCATE(EPS); 
      IF (ALLOCATED(SIG)) DEALLOCATE(SIG); 
      IF (ALLOCATED(VECOD)) DEALLOCATE(VECOD); 
	   ALLOCATE(EPS(0 : NBOU, 0 : MAXS), SIG(0 : NBOU, 0 : MAXS))
	   ALLOCATE(VECOD(0 : NBOU, 0 : MAXS));
         !
	   !EPS = DZERO; SIG = DZERO; VECOD = DZERO
         SIG = SIG0; VECOD = VECOD0
C
	   IF(NBOU > 0) THEN
            IF (ALLOCATED(HOR)) DEALLOCATE(HOR); 
	      ALLOCATE(HOR(NBOU));  !HOR = DZERO
	      !READ(UNIT, *) (HOR(I), I = 1, NBOU)  !
            HOR = HOR0
	      HOR = HOR * FFOOT
         END IF
C
	   !DO I = 0, NBOU; READ(UNIT, *) (VECOD(I, J), J = 1, MAXS); END DO
         VECOD=  VECOD0  
	   VECOD = VECOD * FINCH
C
	   !DO I = 0, NBOU; READ(UNIT, *) (SIG(I, J), J = 0, MAXS); END DO
	   SIG = DONE / SIG;  EPS = DONE;  !EPS = 10.0
C
	   !READ(UNIT, *) ZSTART, ZEND, DELTA
	   !ZSTART = ZSTART * FFOOT; ZEND = ZEND * FFOOT; DELTA = DELTA * FINCH
	   !IMAX = INT(ABS((ZSTART - ZEND) / DELTA)) + 1
	   !SAPOINT = ZSTART
         !DELTA0 = DELTA0 * FINCH;
         !SAPOINT0 = SAPOINT0 * FFOOT
         !DELTA = DELTA0; IMAX = IMAX0; SAPOINT = SAPOINT0;
         !
         DELTA = DELTA0 * FINCH;  IMAX = IMAX0;
         SAPOINT = SAPOINT0 * FFOOT
      CLOSE(UNIT)
C
       END SUBROUTINE INDAT    
 !        
 !
      SUBROUTINE INDAT_wang
C-----------------------------------------------------------------------
	USE GlobalVarbls
	IMPLICIT NONE
C
	INTEGER :: I, J, NumLyrs, NumCldLyrs; REAL*8 :: ZSTART, ZEND 
C
	OPEN(UNIT, FILE = 'IND.INP', STATUS='OLD')
	   READ(UNIT, *) NumLyrs, NumCldLyrs  !!!
	   NBOU = NumLyrs - 1
	   MAXS = NumCldLyrs - 1
	   MAXI = MAXS
C
	   ALLOCATE(EPS(0 : NBOU, 0 : MAXS), SIG(0 : NBOU, 0 : MAXS))
	   ALLOCATE(VECOD(0 : NBOU, 0 : MAXS))
	   EPS = DZERO; SIG = DZERO; VECOD = DZERO
C
	   IF(NBOU > 0) THEN
	      ALLOCATE(HOR(NBOU));  HOR = DZERO
	      READ(UNIT, *) (HOR(I), I = 1, NBOU)  !
	      HOR = HOR * FFOOT
         END IF
C
	   DO I = 0, NBOU; READ(UNIT, *) (VECOD(I, J), J = 1, MAXS); END DO
	   VECOD = VECOD * FINCH
C
	   DO I = 0, NBOU; READ(UNIT, *) (SIG(I, J), J = 0, MAXS); END DO
	   SIG = DONE / SIG;  EPS = DONE
         !EPS = 10.0
C
	   READ(UNIT, *) ZSTART, ZEND, DELTA
	   ZSTART = ZSTART * FFOOT; ZEND = ZEND * FFOOT; DELTA = DELTA * FINCH
	   IMAX = INT(ABS((ZSTART - ZEND) / DELTA)) + 1
	   SAPOINT = ZSTART
      CLOSE(UNIT)
C
         END           
!
!********************************************************************************************
!    The End
!********************************************************************************************
!        