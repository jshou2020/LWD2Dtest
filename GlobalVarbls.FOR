	MODULE GlobalVarbls
C-----------------------------------------------------------------------
C Module to define the parameters, constants and global variables.
C-----------------------------------------------------------------------
	REAL*8, PARAMETER :: PI = 3.1415926536D0
	REAL*8, PARAMETER :: U0 = 1.256637D-6, E0 = 8.85D-12
	REAL*8, PARAMETER :: CURENT = 100.0D0;  !org: CURENT = 10.0D0
	REAL*8, PARAMETER :: DZERO = 0.0D0
!      
	REAL*8, PARAMETER :: DONE = 1.0D0, DTWO = 2.0D0
	COMPLEX*16, PARAMETER :: CZERO = (0.0D0, 0.0D0)
	REAL*8, PARAMETER :: FINCH = 0.0254D0, FFOOT = 12.0D0 * FINCH
      REAL*8, PARAMETER :: INCH2FOOT=1.0/12.0
	INTEGER, PARAMETER :: UNIT = 2, UNIT0 = 200;
C-----------------------------------------------------------------------
!      
      INTEGER :: Itool         ! Itool= -4: EWR-M5; =-5: ARC; =-6: MPR
      CHARACTER*10:: ToolName  ! ToolName   EWR-M5;      ARC;      MPR
C-----------------------------------------------------------------------
C N:		Number of basis functions, which is equal to the number of  
C		grid points, exclusive of the leftmost point (r=0) and  
C		rightmost point (r=r_max).
C MAXORDER / 2 : 
C		Order of basis functions
C Q:		Parameter controling the distribution of grid points. 
C RMAX:	Outmost boundary
!org:
C-----------------------------------------------------------------------
!	INTEGER, PARAMETER :: N = 60, MAXORDER = 6
!	REAL*8, PARAMETER :: Q = 1.215D0, RMAX = 100.0D0
C-----------------------------------------------------------------------
Cnew01:
!-----------------------------------------------------------------------
!	INTEGER, PARAMETER :: N = 100, MAXORDER = 6
!	REAL*8, PARAMETER ::  RMAX = 100.0D0         !meter
!	REAL*8, PARAMETER ::  Q = 1.215D0
C-----------------------------------------------------------------------
Cnew02:
!-----------------------------------------------------------------------
	INTEGER, PARAMETER::  N = 150, MAXORDER = 6
	REAL*8, PARAMETER ::  RMAX = 20.0D0
	REAL*8, PARAMETER ::  Q = 1.1215D0
C-----------------------------------------------------------------------
!
C OMEGA:	Angular frequency
C S0T:	Area of transmitter
C S0R:	Area of receiver
C RADIUST:	Radius of transmitter
C RADIUST:	Radius of receiver
C-----------------------------------------------------------------------
	REAL*8, SAVE :: OMEGA, S0T, S0R, RADIUST, RADIUSR
C-----------------------------------------------------------------------
C EPS:	Relative permittivity of beds
C SIG:	Conductivity of beds
C VECOD:	Borehole radius and invasion front radius of beds
C HOR:	Location of bed boundaries
C NBOU:	Number of boundaries
C MAXS:	Number of radial layers in each bed
C MAXI:	Maximum number of radial layer boundaries in one element
C-----------------------------------------------------------------------
	REAL*8, DIMENSION(:, :), ALLOCATABLE, SAVE :: EPS, SIG, VECOD
	REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: HOR
	INTEGER, SAVE :: NBOU, MAXS, MAXI
C-----------------------------------------------------------------------
C R:		Location of grid points 
C-----------------------------------------------------------------------
	REAL*8, DIMENSION(0 : N + 1), SAVE :: R
C-----------------------------------------------------------------------
C RM:		Location of radial boundaries in an element
C K2:		Wavenumber square of sub-areas in an element
C IN:		Number of radial boundaries in an element
C-----------------------------------------------------------------------
	REAL*8, DIMENSION(:, :, :), ALLOCATABLE, SAVE :: RM
	COMPLEX*16, DIMENSION(:, :, :), ALLOCATABLE, SAVE :: K2
	INTEGER, DIMENSION(:, :), ALLOCATABLE, SAVE :: IN	
C-----------------------------------------------------------------------
C A:		Maxrix A in the eigenvalue problems
C B & B0:	Maxrix B in the eigenvalue problems
C-----------------------------------------------------------------------
	REAL*8, DIMENSION(N, N), SAVE :: A, B, B0
C-----------------------------------------------------------------------
C C:		Eigenvectors of each bed
C V:		Eigenvalues of each bed
C-----------------------------------------------------------------------
	COMPLEX*16, DIMENSION(:, :, :), ALLOCATABLE, SAVE :: C
	COMPLEX*16, DIMENSION(:, :), ALLOCATABLE, SAVE :: V
C-----------------------------------------------------------------------
C QU:		Global reflection matrices for up-going waves
C QD:		Global reflection matrices for down-going waves
C QU:		Global transmission matrices for up-going waves
C QD:		Global transmission matrices for down-going waves
C P:		Coupling matrices
C-----------------------------------------------------------------------
	COMPLEX*16, DIMENSION(:, :, :), ALLOCATABLE, SAVE :: QU, QD
	COMPLEX*16, DIMENSION(:, :, :), ALLOCATABLE, SAVE :: SU, SD, P
	COMPLEX*16, DIMENSION(:, :), ALLOCATABLE, SAVE :: E
C-----------------------------------------------------------------------
C UP:		Up-going waves at the source plane
C UD:		Down-going waves at the source plane
C U:		Up-going waves at the bed boundaries above the source, and
C		Down-going waves at the bed boundaries below the source
C G:		Basis functions at a given point
C-----------------------------------------------------------------------
	COMPLEX*16, DIMENSION(N), SAVE :: UP, UD
	REAL*8, DIMENSION(N), SAVE :: G
	COMPLEX*16, DIMENSION(:, :), ALLOCATABLE, SAVE :: U
C-----------------------------------------------------------------------
C MAXT:	Maximum number of transmiters
C MAXR:	Maximum number of receivers
C NOCL:	Number of transmiters and receivers of an actual tool
C NTT:	Number of transmiters of an actual tool
C NRR:	Number of receivers of an actual tool
C PTT:	Position of transmiters in the tool coordinates 
C PRR:	Position of receivers in the tool coordinates 
C MTT:	Number of turns of transmiters
C MRR:	Number of turns of receivers
C-----------------------------------------------------------------------
	INTEGER, PARAMETER :: MAXT = 20, MAXR = 20
	INTEGER, SAVE :: NOCL, NTT, NRR
	REAL*8, DIMENSION(MAXT), SAVE :: PTT, MMT
	REAL*8, DIMENSION(MAXR), SAVE :: PRR, MMR
C
	REAL*8, SAVE :: TTK
	REAL*8, DIMENSION(MAXT, MAXR), SAVE :: TLK
C
	INTEGER, SAVE :: IMAX
	REAL*8, SAVE :: SAPOINT, DELTA
	COMPLEX*16, DIMENSION(:), ALLOCATABLE, SAVE :: VR
C
	END MODULE GlobalVarbls