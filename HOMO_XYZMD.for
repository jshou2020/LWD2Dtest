c
!******************************************************************************
C
c    File      : HOMO_XYZMD.FOR
!    Purpose   : computation of magnetic components in the formation system
C    SUBROUTINE: XMDINF(), YMDINF(), ZMDINF()
!
c    Date      : 03/18/2009
C    Update    : 09/18/2009
!    Author    : Hou, Junsheng
C    Modified from UH V1D codes
C
!******************************************************************************
!
C  SUBROUTINE XMDINF(XKH,XKV,XL,RO,PR,Z,X0,HXX,HXZ)
C  HORIZONTAL X-ORIENTED MAGNETIC DIPOLE AT THE ORIGIN IN AN
C  INFINITE MEDIUM ; THE HR & HZ ARE CALCULATED AT RO,PHI,Z .
C  PR = PHI IN RADIANS : CONVB. INTO RAD IS DONE IN THE MAIN PROG.
C  MAGNETIC MOMENT IS UNITY.
!
!******************************************************************************
c
      SUBROUTINE XMDINF(XKH,XKV,XL,RO,PR,Z,X0,HXX,HXZ)
	USE CONSTANT
	IMPLICIT NONE
	COMPLEX*16 XKH,XKV,HXX,HXZ
	REAL*8 RO,PR,Z,XL,X0
C
	COMPLEX*16 S
	REAL*8 X,Y,R
C
C	X=RO*DCOS(PR)
	X=X0
	Y=RO*DSIN(PR)
      R=DSQRT(RO*RO+Z*Z)
	S=DSQRT(RO*RO+XL*XL*Z*Z)
! Hxx
      HXX=CDEXP(XJ*XKV*S)*(XKV*XKH/S+XJ*XKH/(RO*RO)-XKV*XKH*DCOS(PR)*
     +   DCOS(PR)/S-2*XJ*XKH*DCOS(PR)*DCOS(PR)/(RO*RO))+CDEXP(XJ*XKH*R)*
     +   (-XJ*XKH/(RO*RO)+XKH*XKH*DCOS(PR)*DCOS(PR)/R+2*XJ*XKH*DCOS(PR)*
     +   DCOS(PR)/(RO*RO)+XJ*XKH/(R*R)-(XKH*XKH*X*X+1)/(R*R*R)-3*XJ*XKH*
     +    X*X/(R*R*R*R)+3*X*X/(R*R*R*R*R))
	HXX=HXX/(4.D0*PI)
c      XHY=CDEXP(XJ*XKV*S)*X*Y*(-XKH*XKV/S-2*XJ*XKH/(RO*RO))+
c     +    CDEXP(XJ*XKH*R)*X*Y*(2*XJ*XKH/(RO*RO*RO*RO)+XKH*XKH/(R*RO*RO)
c     +    -XKH*XKH/(R*R*R)-3*XJ*XKH/(R*R*R*R)+3/(R*R*R*R*R))
c	XHY=XHY/(4.D0*PI)
c      IF(RO<1.0D-4)XHY=0.0
! Hxz
      HXZ=CDEXP(XJ*XKH*R)*X*Z*(-XKH*XKH-3*XJ*XKH/R+3/(R*R))/(R*R*R)
c      IF(RO<1.0D-4)HXZ=0.0
	HXZ=HXZ/(4.D0*PI)
!
      RETURN
      END SUBROUTINE XMDINF
C
!******************************************************************************
!
C  SUBROUTINE YMDINF()
C  HORIZONTAL Y-ORIENTED MAGNETIC DIPOLE AT THE ORIGIN IN AN
C  INFINITE MEDIUM ; THE HR & HZ ARE CALCULATED AT RO,PHI,Z .
C  PR = PHI IN RADIANS : CONVB. INTO RAD IS DONE IN THE MAIN PROG.
C  MAGNETIC MOMENT IS UNITY.
!
!******************************************************************************
C
      SUBROUTINE YMDINF(XKH,XKV,YL,RO,PR,Z,HYY)
	USE CONSTANT
	IMPLICIT NONE
	COMPLEX*16 XKH,XKV,HYY
      REAL*8 S,R,X,Y,RO,PR,Z,YL
C
	X=RO*COS(PR)
	Y=RO*SIN(PR)
      R=DSQRT(RO*RO+Z*Z)
	S=DSQRT(RO*RO+YL*YL*Z*Z)
! Hyy:
      HYY=CDEXP(XJ*XKV*S)*(XKV*XKH/S+XJ*XKH/(RO*RO)-XKV*XKH*DSIN(PR)*
     +   SIN(PR)/S-2*XJ*XKH*DSIN(PR)*DSIN(PR)/(RO*RO))+CDEXP(XJ*XKH*R)*
     +   (-XJ*XKH/(RO*RO)+XKH*XKH*DSIN(PR)*DSIN(PR)/R+2*XJ*XKH*DSIN(PR)*
     +    SIN(PR)/(RO*RO)+XJ*XKH/(R*R)-(XKH*XKH*Y*Y+1)/(R*R*R)-3*XJ*XKH*
     +    Y*Y/(R*R*R*R)+3*Y*Y/(R*R*R*R*R))
      HYY=HYY/(4.0*PI)
!      YHZ=CDEXP(XJ*XKH*R)*Y*Z*(-XKH*XKH-3*XJ*XKH/R+3/(R*R))/(R*R*R)
!      YHZ=YHZ/(4.0*PI)
!
      RETURN
      END SUBROUTINE YMDINF
C
!******************************************************************************
!
C  SUBROUTINE ZMDINF()
C  Z-directional MAGNETIC DIPOLE AT THE ORIGIN IN AN INFINITE MEDIUM.
C  THE HR & HZ FIELDS ARE EVALUATED AT (RO,Z).
C  MAGNETIC MOMENT IS UNITY.
!
!******************************************************************************
C
      SUBROUTINE ZMDINF(XKH,XKV,RO,Z,HZZ)
	USE CONSTANT
	IMPLICIT NONE
      COMPLEX*16 XKH,XKV,HZZ
	REAL*8 RO,Z
	REAL*8 R
C
      R=DSQRT(RO*RO+Z*Z)
! Hzz
      HZZ=CDEXP(XJ*XKH*R)*(XKH*XKH+XJ*XKH/R-(XKH*XKH*Z*Z+1)/(R*R)
     +                    -3*XJ*XKH*Z*Z/(R*R*R)+3*Z*Z/(R*R*R*R))/R
      HZZ=HZZ/(4.0*PI)
!
      RETURN
      END SUBROUTINE ZMDINF
c
!******************************************************************************
c