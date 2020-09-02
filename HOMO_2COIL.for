c
!**********************************************************************
!
c    File      : HOMO_2COIL.FOR
c    SUBROUTINE: HOMO_2COIL(DIP,FREQ,Z0,RH,RV,CXX,CXZ,CYY,CZZ)
!    Purpose   : control of Hxx, Hxz, Hyy, and Hzz calculation   
!                and final apparent conductivity in borehole system
c    Date      : 03/18/2009
C    Update    : 11/23/2009
C    Author    : Hou, Junsheng 
!
!**********************************************************************
c
      SUBROUTINE HOMO_2COIL(DIP,FREQ,Z0,RH,RV,CXX,CXZ,CYY,CZZ)
	USE CONSTANT
	IMPLICIT NONE
!inputs and Outputs::
      REAL(8),INTENT(IN) :: DIP,FREQ,Z0,RH,RV
	REAL(8),INTENT(OUT):: CXX(2),CXZ(2),CYY(2),CZZ(2)
!locals:
	REAL(8) RO,PR,OMEGA,XL,Z,X
	REAL(8) CH,CV
	REAL(8) KXX,KXZ,KZZ
	COMPLEX*16 KH,KV,HXX,HXZ,HYY,HZZ
	REAL(8) CF(3,3),CB(3,3)
      INTEGER I,ICH
      COMMON/COND_HF/ICH
C
c FREQ: unit is kHz
c
	OMEGA=2.0*PI*FREQ*1000.D0
c	PR=1.0D-4 !(0.00001)
c      RO=1.0D-4
	PR=1.0D-5 !(0.00001) better case
      RO=1.0D-5
c	PR=1.0D-6 !(0.00001) no
c      RO=1.0D-6
c	PR=1.0D-7 !(0.00001) error
c      RO=1.0D-7
c	PR=1.0D-8 !(0.00001) error
c      RO=1.0D-8
	CH=1.0/RH
	CV=1.0/RV
	KH=CDSQRT(XJ*OMEGA*UO*CH)
	KV=CDSQRT(XJ*OMEGA*UO*CV)
	XL=CH/CV
	XL=DSQRT(XL)
!
! (1).(x,y,z) coordinate rotation from borehole to formation system: 
!     (xb,yb,zb) -> (xf,yf,zf)
!
!     if we know the tool location in the tool/borehole system (xb,yb,zb), the tool 
!     location should be expressed as the values in the formation system(xf,yf,zf). 
!
      Z=Z0*0.0254
	Z=Z*COSD(DIP) !!NEW
	RO=ABS(-RO-Z0*0.0254*SIND(DIP))  !!NEW
      X=-Z0*0.0254*SIND(DIP)
!
! (2).COMPUTATION of a complex tensor for nine magnetic components 
!     in formation system (xf,yf,zf)
!
    	CALL XMDINF(KH,KV,XL,RO,PR,Z,X,HXX,HXZ) !X-direction dipole
	CALL YMDINF(KH,KV,XL,RO,PR,Z,HYY)       !Y-direction dipole
	CALL ZMDINF(KH,KV,RO,Z,HZZ)             !Z-direction dipole
c
C	write(*,*)'IMAG(HXX)*1000=',IMAG(HXX)*1000
C	write(*,*)'IMAG(HXZ)*1000=',IMAG(HXZ)*1000
C	write(*,*)'IMAG(HYY)*1000=',IMAG(HYY)*1000
C	write(*,*)'IMAG(HZZ)*1000=',IMAG(HZZ)*1000
C	write(*,*)
!
! (3).Rotation (x-z plane) of 2-COIL real fields (or uncalibrated R, X-signal) 
!     from formation system to TOOL/measurement system
!     Note:complex Tensor rotation issue.
!          R.real(C).transpose (R) is equal to real(R.C.transpose (R)).
c
C	KZZ=OMEGA*UO/(4.0*PI*Z)
C	KXX=KZZ/2.0
C	KXZ=KZZ/4.0
C      CXX=IMAG(HXX)
C      CXZ=IMAG(HXZ)
C      CYY=IMAG(HYY)
C      CZZ=IMAG(HZZ)
C	CXX=CXX/KXX*1000.   !XX(mS/m)
C	CXZ=CXZ/KXZ*1000.   !XZ
C	CYY=CYY/KXX*1000.   !YY
C	CZZ=CZZ/KZZ*1000.   !ZZ(mS/m)
C R-signal:
      CF=0.0D0
      CF(1,1)=IMAG(HXX)
      CF(1,3)=IMAG(HXZ)
      CF(2,2)=IMAG(HYY)
      CF(3,1)=IMAG(HXZ)
      CF(3,3)=IMAG(HZZ)
      CALL CFB_xzROTATION(DIP,CF,CB)
      CXX(1)=CB(1,1)
      CXZ(1)=CB(1,3)
      CYY(1)=CB(2,2)
      CZZ(1)=CB(3,3)
C X-signal:
      CF=0.0D0
      CF(1,1)=REAL(HXX)
      CF(1,3)=REAL(HXZ)
      CF(2,2)=REAL(HYY)
      CF(3,1)=REAL(HXZ)
      CF(3,3)=REAL(HZZ)
      CALL CFB_xzROTATION(DIP,CF,CB)
      CXX(2)=CB(1,1)
      CXZ(2)=CB(1,3)
      CYY(2)=CB(2,2)
      CZZ(2)=CB(3,3)
c
C	RETURN
!
C     CXX2=HXX*COSD(DIP)**2+HZZ*SIND(DIP)**2-2*HXZ*SIND(DIP)*COSD(DIP)
C     CZZ2=HXX*SIND(DIP)**2+HZZ*COSD(DIP)**2+2*HXZ*SIND(DIP)*COSD(DIP)
!
C	 Write(*,*) "IMAG(Hxx)*1000=", (CXX)*1000
C	 Write(*,*) "IMAG(Hxy)*1000=", (CXy)*1000
C	 Write(*,*) "IMAG(Hxz)*1000=", (CXZ)*1000
C	 Write(*,*) "IMAG(Hyx)*1000=", (Cyx)*1000
C	 Write(*,*) "IMAG(Hyy)*1000=", (Cyy)*1000
C	 Write(*,*) "IMAG(Hyz)*1000=", (Cyz)*1000
C	 Write(*,*) "IMAG(Hzx)*1000=", (Czx)*1000
C	 Write(*,*) "IMAG(Hzy)*1000=", (Czy)*1000
C	 Write(*,*) "IMAG(Hzz)*1000=", (Czz)*1000
c
!  (4).Calibration for 2-coil apparent conductivity in TOOL/measurement system
!
!=============================================================================
      IF(ICH==1)THEN
	 KZZ=OMEGA*UO/(4.0*PI*Z0*0.0254)  !ZZ
	 KXX=KZZ/2.0                      !XX,YY,xy,yx
	 KXZ=-KZZ/4.0                      !XZ,ZX,zy,yz(==new==)
c	 KXZ= KZZ/4.0                     !XZ,ZX,zy,yz(original)
!=============================================================================
	ELSE
	 KZZ=1.0
	 KXX=1.0
	 KXZ=1.0
	ENDIF
c
      DO I=1,2
	 CXX(I)=CXX(I)/KXX*1000.   !XX(conductivity:mS/m or H field: mA/m)
	 CXZ(I)=CXZ(I)/KXZ*1000.   !XZ
	 CYY(I)=CYY(I)/KXX*1000.   !YY
	 CZZ(I)=CZZ(I)/KZZ*1000.   !ZZ(mS/m)
      ENDDO
!
	RETURN 
	END SUBROUTINE HOMO_2COIL
c
!**********************************************************************
c
!
!**********************************************************************
c
      SUBROUTINE HOMO_2COIL02(DIP,FREQ,Z0,RH,RV,HZZ)
	USE CONSTANT
	IMPLICIT NONE
!inputs and Outputs::
      REAL(8),INTENT(IN) :: DIP,FREQ,Z0,RH,RV
	!REAL(8),INTENT(OUT):: CXX(2),CXZ(2),CYY(2),CZZ(2)
	REAL(8):: CXX(2),CXZ(2),CYY(2),CZZ(2)
!locals:
	REAL(8) RO,PR,OMEGA,XL,Z,X
	REAL(8) CH,CV
	REAL(8) KXX,KXZ,KZZ
	COMPLEX*16 KH,KV,HXX,HXZ,HYY,HZZ
	REAL(8) CF(3,3),CB(3,3)
      INTEGER I,ICH
      COMMON/COND_HF/ICH
C
c FREQ: unit is kHz
c
	OMEGA=2.0*PI*FREQ*1000.D0
c	PR=1.0D-4 !(0.00001)
c      RO=1.0D-4
	PR=1.0D-5 !(0.00001) better case
      RO=1.0D-5
c	PR=1.0D-6 !(0.00001) no
c      RO=1.0D-6
c	PR=1.0D-7 !(0.00001) error
c      RO=1.0D-7
c	PR=1.0D-8 !(0.00001) error
c      RO=1.0D-8
	CH=1.0/RH
	CV=1.0/RV
	KH=CDSQRT(XJ*OMEGA*UO*CH)
	KV=CDSQRT(XJ*OMEGA*UO*CV)
	XL=CH/CV
	XL=DSQRT(XL)
!
! (1).(x,y,z) coordinate rotation from borehole to formation system: 
!     (xb,yb,zb) -> (xf,yf,zf)
!
!     if we know the tool location in the tool/borehole system (xb,yb,zb), the tool 
!     location should be expressed as the values in the formation system(xf,yf,zf). 
!
      Z=Z0*0.0254
	Z=Z*COSD(DIP) !!NEW
	RO=ABS(-RO-Z0*0.0254*SIND(DIP))  !!NEW
      X=-Z0*0.0254*SIND(DIP)
!
! (2).COMPUTATION of a complex tensor for nine magnetic components 
!     in formation system (xf,yf,zf)
!
    	CALL XMDINF(KH,KV,XL,RO,PR,Z,X,HXX,HXZ) !X-direction dipole
	CALL YMDINF(KH,KV,XL,RO,PR,Z,HYY)       !Y-direction dipole
	CALL ZMDINF(KH,KV,RO,Z,HZZ)             !Z-direction dipole
c
C	write(*,*)'IMAG(HXX)*1000=',IMAG(HXX)*1000
C	write(*,*)'IMAG(HXZ)*1000=',IMAG(HXZ)*1000
C	write(*,*)'IMAG(HYY)*1000=',IMAG(HYY)*1000
C	write(*,*)'IMAG(HZZ)*1000=',IMAG(HZZ)*1000
C	write(*,*)
!
! (3).Rotation (x-z plane) of 2-COIL real fields (or uncalibrated R, X-signal) 
!     from formation system to TOOL/measurement system
!     Note:complex Tensor rotation issue.
!          R.real(C).transpose (R) is equal to real(R.C.transpose (R)).
c
C	KZZ=OMEGA*UO/(4.0*PI*Z)
C	KXX=KZZ/2.0
C	KXZ=KZZ/4.0
C      CXX=IMAG(HXX)
C      CXZ=IMAG(HXZ)
C      CYY=IMAG(HYY)
C      CZZ=IMAG(HZZ)
C	CXX=CXX/KXX*1000.   !XX(mS/m)
C	CXZ=CXZ/KXZ*1000.   !XZ
C	CYY=CYY/KXX*1000.   !YY
C	CZZ=CZZ/KZZ*1000.   !ZZ(mS/m)
C R-signal:
      CF=0.0D0
      CF(1,1)=IMAG(HXX)
      CF(1,3)=IMAG(HXZ)
      CF(2,2)=IMAG(HYY)
      CF(3,1)=IMAG(HXZ)
      CF(3,3)=IMAG(HZZ)
      CALL CFB_xzROTATION(DIP,CF,CB)
      CXX(1)=CB(1,1)
      CXZ(1)=CB(1,3)
      CYY(1)=CB(2,2)
      CZZ(1)=CB(3,3)
C X-signal:
      CF=0.0D0
      CF(1,1)=REAL(HXX)
      CF(1,3)=REAL(HXZ)
      CF(2,2)=REAL(HYY)
      CF(3,1)=REAL(HXZ)
      CF(3,3)=REAL(HZZ)
      CALL CFB_xzROTATION(DIP,CF,CB)
      CXX(2)=CB(1,1)
      CXZ(2)=CB(1,3)
      CYY(2)=CB(2,2)
      CZZ(2)=CB(3,3)
c
C	RETURN
!
C     CXX2=HXX*COSD(DIP)**2+HZZ*SIND(DIP)**2-2*HXZ*SIND(DIP)*COSD(DIP)
C     CZZ2=HXX*SIND(DIP)**2+HZZ*COSD(DIP)**2+2*HXZ*SIND(DIP)*COSD(DIP)
!
C	 Write(*,*) "IMAG(Hxx)*1000=", (CXX)*1000
C	 Write(*,*) "IMAG(Hxy)*1000=", (CXy)*1000
C	 Write(*,*) "IMAG(Hxz)*1000=", (CXZ)*1000
C	 Write(*,*) "IMAG(Hyx)*1000=", (Cyx)*1000
C	 Write(*,*) "IMAG(Hyy)*1000=", (Cyy)*1000
C	 Write(*,*) "IMAG(Hyz)*1000=", (Cyz)*1000
C	 Write(*,*) "IMAG(Hzx)*1000=", (Czx)*1000
C	 Write(*,*) "IMAG(Hzy)*1000=", (Czy)*1000
C	 Write(*,*) "IMAG(Hzz)*1000=", (Czz)*1000
c
!  (4).Calibration for 2-coil apparent conductivity in TOOL/measurement system
!
!=============================================================================
      IF(ICH==1)THEN
	 KZZ=OMEGA*UO/(4.0*PI*Z0*0.0254)  !ZZ
	 KXX=KZZ/2.0                      !XX,YY,xy,yx
	 KXZ=-KZZ/4.0                      !XZ,ZX,zy,yz(==new==)
c	 KXZ= KZZ/4.0                     !XZ,ZX,zy,yz(original)
!=============================================================================
	ELSE
	 KZZ=1.0
	 KXX=1.0
	 KXZ=1.0
	ENDIF
c
      DO I=1,2
	 CXX(I)=CXX(I)/KXX*1000.   !XX(conductivity:mS/m or H field: mA/m)
	 CXZ(I)=CXZ(I)/KXZ*1000.   !XZ
	 CYY(I)=CYY(I)/KXX*1000.   !YY
	 CZZ(I)=CZZ(I)/KZZ*1000.   !ZZ(mS/m)
      ENDDO
!
	RETURN 
	END SUBROUTINE HOMO_2COIL02
c
!**********************************************************************
c