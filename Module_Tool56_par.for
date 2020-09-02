c
!***********************************************************************************
c    File name   :  Module_tool56_par.for
c    Module name :  Tool5_par, Tool52_par, Tool6_par,  Tool7_par,
!    Purpose     :  Define TOOL 5, 52, and 6, 7 pars - OD, frequency,spacing; turns,...
c    Date        :  Wed, 10/09/2013
C    Update      :  Wed, 09/02/2020
!***********************************************************************************
!
!    MODULE Tool5_par    !halLWD tool Itool = -4 AND EWR-M5 
!***********************************************************************************
c
      MODULE Tool5_par; IMPLICIT NONE; SAVE
!  frequency: 250K, 500K, 2000K(=2M); 
      INTEGER, PARAMETER ::  NFREQ35=3; 
!  sub-array: 16,24, 32,40 48
      INTEGER, PARAMETER ::  NARRAYzz5=5  !NARRAYzz5=5
      !INTEGER, PARAMETER ::  NReceiver=6  !TOTAL receiver num for NARRAYzz5 subs
      INTEGER, PARAMETER ::  NReceiver5=6  !TOTAL receiver num for NARRAYzz5 subs
      !REAL(8), PARAMETER ::  ODtool5 = 5.0;  !Tool OD, unit: inch
      REAL(8), PARAMETER ::  ODtool5 = 1.0; 
      REAL(8),DIMENSION(NFREQ35):: FREQ35=(/250.,500.,2000./) !kHz
c 6-array spacing data:
      REAL(8):: Rspacing5(1:NReceiver5)=(/12,20,28,36,44,52/)  !inch)
!      
C Z-transmitter turn
      INTEGER:: NZTturn5=5    !T-coil turn
      REAL(8):: ZTClen5=0.5   !T-coil length, inch
      !INTEGER:: NZTturn5=3    !T-coil turn
      !REAL(8):: ZTClen5=0.25   !T-coil length, inch
      !
      !REAL(8):: ZTCdia5=6     !T-coil diameter, inch
      REAL(8):: ZTCdia5=2.0  
c 3-array ZZ-R-channel turn:
      !INTEGER:: NZRturn5(1:NReceiver)=(/5,5,5,5,5,5/) !Near -> far sub
      INTEGER:: NZRturn5(1:NReceiver5)=(/5,5,5,5,5,5/) !Near -> far sub
      !REAL(8):: ZRClen5(1:NReceiver)=(/0.5,0.5,0.5,0.5,0.5,0.5/) 
      REAL(8):: ZRClen5(1:NReceiver5)=(/0.5,0.5,0.5,0.5,0.5,0.5/) 
!      
      !REAL(8):: ZRClen5(1:2*NARRAYzz5)=(/0.25,0.25,0.25,0.25,0.25,0.25/)
                             !R-coil length,inch, Near -> far sub
      !REAL(8):: ZRCdia5(1:NReceiver)=(/6,6,6,6,6,6/)
      REAL(8):: ZRCdia5(1:NReceiver5)=(/6,6,6,6,6,6/)
                             !R-coil diameter,inch, Near -> far sub
      REAL(8):: rPipe5=ODtool5/2.0  !Pipe radius, inch
      END MODULE Tool5_par
      
!***********************************************************************************
!    MODULE Tool52_par; !halLWD tool Itool = -4 and EWR-P4 
!    NOTE: In EWR-phase 4 (EWR-P4), the transmitter at 39" operates at 1MHz 
!                                   while the other transmitters operate at 2MHz.
!***********************************************************************************
c
      MODULE Tool52_par; IMPLICIT NONE; SAVE
      INTEGER, PARAMETER ::  NFREQ352=2;   !2 frequency: 1000K, 2000K(=2M);
      !
      INTEGER, PARAMETER ::  NARRAYzz52=4  !4 T locations:  9,15,27, 39
      INTEGER, PARAMETER ::  NReceiver52=8  
      REAL(8), PARAMETER ::  ODtool52 = 1.0; 
      REAL(8),DIMENSION(NFREQ352):: FREQ352=(/1000.,2000./) !kHz
c 8 RECEIVERS:
      REAL(8):: Rspacing52(1:NReceiver52)=(/6,12,12,18,24,30,36,42/)  !inch)
!      
C Z-transmitter turn
      INTEGER:: NZTturn52=5    !T-coil turn
      REAL(8):: ZTClen52=0.5   !T-coil length, inch
      REAL(8):: ZTCdia52=2.0  !T-coil diameter, inch
c 3-array ZZ-R-channel turn:
      INTEGER:: NZRturn52(1:NReceiver52)=(/5,5,5,5,5,5,5,5/) !Near -> far sub
      REAL(8):: ZRClen52(1:NReceiver52)=(/0.5,.5,.5,.5,.5,.5,.5,0.5/) 
      REAL(8):: ZRCdia52(1:NReceiver52)=(/6,6,6,6,6,6,6,6/); !R-coil length,inch, Near -> far sub
                             !R-coil diameter,inch, Near -> far sub
      REAL(8):: rPipe52=ODtool52/2.0  !Pipe radius, inch
      END MODULE Tool52_par
c
!***********************************************************************************
c    MODULE Tool6_par; Itool = -5: 'ARC675,ARC825,ARC900', 'ARC312,ARC475','IMPulse'  
!***********************************************************************************
c
      MODULE Tool6_par; IMPLICIT NONE; SAVE
!  frequency: 400K, 2000K(=2M); 
      INTEGER, PARAMETER ::  NFREQ36=2; 
!  sub-array: /13,19,25,31,37,43/
      INTEGER, PARAMETER ::  NARRAYzz6=5   !sub number
      INTEGER, PARAMETER ::  NReceiver6=6  !TOTAL receiver num for NARRAYzz5 subs
      REAL(8), PARAMETER ::  ODtool6 = 1.0; !Tool OD, unit: inch
      REAL(8),DIMENSION(NFREQ36):: FREQ36=(/400.,2000./) !kHz
c spacing data  for 6 RECEIVERs 
      REAL(8):: Rspacing6(1:NReceiver6)=(/13,19,25,31,37,43/)
                         !inch)
C Z-transmitter turn
      INTEGER:: NZTturn6=5    !T-coil turn
      REAL(8):: ZTClen6=0.5   !T-coil length, inch
      REAL(8):: ZTCdia6=2.0  !T-coil diameter, inch
c 3-array ZZ-R-channel turn:
      INTEGER:: NZRturn6(1:NReceiver6)=(/5,5,5,5,5,5/) !Near -> far sub
      REAL(8):: ZRClen6(1:NReceiver6)=(/0.5,0.5,0.5,0.5,0.5,0.5/) 
      !REAL(8):: ZRClen5(1:2*NARRAYzz5)=(/0.25,0.25,0.25,0.25,0.25,0.25/)
                             !R-coil length,inch, Near -> far sub
      REAL(8):: ZRCdia6(1:NReceiver6)=(/6,6,6,6,6,6/)
                             !R-coil diameter,inch, Near -> far sub
      REAL(8):: rPipe6=ODtool6/2.0  !Pipe radius, inch
      !
	END MODULE Tool6_par
c
!***********************************************************************************
c    MODULE Tool7_par;  Itool = -6 - MPR, MPR95
!***********************************************************************************
c
      MODULE Tool7_par; IMPLICIT NONE; SAVE
!  frequency: 400K, 2000K(=2M); 
      INTEGER, PARAMETER ::  NFREQ37=2; 
!  sub-array: /22,35/
      INTEGER, PARAMETER ::  NARRAYzz7=2   !sub number
      INTEGER, PARAMETER ::  NReceiver7=4  !TOTAL receiver num for NARRAYzz5 subs
      REAL(8), PARAMETER ::  ODtool7 = 1.0; !Tool OD, unit: inch
      REAL(8),DIMENSION(NFREQ37):: FREQ37=(/400.,2000./) !kHz
c spacing data  for 4 RECEIVERs 
      REAL(8):: Rspacing7(1:NReceiver7)=(/18.375,26.375,31.625,39.625/)
                         !inch)
C Z-transmitter turn
      INTEGER:: NZTturn7=5    !T-coil turn
      REAL(8):: ZTClen7=0.5   !T-coil length, inch
      REAL(8):: ZTCdia7=2.0  !T-coil diameter, inch
c 3-array ZZ-R-channel turn:
      INTEGER:: NZRturn7(1:NReceiver7)=(/5,5,5,5/) !Near -> far sub
      REAL(8):: ZRClen7(1:NReceiver7)=(/0.5,0.5,0.5,0.5/) 
                             !R-coil length,inch, Near -> far sub
      REAL(8):: ZRCdia7(1:NReceiver7)=(/6,6,6,6/)
                             !R-coil diameter,inch, Near -> far sub
      REAL(8):: rPipe7=ODtool7/2.0  !Pipe radius, inch
      !
	END MODULE Tool7_par
!
!***********************************************************************************
c   The End 
!***********************************************************************************
!