!
!**********************************************************************
!
!   File name  :  COND_xzRotation.f
!   SUBROUTINE :  Cfb_xzRotation(),Cbf_xzRotation()
!   Purpose    :  Rotation of apparent conductivity tensors
!                 in the X-Z plane
!   Date       :  09/02/2009
!   Update     :  09/24/2009
!   Author     :  Hou, Junsheng
!
!**********************************************************************
!
!    SUBROUTINE :  Cfb_xzRotation(dip,Cf,Cb)
!    Purpose    :  Rotation of (apparent conductivity) tensor
!                  from formation system (Cf) to borehole system (Cb)
!                  (Clockwise rotation)
!    transformation (or rotation) matrix
!
!**********************************************************************
!
      SUBROUTINE Cfb_xzRotation(dip,Cf,Cb)
      IMPLICIT NONE
      REAL(8),INTENT(IN)    :: dip     !(unit: degree)
      REAL(8),INTENT(IN) :: Cf(3,3)
      REAL(8),INTENT(OUT):: Cb(3,3)
!local:
      REAL(8) RFB(3,3)
!
      Cb=0.0D0
	RFB(1,1)=COSD(dip)
	RFB(1,2)=0.D0
	RFB(1,3)=SIND(dip)
	RFB(2,1)=0.D0
	RFB(2,2)=1.D0
	RFB(2,3)=0.D0
	RFB(3,1)=-SIND(dip)
	RFB(3,2)=0.D0
	RFB(3,3)=COSD(dip)
!
	Cb=MATMUL(MATMUL(RFB,Cf),TRANSPOSE(RFB))
!
      RETURN
      END SUBROUTINE CFB_xzRotation
!
!**********************************************************************
!
!    SUBROUTINE: Cbf_xzRotation(dip,Cb,Cf)
!    Purpose   : Rotation of apparent conductivity tensor
!                from borehole system (Cb) to formation system (Cf)
!                (Counter-clockwise rotation)
!
!**********************************************************************
!
      SUBROUTINE Cbf_xzRotation(dip,Cb,Cf)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: dip       !(unit: degree)
      REAL(8),INTENT(IN) :: Cb(3,3)
      REAL(8),INTENT(OUT):: Cf(3,3)
!
      REAL(8) RBF(3,3)
!
      Cf=0.0D0
	  RBF(1,1)=COSD(dip)
	  RBF(1,2)=0.D0
	  RBF(1,3)=-SIND(dip)
	  RBF(2,1)=0.D0
	  RBF(2,2)=1.D0
	  RBF(2,3)=0.D0
	  RBF(3,1)=SIND(dip)
	  RBF(3,2)=0.D0
	  RBF(3,3)=COSD(dip)
	Cf=MATMUL(MATMUL(RBF,Cb),TRANSPOSE(RBF))
!
      RETURN
      END SUBROUTINE Cbf_xzRotation
!
!**********************************************************************
!