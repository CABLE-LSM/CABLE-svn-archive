! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
!    SUBROUTINE DARCY_VG-----------------------------------------------

! Description:
!     Calculates the Darcian fluxes between adjacent soil layers
!     using Van Genuchten formulation

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE darcy_vg(npnts,soil_pts,soil_index                     &
,                 b,ks,sathh                                      &
,                 sthu1,dz1,sthu2,dz2,wflux                       &
,                 dwflux_dsthu1,dwflux_dsthu2                     &
! LOGICAL LTIMER
,ltimer                                                           &
)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Global variables:

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 npnts                                                            &
                      ! IN Number of gridpoints.
,soil_pts             ! IN Number of soil points.


!   Array arguments with intent(IN) :
INTEGER                                                           &
 soil_index(npnts)    ! IN Array of soil points.

REAL                                                              &
 b(npnts,2)                                                       &
                      ! IN Clapp-Hornberger exponent for PSI.
,dz1                                                              &
                      ! IN Thickness of the upper layer (m).
,dz2                                                              &
                      ! IN Thickness of the lower layer (m).
,ks(npnts)                                                        &
                      ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
,sathh(npnts,2)                                                   &
                      ! IN Saturated soil water pressure (m).
,sthu1(npnts)                                                     &
                      ! IN Unfrozen soil moisture content of upper
!                           !    layer as a fraction of saturation.

,sthu2(npnts)         ! IN Unfrozen soil moisture content of lower
!                           !    layer as a fraction of saturation.

LOGICAL ltimer        ! Logical switch for TIMER diags


!   Array arguments with intent(OUT) :
REAL                                                              &
 wflux(npnts)                                                     &
                      ! OUT The flux of water between layers
!                           !     (kg/m2/s).
,dwflux_dsthu1(npnts)                                             &
                      ! OUT The rate of change of the explicit
!                           !     flux with STHU1 (kg/m2/s).
,dwflux_dsthu2(npnts) ! OUT The rate of change of the explicit
!                           !     flux with STHU2 (kg/m2/s).

! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters.

REAL                                                              &
 dthk_dth1,dthk_dth2                                              &
                      ! WORK DTHETAK/DTHETA(1:2).
,dk_dth1,dk_dth2                                                  &
                      ! WORK DK/DTHETA(1:2) (kg/m2/s).
,pd                   ! WORK Hydraulic potential difference (m).

! Local arrays:
REAL                                                              &
 sdum(npnts,2)                                                    &
                      ! WORK bounded values of THETAK.
,bk(npnts)                                                        &
                      ! WORK Clapp-Hornberger exponent at the
!                           !      layer boundary
,theta(npnts,2)                                                   &
                      ! WORK Fractional saturation of the upper
!                           !      and lower layer respectively.
,thetak(npnts)                                                    &
                      ! WORK Fractional saturation at the layer
!                           !      boundary.
,k(npnts)                                                         &
                      ! WORK The hydraulic conductivity between
!                           !      layers (kg/m2/s).
,dk_dthk(npnts)                                                   &
                      ! WORK The rate of change of K with THETAK
!                           !      (kg/m2/s).
,psi(npnts,2)                                                     &
                      ! WORK The soil water suction of the upper
!                           !      and lower layer respectively (m).
,dpsi_dth(npnts,2)                                                &
                      ! WORK The rate of change of PSI with
!                           !      THETA(1:2) (m).
,bracket(npnts,2)                                                 &
                      ! WORK 1-S^(-b-1)
,dbr_dth(npnts,2)     ! WORK The rate of change of BRACKET with
!                           !      THETA(1:2) (m).

REAL                                                              &
 theta_min                                                        &
                  ! Minimum value of THETAK for K calculation.
,theta_max        ! Maximum value of THETAK for K calculation.
PARAMETER (theta_min=0.05, theta_max=0.95)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('DARCY_VG',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Calculate the fractional saturation of the layers
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  theta(i,1)=sthu1(i)
  theta(i,2)=sthu2(i)

! Assume values at layer boundaries are equal to those at the centre of
! the layer above (Thomas Mayr)
  bk(i)=b(i,1)

  sdum(i,1)=MAX(theta(i,1),theta_min)
  sdum(i,1)=MIN(sdum(i,1),theta_max)
  sdum(i,2)=MAX(theta(i,2),theta_min)
  sdum(i,2)=MIN(sdum(i,2),theta_max)
END DO

!-----------------------------------------------------------------------
! Calculate the soil water suction of the layers.
!-----------------------------------------------------------------------
DO n=1,2
!CDIR NODEP
  DO j=1,soil_pts
    i=soil_index(j)
    bracket(i,n)=-1+sdum(i,n)**(-b(i,n)-1)
    psi(i,n)=sathh(i,n)*bracket(i,n)**(b(i,n)/(b(i,n)+1))

!----------------------------------------------------------------------
! To avoid blow-up of implicit increments approximate by piecewise
! linear functions
! (a) for THETA>THETA_MAX  (ensuring that PSI=0 at THETA=1)
! (b) for THETA<THETA_MIN  (extrapolating using DPSI_DTH at THETA_MIN)
!----------------------------------------------------------------------
    IF (theta(i,n).gt.theta_max) THEN
      dpsi_dth(i,n)=-psi(i,n)/(1.0-theta_max)
      psi(i,n)=psi(i,n)                                           &
              +dpsi_dth(i,n)*(MIN(theta(i,n),1.0)-theta_max)
    ELSE
      dbr_dth(i,n)=(-b(i,n)-1)*sdum(i,n)**(-b(i,n)-2)
      dpsi_dth(i,n)=sathh(i,n)*b(i,n)/(b(i,n)+1)                  &
                 *bracket(i,n)**(-1.0/(b(i,n)+1))                 &
                 *dbr_dth(i,n)
      IF (theta(i,n).lt.theta_min) THEN
        psi(i,n)=psi(i,n)                                         &
                +dpsi_dth(i,n)*(MAX(theta(i,n),0.0)-theta_min)
      END IF
    END IF

    IF ((theta(i,n).gt.1.0).OR.(theta(i,n).lt.0.0)) THEN
      dpsi_dth(i,n)=0.0
      psi(i,n)=MAX(psi(i,n),0.0)
    END IF

  END DO
END DO

!-----------------------------------------------------------------------
! Estimate the fractional saturation at the layer boundary by
! interpolating the soil moisture.
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  thetak(i)=(dz2*theta(i,1)+dz1*theta(i,2))/(dz2+dz1)
  bk(i)=(dz2*b(i,1)+dz1*b(i,2))/(dz2+dz1)
END DO
dthk_dth1=dz2/(dz1+dz2)
dthk_dth2=dz1/(dz2+dz1)

!-----------------------------------------------------------------------
! Calculate the hydraulic conductivities for transport between layers.
!-----------------------------------------------------------------------
! DEPENDS ON: hyd_con_vg
CALL hyd_con_vg(npnts,soil_pts,soil_index,bk,ks,thetak,k          &
,             dk_dthk,ltimer)

!-----------------------------------------------------------------------
! Calculate the Darcian flux from the upper to the lower layer.
!-----------------------------------------------------------------------
DO j=1,soil_pts
  i=soil_index(j)
  pd=(2.0*(psi(i,2)-psi(i,1))/(dz2+dz1)+1)
  wflux(i)=k(i)*pd

!-----------------------------------------------------------------------
! Calculate the rate of change of WFLUX with respect to the STHU1 and
! STHU2.
!-----------------------------------------------------------------------
  dk_dth1=dk_dthk(i)*dthk_dth1
  dk_dth2=dk_dthk(i)*dthk_dth2
  dwflux_dsthu1(i)=dk_dth1*pd-2*k(i)*dpsi_dth(i,1)/(dz1+dz2)
  dwflux_dsthu2(i)=dk_dth2*pd+2*k(i)*dpsi_dth(i,2)/(dz1+dz2)
END DO

IF (lhook) CALL dr_hook('DARCY_VG',zhook_out,zhook_handle)
RETURN
END SUBROUTINE darcy_vg
