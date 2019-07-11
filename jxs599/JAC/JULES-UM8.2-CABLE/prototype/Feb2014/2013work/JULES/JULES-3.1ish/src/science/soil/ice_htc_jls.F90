! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ICE_HTC------------------------------------------------

! Description:
!     Updates deep soil temperatures for ice. No external subroutines
!     are called.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
SUBROUTINE ice_htc (                                              &
 npnts,nshyd,lice_pts,lice_index,dz                               &
,surf_ht_flux,timestep                                            &
,tsoil                                                            &
,ltimer                                                           &
)

USE snow_param, ONLY :                                            &
!      imported scalars with intent(in)
 snow_hcap,snow_hcon

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER                                                           &
 lice_pts                                                         &
                      ! IN Number of land ice points.
,npnts                                                            &
                      ! IN Number of gridpoints.
,nshyd                ! IN Number of soil moisture levels.

REAL                                                              &
 timestep             ! IN Model timestep (s).


!   Array arguments with intent(IN) :
INTEGER                                                           &
 lice_index(npnts)    ! IN Array of ice points.

REAL                                                              &
 dz(nshyd)                                                        &
                      ! IN Thicknesses of the soil layers (m).
,surf_ht_flux(npnts)  ! IN Net downward surface heat flux (W/m2).

LOGICAL ltimer        ! Logical switch for TIMER diags

!   Array arguments with intent(INOUT) :
REAL                                                              &
 tsoil(npnts,nshyd)   ! INOUT Sub-surface temperatures (K).

! Local scalars:
INTEGER                                                           &
 i,j,n                ! WORK Loop counters.

! Local arrays:
REAL                                                              &
 h_flux(npnts,0:nshyd)! WORK The fluxes of heat between layers
!                           !      (W/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ICE_HTC',zhook_in,zhook_handle)

!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------
DO n=1,nshyd-1
  DO j=1,lice_pts
    i=lice_index(j)
    h_flux(i,n)=-snow_hcon*2.0*(tsoil(i,n+1)-tsoil(i,n))          &
                             /(dz(n+1)+dz(n))
  END DO
END DO

!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
DO j=1,lice_pts
  i=lice_index(j)
  h_flux(i,nshyd)=0.0
  h_flux(i,0)=surf_ht_flux(i)
END DO


!--------------------------------------------------------------------
! Update the sub-surface temperatures
!--------------------------------------------------------------------
DO n=1,nshyd
!CDIR NOVECTOR
! CDIR$ IVDEP here would force vectorization but changes results!
  DO j=1,lice_pts
    i=lice_index(j)

    tsoil(i,n)=tsoil(i,n)                                         &
     +1.0/(snow_hcap*dz(n))*(h_flux(i,n-1)                        &
     -h_flux(i,n))*timestep

  END DO
END DO

IF (lhook) CALL dr_hook('ICE_HTC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ice_htc
