#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/initial_conditions/topmodel_init.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE topmodel_init()

  USE ancil_info, ONLY : land_pts, sm_levels, soil_pts, soil_index

  USE p_s_parms, ONLY : b, satcon, sthf, sthu

  USE soil_param, ONLY : dzsoil

  USE top_pdm, ONLY : fexp, gamtot, ti_mean, ti_sig, zw, fsat, fwetl, qbase

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Finish initialising TOPMODEL by calculating surface saturated fraction
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  LOGICAL, PARAMETER :: l_gamtot = .TRUE.  !  switch for calculation of gamtot

  INTEGER :: i, j, n  ! Index variables

  REAL :: qbase_l(land_pts,sm_levels+1) !  Base flow from each layer (kg/m2/s).
  REAL :: top_crit(land_pts)            !  Critical topographic index required
                                        !    to calculate the surface saturation
                                        !    fraction.
  REAL :: zdepth(0:sm_levels)           !  Lower soil layer boundary depth (m).
  REAL :: wutot(land_pts)               !  UNFROZEN to TOTAL fraction at ZW.
  REAL :: ksz(land_pts,0:sm_levels)     !  Saturated hydraulic conductivity for
                                        !    each layer (kg/m2/s).

!-------------------------------------------------------------------------------
  zdepth(0) = 0.0
  DO n = 1,sm_levels
    zdepth(n) = zdepth(n - 1) + dzsoil(n)
  ENDDO

! Set value that is retained at non-soil points.
  fsat(:) = 0.0
  fwetl(:) = 0.0

  IF ( soil_pts /= 0 ) THEN
    DO j = 1,soil_pts
      i = soil_index(j)
      DO n = 0,sm_levels
        ksz(i,n) = satcon(i,n)
      ENDDO
    ENDDO

!   We need top_crit to calculate fsat - get this from calc_baseflow.
    CALL calc_baseflow_jules(soil_pts, soil_index, land_pts, sm_levels,       &
                             zdepth, ksz, b, fexp, ti_mean, zw, sthf, sthu,   &
                             wutot, top_crit, qbase, qbase_l)

!   Call calc_fsat with 1st argument (l_gamtot)=.FALSE. so as get fsat.
    CALL calc_fsat(l_gamtot, soil_pts, soil_index, land_pts, ti_mean, ti_sig, &
                   wutot, top_crit, gamtot, fsat, fwetl)
  ENDIF  !  soil_pts

  RETURN

END SUBROUTINE topmodel_init
#endif
