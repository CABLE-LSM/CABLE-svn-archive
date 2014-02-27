#if !defined(UM_JULES)
!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/var/p_s_parms.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

MODULE p_s_parms

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing plant and soil variables (plus a few others)
!-----------------------------------------------------------------------------

  REAL, ALLOCATABLE ::                                                        &
    ALBSOIL(:),                                                               &
      !  Soil albedo
    B(:,:),                                                                   &
      !  Exponent for soil moisture characteristic functions
      !    Clapp-Hornberger model: b is the Clapp-Hornberger exponent
      !    van Genuchten model: b=1/(n-1)  (metres)
    CATCH(:,:),                                                               &
      !  Surface/canopy water capacity of snow-free land tiles (kg/m2)
    CATCH_SNOW(:,:),                                                          &
      !  Snow interception capacity (kg/m2)
    COSZ(:),                                                                  &
      !  Cosine of the zenith angle
    SATHH(:,:),                                                               &
      !  Parameter for soil moisture characteristic functions
      !    Clapp-Hornberger model: sathh is the saturated soil water pressure (m)
      !    van Genuchten model: sathh=1/alpha
    HCON(:,:),                                                                &
      !  Soil thermal conductivity (W/m/K)
    SMVCCL(:,:),                                                              &
      !  Critical volumetric SMC (cubic m per cubic m of soil)
    SMVCST(:,:),                                                              &
      !  Volumetric saturation point (m3/m3 of soil)
    SMVCWT(:,:),                                                              &
      !  Volumetric wilting point (cubic m per cubic m of soil)
    HCAP(:,:),                                                                &
      !  Soil heat capacity (J/K/m3)
    SATCON(:,:),                                                              &
      !  Saturated hydraulic conductivity (kg/m2/s)
    INFIL_TILE(:,:),                                                          &
      !  Maximum possible surface infiltration for tiles (kg/m2/s)
    Z0_TILE(:,:),                                                             &
      !  Surface roughness on tiles (m).
    STHU(:,:),                                                                &
      !  Unfrozen soil moisture content of the layers as a fraction of saturation.
    STHF(:,:),                                                                &
      !  Frozen soil moisture content of the layers as a fraction of saturation.
    SOIL_CLAY(:,:)
      !  Soil clay fraction

END MODULE p_s_parms
#endif
