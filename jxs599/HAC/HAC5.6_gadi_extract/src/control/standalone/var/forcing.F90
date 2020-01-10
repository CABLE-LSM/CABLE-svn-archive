#if !defined(UM_JULES)
! Module containing all of the driving (atmospheric forcing) variables.

MODULE forcing

IMPLICIT NONE

!-------------------------------------------------------------------------------

! The forcing variables.
REAL, ALLOCATABLE ::                                                          &
  qw_1_ij(:,:),                                                               &
                      !  Total water content (Kg/Kg)
  tl_1_ij(:,:),                                                               &
                      ! Ice/liquid water temperature (k)
  u_0_ij(:,:),                                                                &
                      ! W'ly component of surface current (m/s)
  v_0_ij(:,:),                                                                &
                      ! S'ly component of surface current (m/s)
  u_1_ij(:,:),                                                                &
                      ! W'ly wind component (m/s)
  v_1_ij(:,:),                                                                &
                      ! S'ly wind component (m/s)
  pstar_ij(:,:),                                                              &
                      ! Surface pressure (Pascals)
  ls_rain_ij(:,:),                                                            &
                      ! Large-scale rain (kg/m2/s)
  con_rain_ij(:,:),                                                           &
                      ! Convective rain (kg/m2/s)
  ls_snow_ij(:,:),                                                            &
                      ! Large-scale snowfall (kg/m2/s)
  con_snow_ij(:,:),                                                           &
                      ! Convective snowfall (kg/m2/s)
  sw_down_ij(:,:),                                                            &
                      ! Surface downward SW radiation (W/m2)
  lw_down_ij(:,:),                                                            &
                      ! Surface downward LW radiation (W/m2)
  diurnal_temperature_range_ij(:,:)
                      ! diurnal temperature range (K), used when
                      !l_dailydisagg=T

! Variables that aid in the calculation of the actual forcing variables
REAL, ALLOCATABLE ::                                                          &
  diff_rad_ij(:,:)       ! Input diffuse radiation (W/m2)

!-------------------------------------------------------------------------------

END MODULE forcing
#endif
