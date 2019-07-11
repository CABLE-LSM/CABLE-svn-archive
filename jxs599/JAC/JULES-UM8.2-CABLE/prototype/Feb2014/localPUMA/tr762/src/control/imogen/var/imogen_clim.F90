#if !defined(UM_JULES)
MODULE imogen_clim

  IMPLICIT NONE


! Maximum and minimum latitude and longitude for the climatology
! variables
  REAL :: LATMIN_CLIM,LATMAX_CLIM,LONGMIN_CLIM,LONGMAX_CLIM

!-----------------------------------------------------------------
! Driving "control" climatology
!-----------------------------------------------------------------
  REAL, ALLOCATABLE ::                                            &
    T_CLIM(:,:),                                                  &
                ! Control climate temperature
    RAINFALL_CLIM(:,:),                                           &
                ! Control climate rainfall
    SNOWFALL_CLIM(:,:),                                           &
                ! Control climate snowfall
    RH15M_CLIM(:,:),                                              &
                ! Control climate relative humidity
    UWIND_CLIM(:,:),                                              &
                ! Control climate u-wind
    VWIND_CLIM(:,:),                                              &
                ! Control climate v-wind
    DTEMP_CLIM(:,:),                                              &
                ! Control climate diurnal Temp
    PSTAR_HA_CLIM(:,:),                                           &
                ! Control climate pressure
    SW_CLIM(:,:),                                                 &
                ! Control climate shortwave radiation
    LW_CLIM(:,:),                                                 &
                ! Control climate longwave radiation
    F_WET_CLIM(:,:)
                ! Control climate fraction we

  REAL, ALLOCATABLE ::                                            &
    LAT(:),                                                       &
                  ! WORK Latitudinal position of land
    LONG(:),                                                      &
                  ! WORK Longitudinal position of land
    DCTOT(:)

END MODULE imogen_clim
#endif
