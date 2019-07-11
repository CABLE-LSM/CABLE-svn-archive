#if !defined(UM_JULES)
MODULE imogen_drive_vars

  IMPLICIT NONE

!-----------------------------------------------------------------
! Fine temperal resolution year of climatology (to be used to
! drive JULES)
!-----------------------------------------------------------------
    REAL, ALLOCATABLE ::                                          &
      T_OUT(:,:,:,:),                                             &
                  ! Calculated temperature (K)
      CONV_RAIN_OUT(:,:,:,:),                                     &
                  ! Calculated convective rainfall (mm/day)
      CONV_SNOW_OUT(:,:,:,:),                                     &
                  ! Calculated convective snowfall (mm/day)
      LS_RAIN_OUT(:,:,:,:),                                       &
                  ! Calculated large scale rainfall (mm/day)
      LS_SNOW_OUT(:,:,:,:),                                       &
                  ! Calculated large scale snowfall (mm/day)
      QHUM_OUT(:,:,:,:),                                          &
                  ! Calculated humidity
      WIND_OUT(:,:,:,:),                                          &
                  ! Calculated wind  (m/s)
      PSTAR_OUT(:,:,:,:),                                         &
                  ! Calculated pressure (Pa)
      SW_OUT(:,:,:,:),                                            &
                  ! Calculated shortwave radiation
      LW_OUT(:,:,:,:)
                  ! Calculated longwave radiation

END MODULE imogen_drive_vars
#endif
