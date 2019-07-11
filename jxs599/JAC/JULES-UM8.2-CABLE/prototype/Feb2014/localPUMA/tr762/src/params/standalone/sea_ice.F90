#if !defined(UM_JULES)
! Module containing all of the sea-ice parameter variables

  MODULE sea_ice

  REAL :: SW_alpham
!                    ! if l_ssice_albedo then
!                    !   Albedo of snow on sea-ice at melting point (TM)
!                    ! if .not. l_ssice_albedo then
!                    !   Albedo of sea-ice at melting point (TM) 

  REAL :: SW_alphac
!                    ! if l_ssice_albedo then
!                    !    Albedo of snow on sea-ice at and below TM-DTICE
!                    ! if .not. l_ssice_albedo then
!                    !    Albedo of sea-ice at and below TM-DTICE

  REAL :: SW_alphab
!                    ! if l_ssice_albedo then
!                    !   Albedo of snow-free sea-ice
!                    ! if .not. l_ssice_albedo then
!                    !   Not used

  REAL :: SW_dtice
!                    ! if l_ssice_albedo then
!                    !    Temperature range in which albedo of 
!                    !      snow varies between its limits
!                    ! if .not. l_ssice_albedo then
!                    !    Temperature range in which albedo of 
!                    !      sea-ice varies between its limits

!-----------------------------------------------------------------------
! New variables for sea-ice albedo calculations
!-----------------------------------------------------------------------
  REAL, PARAMETER ::                                                &
 & DT_BARE       = 1.0                                              &
      ! Temperature range below TM over which meltponds form if
      ! l_sice_meltponds and l_ssice_albedo
 &,DALB_BARE_WET = -0.075                                           &
      ! Increment to albedo for each degree temperature rises above
      ! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
 &,PEN_RAD_FRAC  = 0.17                                             &
      ! Fraction of SW radiation that penetrates seaice and scatters
      ! back causing an addition to the albedo. Only active if
      ! l_ssice_albedo and l_sice_scattering
 &,BETA          = 0.4
      ! attenutation parameter for SW in seaice which controls the
      ! additional albedo due to internal scattering

  END MODULE sea_ice
#endif
