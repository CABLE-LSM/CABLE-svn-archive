#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/dump/required_vars_for_configuration.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE required_vars_for_configuration(nvars, identifiers, total_snow,    &
                                                          include_imogen_vars)

  USE switches, ONLY : can_model, l_imogen, l_phenol, l_spec_albedo, l_top,   &
                       l_triffid, l_veg_compete

  USE ancil_info, ONLY : nsmax

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns the identifiers of the prognostic variables for the current
!   configuration of the model
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(OUT) :: nvars  ! The number of variables
  CHARACTER(len=*), INTENT(OUT) :: identifiers(:)
                                 ! The model identifiers of the required
                                 ! variables
  LOGICAL, INTENT(IN), OPTIONAL :: total_snow
                                 ! T - return only the variables required if
                                 !     snow layer values are to be set by the
                                 !     model
                                 ! F - return all snow variables required
                                 !     by the current configuration
  LOGICAL, INTENT(IN), OPTIONAL :: include_imogen_vars
                                 ! T - include IMOGEN prognostics in the
                                 !     returned list if they are needed
                                 ! F - do not include IMOGEN prognostics in
                                 !     the returned list, even if they are
                                 !     needed

! Work variables
  LOGICAL :: total_snow_local  ! Local version of total_snow
                               ! Defaults to FALSE if not present (i.e.
                               ! return all required snow vars)
  LOGICAL :: inc_imogen_local  ! Local version of include_imogen_vars
                               ! Defaults to TRUE if not present, i.e.
                               ! return IMOGEN vars if they are required


!-----------------------------------------------------------------------------

! Check for presence of optional arguments
  total_snow_local = .FALSE.
  inc_imogen_local = .TRUE.
  IF ( PRESENT(total_snow) ) total_snow_local = total_snow
  IF ( PRESENT(include_imogen_vars) ) inc_imogen_local = include_imogen_vars


! First, set up the array with variables that are required with every
! configuration
! Note that for total_snow = T, snow_tile is the only required snow variable
  nvars = 7
  identifiers(1:nvars) = (/ 'canopy    ', 'cs        ', 'gs        ',         &
                            'snow_tile ', 'sthuf     ', 't_soil    ',         &
                            'tstar_tile' /)

!-------------------------------------------------------------------------------
!   Variables that are needed only if certain options set.
!-------------------------------------------------------------------------------
  IF ( l_phenol ) THEN
! With phenology on, LAI is prognostic
    nvars = nvars + 1
    identifiers(nvars) = 'lai'
  ENDIF

  IF ( l_triffid ) THEN
! With TRIFFID on, canopy height is prognostic
    nvars = nvars + 1
    identifiers(nvars) = 'canht'
  ENDIF

  IF ( l_veg_compete ) THEN
! With competing vegetation on, frac is prognostic
    nvars = nvars + 1
    identifiers(nvars) = 'frac'
  END IF

! TOPMODEL variables.
  IF ( l_top ) THEN
! Wetness in deep layer.
    nvars = nvars + 1
    identifiers(nvars) = 'sthzw'

! Depth to water table.
    nvars = nvars + 1
    identifiers(nvars) = 'zw'
  ENDIF

  IF ( l_spec_albedo ) THEN
    nvars = nvars + 1
    identifiers(nvars) = 'rgrain'
  ENDIF

!-----------------------------------------------------------------------------
! Add IMOGEN variables if IMOGEN is on and they are requested
!-----------------------------------------------------------------------------
  IF ( l_imogen .AND. inc_imogen_local ) THEN
    nvars = nvars + 1
    identifiers(nvars) = 'co2_ppmv'

    nvars = nvars + 1
    identifiers(nvars) = 'co2_change_ppmv'

    nvars = nvars + 1
    identifiers(nvars) = 'dtemp_o'

    nvars = nvars + 1
    identifiers(nvars) = 'fa_ocean'

    nvars = nvars + 1
    identifiers(nvars) = 'seed_rain'

! This should possibly be under an l_triffid switch, but is here for now
! so as not to change behaviour for non-TRIFFID runs
    nvars = nvars + 1
    identifiers(nvars) = 'cv'
  ENDIF

!-----------------------------------------------------------------------------
! Work out what snow variables are required
! If total_snow = T, only snow_tile is required, and is always required
! so if already in the list
! We just need to add the other required snow variables based on the scheme
! enabled if total_snow = F
!-----------------------------------------------------------------------------
  IF ( .NOT. total_snow_local ) THEN
! Snow variables not specifically for the multi-layer model
    nvars = nvars + 1
    identifiers(nvars) = 'rho_snow'

    nvars = nvars + 1
    identifiers(nvars) = 'snow_depth'

    IF ( can_model == 4 ) THEN
      nvars = nvars + 1
      identifiers(nvars) = 'snow_grnd'
    ENDIF

! Variables for the multi-layer snow model.
    IF ( nsmax > 0 ) THEN
      nvars = nvars + 1
      identifiers(nvars) = 'nsnow'

      nvars = nvars + 1
      identifiers(nvars) = 'snow_ds'

      nvars = nvars + 1
      identifiers(nvars) = 'snow_ice'

      nvars = nvars + 1
      identifiers(nvars) = 'snow_liq'

      nvars = nvars + 1
      identifiers(nvars) = 'tsnow'

      IF ( l_spec_albedo ) THEN
        nvars = nvars + 1
        identifiers(nvars) = 'rgrainl'
      ENDIF
    ENDIF   ! nsmax

  END IF  ! total_snow

  RETURN

END SUBROUTINE required_vars_for_configuration
#endif
