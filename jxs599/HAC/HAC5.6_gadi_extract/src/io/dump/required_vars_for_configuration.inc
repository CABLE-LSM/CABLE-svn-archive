#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE required_vars_for_configuration(nvars, identifiers,                &
                                           nvars_from_ancil, vars_from_ancil, &
                                           l_output_mode,                     &
                                           total_snow, include_imogen_vars,   &
                                           read_or_write_dump)

USE update_mod, ONLY: l_imogen, have_prescribed_sthuf

USE jules_hydrology_mod, ONLY: l_top

USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow

USE jules_soil_biogeochem_mod, ONLY: soil_model_ecosse, soil_model_rothc,     &
                                      soil_bgc_model

USE jules_soil_ecosse_mod, ONLY: l_soil_N

USE jules_vegetation_mod, ONLY: can_model, can_rad_mod,                       &
                                l_triffid, l_phenol, l_veg_compete,           &
                                l_crop, l_irrig_dmd, l_landuse,               &
                                l_nitrogen, l_prescsow, l_trif_crop

USE jules_snow_mod, ONLY: nsmax

USE jules_surface_mod, ONLY: l_elev_land_ice

USE jules_soil_mod, ONLY: l_bedrock, sm_levels, l_tile_soil

USE fire_mod,     ONLY: fire_cntl, l_fire

USE metstats_mod, ONLY: metstats_flag, l_metstats

USE update_mod, ONLY: l_daily_disagg, precip_disagg_method,                   &
                      prescribed_sthuf_levels

USE jules_rivers_mod, ONLY: l_rivers, i_river_vn, rivers_trip, rivers_rfm

USE logging_mod, ONLY: log_warn
  
USE imogen_run, ONLY: land_feed_ch4

USE prognostics, ONLY:                                                        &
  l_broadcast_soilt


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns the identifiers of the prognostic variables for the current
!   configuration of the model
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(OUT) :: nvars  ! The number of variables
CHARACTER(LEN=*), INTENT(OUT) :: identifiers(:)
                               ! The model identifiers of the required
                               ! variables

INTEGER, INTENT(OUT)          :: nvars_from_ancil  !Number of ancil variables
CHARACTER(LEN=*), INTENT(OUT) :: vars_from_ancil(:)
                               ! The model identifiers of the required
                               ! ancil variables

LOGICAL, INTENT(IN)           :: l_output_mode
                               ! True if the intention of the call to this
                               ! routine is to get the variables to output
                               ! False if for input
                               ! This enables the broadcasting of non-soil
                               ! tiled initial conditions for soil tiling
                               ! runs

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
LOGICAL, INTENT(IN), OPTIONAL :: read_or_write_dump
                               ! T - include all variables that are needed
                               !     when reading or writing to dump files
                               ! F - don't include variables that are only
                               !     inputted or outputted in dump files

! Work variables
LOGICAL :: total_snow_local  ! Local version of total_snow
                             ! Defaults to FALSE if not present (i.e.
                             ! return all required snow vars)
LOGICAL :: inc_imogen_local  ! Local version of include_imogen_vars
                             ! Defaults to TRUE if not present, i.e.
                             ! return IMOGEN vars if they are required
LOGICAL :: read_or_write_dump_local
                             ! Local version of read_or_write_dump
                             ! Defaults to TRUE if not present
LOGICAL :: l_append_soilt
                             ! Controls whether vars are appended with _soilt
INTEGER :: num_prescribed_sthuf_levels
                             ! Number of levels in the prescibed sthuf field
!-----------------------------------------------------------------------------

! Check for presence of optional arguments
total_snow_local          = .FALSE.
inc_imogen_local          = .TRUE.
read_or_write_dump_local  = .TRUE.

! Decide whether we need to append _soilt to variables
l_append_soilt = .FALSE.
IF (l_tile_soil) THEN
  IF (l_output_mode) THEN
    l_append_soilt = .TRUE.
  ELSE
    IF ( .NOT. l_broadcast_soilt) THEN
      l_append_soilt = .TRUE.
    END IF
  END IF
END IF

IF ( PRESENT(total_snow) )                                                    &
  total_snow_local         = total_snow
IF ( PRESENT(include_imogen_vars) )                                           &
  inc_imogen_local         = include_imogen_vars
IF ( PRESENT(read_or_write_dump) )                                            &
  read_or_write_dump_local = read_or_write_dump

! First, set up the array with variables that are required with every
! configuration.
nvars = 0

nvars = nvars + 1
identifiers(nvars) = 'canopy'

nvars = nvars + 1
identifiers(nvars) = 'cs'
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

! Note that for total_snow = T, snow_surft is the only required snow variable
nvars = nvars + 1
identifiers(nvars) = 'snow_tile'

nvars = nvars + 1
identifiers(nvars) = 't_soil'
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'tstar_tile'

!-----------------------------------------------------------------------------
! Variables that are needed only if certain options set.
!-----------------------------------------------------------------------------
IF ( can_rad_mod == 1 ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'gs'
END IF

IF (ALLOCATED(prescribed_sthuf_levels)) THEN
  num_prescribed_sthuf_levels = SIZE(prescribed_sthuf_levels)
ELSE
  num_prescribed_sthuf_levels = 0
END IF

IF ( ( .NOT. have_prescribed_sthuf ) .OR.                                     &
     ( num_prescribed_sthuf_levels < sm_levels) ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'sthuf'

  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF
END IF

IF ( l_phenol ) THEN
  ! With phenology on, LAI is prognostic for all PFTs
  nvars = nvars + 1
  identifiers(nvars) = 'lai'
ELSE IF ( l_crop ) THEN
  ! Otherwise, if the crop scheme is on, LAI is prognostic for crop PFTs only
  nvars = nvars + 1
  identifiers(nvars) = 'croplai'
END IF

IF ( l_triffid ) THEN
  ! With TRIFFID on, canopy height is prognostic for all PFTs
  nvars = nvars + 1
  identifiers(nvars) = 'canht'
ELSE IF ( l_crop ) THEN
  ! Otherwise, if the crop scheme is on, canopy height is prognostic for crop
  ! PFTs only
  nvars = nvars + 1
  identifiers(nvars) = 'cropcanht'
END IF

IF ( l_triffid .AND. l_landuse ) THEN
  ! With TRIFFID on,  agricultural frac and WP pools are required for LUC
  nvars = nvars + 1
  identifiers(nvars) = 'frac_agr_prev'
  nvars = nvars + 1
  identifiers(nvars) = 'wood_prod_fast'
  nvars = nvars + 1
  identifiers(nvars) = 'wood_prod_med'
  nvars = nvars + 1
  identifiers(nvars) = 'wood_prod_slow'
  IF (l_trif_crop) THEN
    nvars = nvars + 1
    identifiers(nvars) = 'frac_past_prev' 
  END IF
END IF

IF ( soil_bgc_model == soil_model_rothc .AND. l_nitrogen ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'n_inorg'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

  nvars = nvars + 1
  identifiers(nvars) = 'ns'
END IF

IF ( soil_bgc_model == soil_model_ecosse .AND. l_soil_N ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'ns'

  nvars = nvars + 1
  identifiers(nvars) = 'n_amm'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

  nvars = nvars + 1
  identifiers(nvars) = 'n_nit'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

END IF

! TOPMODEL variables.
IF ( l_top ) THEN
  ! Wetness in deep layer.
  nvars = nvars + 1
  identifiers(nvars) = 'sthzw'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

  ! Depth to water table.
  nvars = nvars + 1
  identifiers(nvars) = 'zw'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

END IF

IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'rgrain'
END IF

! Additional crop scheme prognostics, required if crop scheme is on
IF ( l_crop ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'cropdvi'

  nvars = nvars + 1
  identifiers(nvars) = 'croprootc'

  nvars = nvars + 1
  identifiers(nvars) = 'cropharvc'

  nvars = nvars + 1
  identifiers(nvars) = 'cropreservec'
END IF

IF ( l_irrig_dmd ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'sthu_irr'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF
END IF

! Additional deep soil temperature if bedrock is on
IF ( l_bedrock ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'tsoil_deep'
END IF

! River storage if river routing on and using TRIP
IF ( l_rivers .AND. ( i_river_vn == rivers_trip ) ) THEN
  IF ( read_or_write_dump_local ) THEN
    nvars = nvars + 1
    identifiers(nvars) = 'rivers_sto_rp'
  ELSE
    CALL log_warn("required_vars_for_configuration",                          &
                  "rivers_sto_rp will be initialised to zero.")      
  END IF
END IF

! Surface and subsurface stores and flows if river routing on and using RFM
IF ( l_rivers .AND. (i_river_vn == rivers_rfm) ) THEN 
  IF ( read_or_write_dump_local ) THEN 
    nvars = nvars + 1
    identifiers(nvars) = 'rfm_surfstore_rp'
  
    nvars = nvars + 1
    identifiers(nvars) = 'rfm_substore_rp'

    nvars = nvars + 1
    identifiers(nvars) = 'rfm_flowin_rp'

    nvars = nvars + 1
    identifiers(nvars) = 'rfm_bflowin_rp'
  ELSE
    CALL log_warn("required_vars_for_configuration",                          &
                  "RFM river prognostics will be initialised to zero.")
  END IF
END IF

!-----------------------------------------------------------------------------
! Add metstats prognostic variables if switched on
!-----------------------------------------------------------------------------
IF ( l_metstats ) THEN
  IF (metstats_flag%temp_max_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'temp_max_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'temp_max_00h_r'
  END IF

  IF (metstats_flag%temp_ave_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'temp_ave_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'temp_ave_00h_r'
  END IF

  IF (metstats_flag%temp_pnt_12h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'temp_pnt_12h'
  END IF

  IF (metstats_flag%prec_tot_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'prec_tot_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'prec_tot_00h_r'
  END IF

  IF (metstats_flag%prec_tot_12h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'prec_tot_12h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'prec_tot_12h_r'
  END IF

  IF (metstats_flag%rhum_min_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'rhum_min_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'rhum_min_00h_r'
  END IF

  IF (metstats_flag%rhum_pnt_12h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'rhum_pnt_12h'
  END IF

  IF (metstats_flag%dewp_ave_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'dewp_ave_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'dewp_ave_00h_r'
  END IF

  IF (metstats_flag%wind_ave_00h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'wind_ave_00h'
    nvars              =  nvars + 1
    identifiers(nvars) = 'wind_ave_00h_r'
  END IF

  IF (metstats_flag%wind_pnt_12h) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'wind_pnt_12h'
  END IF
END IF

!-----------------------------------------------------------------------------
! Add fire prognostic variables if switched on
!-----------------------------------------------------------------------------
IF ( l_fire ) THEN
  IF (fire_cntl%mcarthur%flag) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_mcarthur_r_dr'
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_mcarthur_n_dr'
  END IF

  IF (fire_cntl%canadian%flag) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_canadian_ffmc'
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_canadian_ffmc_mois'
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_canadian_dmc'
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_canadian_dc'
  END IF

  IF (fire_cntl%nesterov%flag) THEN
    nvars              =  nvars + 1
    identifiers(nvars) = 'fire_nesterov'
  END IF
END IF !Fire

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

  IF (land_feed_ch4) THEN
    nvars = nvars + 1
    identifiers(nvars) = 'ch4_ppbv'
  END IF
END IF

IF ( l_daily_disagg .AND. read_or_write_dump_local .AND.                      &
     (precip_disagg_method > 1) ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'seed_rain'
END IF

!-----------------------------------------------------------------------------
! Work out what snow variables are required
! If total_snow = T, only snow_surft is required, and is always required
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
  END IF

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

    IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
      nvars = nvars + 1
      identifiers(nvars) = 'rgrainl'
    END IF
  END IF   ! nsmax

END IF  ! total_snow

IF ( l_elev_land_ice) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'tsurf_elev_surft'
END IF

!-----------------------------------------------------------------------------
! Work out which ancillary fields need to be in the dump file
! The logic here is slightly different from normal prognostic variables
!
! There are several cases to consider:
! -The var is not needed
! -The var should be populated by ancil and written to the dump (ie write only
!  as far as this module is concerned)
! -The var should be both read from the initial condition namelist and written
!  to the dump
!
! To get the information about which ancils will be read from file, we use
! ancil_dump_read%xyz:  Defaulting to false, set in init_ancillaries_mod
!-----------------------------------------------------------------------------

  !Initialise
nvars_from_ancil = 0

!-----------------------------------------------------------------------------
!Set up read/write details for frac ancil
!This is a special case: with competing vegetation on, frac is prognostic

! frac should always be written to the dump
nvars = nvars + 1
identifiers(nvars) = 'frac'

IF ( l_veg_compete ) THEN
  ! With competing veg on, frac is a prognostic and should always be read
  ! as an initial condition
  ancil_dump_read%frac  = .TRUE.
END IF

IF ( .NOT. ancil_dump_read%frac ) THEN
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'frac'
END IF

!---------------------------------------------------------------------------
! Set up read/write details for soil properites ancils
! Most soil properties are always needed.
nvars = nvars + 1
identifiers(nvars) = 'b      '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'sathh  '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'satcon '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'sm_sat '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'sm_crit'
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'sm_wilt'
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'hcap   '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'hcon   '
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

nvars = nvars + 1
identifiers(nvars) = 'albsoil'
IF (l_append_soilt) THEN
  identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
END IF

IF ( soil_bgc_model == soil_model_rothc .OR.                                  &
     soil_bgc_model == soil_model_ecosse) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'clay   '
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF
END IF

IF ( soil_bgc_model == soil_model_ecosse ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'soil_ph'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF
END IF

IF ( .NOT. ancil_dump_read%soil_props ) THEN
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'b      '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'sathh  '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'satcon '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'sm_sat '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'sm_crit'
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'sm_wilt'
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'hcap   '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'hcon   '
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'albsoil'
  IF ( soil_bgc_model == soil_model_rothc .OR.                                &
       soil_bgc_model == soil_model_ecosse ) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'clay   '
  END IF
  IF ( soil_bgc_model == soil_model_ecosse ) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'soil_ph'
  END IF
END IF

!---------------------------------------------------------------------------
!Set up read/write details for topmodel ancils.
!Depends on the value of l_top

IF ( l_top ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'fexp   '
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

  nvars = nvars + 1
  identifiers(nvars) = 'ti_mean'
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF

  nvars = nvars + 1
  identifiers(nvars) = 'ti_sig '
  IF (l_append_soilt) THEN
    identifiers(nvars) = TRIM(identifiers(nvars)) // "_soilt"
  END IF


  IF ( .NOT. ancil_dump_read%top) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'fexp   '
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'ti_mean'
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'ti_sig '
  END IF
END IF !l_top

!---------------------------------------------------------------------------
!Set up read/write details for agric ancil namelist
!This is always needed

nvars = nvars + 1
identifiers(nvars) = 'frac_agr'

IF ( .NOT. ancil_dump_read%agric ) THEN
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'frac_agr'
END IF


!---------------------------------------------------------------------------
!Set up read/write details for soil properites ancil namelist
!The vars required depend on a couple of switches

IF ( l_crop ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'cropttveg'
  nvars = nvars + 1
  identifiers(nvars) = 'cropttrep'
  IF ( l_prescsow ) THEN  !we additionally need the crop sowing date
    nvars = nvars + 1
    identifiers(nvars) = 'cropsowdate'
    nvars = nvars + 1
    identifiers(nvars) = 'croplatestharvdate'
  END IF
    
  IF ( .NOT. ancil_dump_read%crop_props ) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'cropttveg'
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'cropttrep'
    IF ( l_prescsow ) THEN  !we additionally need the crop sowing date
      nvars_from_ancil = nvars_from_ancil + 1
      vars_from_ancil(nvars_from_ancil) = 'cropsowdate'
      nvars_from_ancil = nvars_from_ancil + 1
      vars_from_ancil(nvars_from_ancil) = 'croplatestharvdate'
    END IF
  END IF

END IF !l_crop

!---------------------------------------------------------------------------
!Set up read/write details for irrigation ancil namelist
!The vars required depend on l_irrig_dmd
IF ( l_irrig_dmd ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'frac_irrig'
  nvars = nvars + 1
  identifiers(nvars) = 'frac_irr_all_tiles'
  !These only matter if frac_irr_all_tiles is false, but can't gurantee that
  !that its value will be known when reading that variable from the dump
  nvars = nvars + 1
  identifiers(nvars) = 'irrtiles'
  nvars = nvars + 1
  identifiers(nvars) = 'nirrtile'

  IF ( .NOT. ancil_dump_read%irrig ) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'frac_irrig'
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'frac_irr_all_tiles'
    !These only matter if frac_irr_all_tiles is false, but can't gurantee that
    !that its value will be known when reading that variable from the dump
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'irrtiles'
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'nirrtile'
  END IF

END IF

!Rivers need further work before being included
!IF ( ancil_dump_read%rivers_props ) THEN
!END IF

  !---------------------------------------------------------------------------
  !Set up read/write details for co2 ancil namelist
  !Always required, even if it just carries around the default value set in
  !the aero module
nvars = nvars + 1
identifiers(nvars) = 'co2_mmr'

IF ( .NOT. ancil_dump_read%co2 ) THEN
  nvars_from_ancil = nvars_from_ancil + 1
  vars_from_ancil(nvars_from_ancil) = 'co2_mmr'
END IF

RETURN

END SUBROUTINE required_vars_for_configuration
#endif
