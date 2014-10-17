#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_model_levels(nml_dir)

  USE io_constants, ONLY : NAMELIST_UNIT

  USE max_dimensions, ONLY : nnvg_max, npft_max, ntiles_max, ntype_max,       &
                             sm_levels_max, snow_layers_max

  USE string_utils_mod, ONLY : to_string

  USE switches, ONLY : l_aggregate, l_triffid, l_trif_eq

  USE nstypes, ONLY : npft,nnvg,urban,lake,soil,ice,urban_canyon,urban_roof
  USE ancil_info, ONLY : sm_levels,nsmax

  USE nstypes, ONLY : ntype
  USE ancil_info, ONLY : ntiles, npft_trif, dim_cs1, nice, nice_use

  USE model_interface_mod, ONLY : pft_dim_size, nvg_dim_size, type_dim_size,  &
                                  tile_dim_size, snow_dim_size,               &
                                  soil_dim_size, scpool_dim_size

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the vertical levels and checks them for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of jules_model_levels namelist - this combines variables in
! different modules into a logically coherent namelist
!-----------------------------------------------------------------------------
  NAMELIST /jules_model_levels/ npft,nnvg,urban,lake,soil,ice,urban_canyon,   &
                                urban_roof,sm_levels,nsmax


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! First, we read the model levels namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_model_levels", "Reading JULES_MODEL_LEVELS namelist...")

  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'model_levels.nml'),      &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_model_levels",                                       &
                   "Error opening namelist file model_levels.nml " //         &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_model_levels, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_model_levels",                                       &
                   "Error reading namelist JULES_MODEL_LEVELS " //            &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_model_levels",                                       &
                   "Error closing namelist file model_levels.nml " //         &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Set derived variables and check for consistency
!-----------------------------------------------------------------------------
! Calculate number of surface types.
  ntype = npft + nnvg
  IF ( l_aggregate ) THEN
    ntiles = 1
  ELSE
    ntiles = ntype
  ENDIF

! TRIFFID needs 5 pfts
  IF ( ( l_triffid .OR. l_trif_eq ) .AND. npft /= 5 )                         &
    CALL log_fatal("init_model_levels",                                       &
                   'TRIFFID is hardwired to expect certain PFTs - ' //        &
                   'npft must be 5, and they must be the "expected" PFTs')

! Only one ice category
  nice = 1
! nice_use MUST EQUAL NICE FOR STANDALONE JULES
  nice_use = nice

! Set dimensions that are dependent on TRIFFID
  IF ( l_triffid ) THEN
    npft_trif = npft
!   Soil carbon dimensions
!   Set dim_cs1=4 for the 4 pools of RothC.
    dim_cs1 = 4
  ELSE
    npft_trif = 1
!   Soil carbon dimensions
!   Set dim_cs1=1 to use a single soil C pool.
    dim_cs1 = 1
  END IF
    npft_trif = npft

! Set the dimension sizes for IO
  pft_dim_size    = npft
  nvg_dim_size    = nnvg
  type_dim_size   = ntype
  tile_dim_size   = ntiles
  snow_dim_size   = nsmax
  soil_dim_size   = sm_levels
  scpool_dim_size = dim_cs1

!-----------------------------------------------------------------------------
! Check that fixed-length i/o dimensions are large enough.
!-----------------------------------------------------------------------------
  IF ( nnvg > nnvg_max )                                                      &
    CALL log_fatal("init_model_levels",                                       &
                   "nnvg > nnvg_max - increase nnvg_max and recompile")

  IF ( npft > npft_max )                                                      &
    CALL log_fatal("init_model_levels",                                       &
                   "npft > npft_max - increase npft_max and recompile")

  IF ( ntiles > ntiles_max )                                                  &
    CALL log_fatal("init_model_levels",                                       &
                   "ntiles > ntiles_max - increase ntiles_max and recompile")

  IF ( ntype > ntype_max )                                                    &
    CALL log_fatal("init_model_levels",                                       &
                   "ntype > ntype_max - increase ntype_max and recompile")

  IF ( sm_levels > sm_levels_max )                                            &
    CALL log_fatal("init_model_levels",                                       &
                   "sm_levels > sm_levels_max - increase sm_levels_max " //   &
                   "and recompile")

  IF ( nsmax > snow_layers_max )                                              &
    CALL log_fatal("init_model_levels",                                       &
                   "nsmax > snow_layers_max - increase snow_layers_max " //   &
                   "and recompile")

!-----------------------------------------------------------------------------
! Check values for the specific surface types are sensible
!-----------------------------------------------------------------------------
! Check that we have soil type, and also check if lake and ice types available.
! We need the soil type since this is used for the bare soil under PFTs.

! Note that the non-veg surface types must come after the PFTs since
! (a) PHYSIOL refers to type soil-npft
! (b) INIT_NONVEG refers to z0h_z0m(npft+i)

! Check that all the given specific surface types are within range
  IF ( urban > ntype )                                                        &
    CALL log_fatal("init_model_levels",                                       &
                   "urban tile is out of range (> ntype)")

  IF ( lake > ntype )                                                         &
    CALL log_fatal("init_model_levels",                                       &
                   "lake tile is out of range (> ntype)")

  IF ( soil > ntype )                                                         &
    CALL log_fatal("init_model_levels",                                       &
                   "soil tile is out of range (> ntype)")

  IF ( ice > ntype )                                                          &
    CALL log_fatal("init_model_levels",                                       &
                   "ice tile is out of range (> ntype)")

  IF ( urban_canyon > ntype )                                                 &
    CALL log_fatal("init_model_levels",                                       &
                   "urban_canyon tile is out of range (> ntype)")

  IF ( urban_roof > ntype )                                                   &
    CALL log_fatal("init_model_levels",                                       &
                   "urban_roof tile is out of range (> ntype)")

  IF ( soil < 1 )                                                             &
    CALL log_fatal("init_model_levels", "No soil surface type specified")


!-----------------------------------------------------------------------------
! Log some info about the model levels
!-----------------------------------------------------------------------------
  CALL log_info("init_model_levels",                                          &
                "Using values - " //                                          &
                "npft = " // TRIM(to_string(npft)) // "; " //                 &
                "nnvg = " // TRIM(to_string(nnvg)) // "; " //                 &
                "ntiles = " // TRIM(to_string(ntiles)) // "; " //             &
                "sm_levels = " // TRIM(to_string(sm_levels)))
  IF ( nsmax > 0 ) THEN
    CALL log_info("init_model_levels",                                        &
                  "Using 'new' snow scheme with nsmax = " //                  &
                  TRIM(to_string(nsmax)))
  ELSE
    CALL log_info("init_model_levels", "Using 'old' snow scheme")
  END IF


  CALL log_info("init_model_levels",                                          &
                "Soil is surface type #" // TRIM(to_string(soil)))

  IF ( lake > 0 ) THEN
    CALL log_info("init_model_levels",                                          &
                  "Lake (inland water) is surface type #" // TRIM(to_string(lake)))
  ELSE
    CALL log_info("init_model_levels",                                          &
                  "No lake (inland water) type has been indicated")
  END IF

  IF ( ice > 0 ) THEN
    CALL log_info("init_model_levels",                                          &
                  "Land ice is surface type #" // TRIM(to_string(ice)))
  ELSE
    CALL log_info("init_model_levels",                                          &
                  "No land ice type has been indicated")
  END IF

  IF ( urban > 0 ) THEN
    CALL log_info("init_model_levels",                                          &
                  "Urban is surface type #" // TRIM(to_string(urban)))
  ELSE
    CALL log_info("init_model_levels",                                          &
                  "No urban type has been indicated (URBAN-1T)")
  END IF

  IF ( urban_canyon > 0 ) THEN
    CALL log_info("init_model_levels",                                          &
                  "Urban canyon is surface type #" // TRIM(to_string(urban_canyon)))
  ELSE
    CALL log_info("init_model_levels",                                          &
                  "No canyon type has been indicated (URBAN-2T or MORUSES)")
  END IF

  IF ( urban_roof > 0 ) THEN
    CALL log_info("init_model_levels",                                          &
                  "Urban roof is surface type #" // TRIM(to_string(urban_roof)))
  ELSE
    CALL log_info("init_model_levels",                                          &
                  "No roof type has been indicated (URBAN-2T or MORUSES)")
  END IF

  RETURN

END SUBROUTINE init_model_levels
#endif
