! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_vegetation_mod

USE max_dimensions, ONLY: npft_max, nsurft_max
USE missing_data_mod, ONLY: rmdi, imdi
USE crop_vars_mod, ONLY: nirrtile

!-----------------------------------------------------------------------------
! Description:
!   Contains vegetation options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
INTEGER, PARAMETER ::                                                         &
  i_veg_vn_1b = 1,                                                            &
      ! Constant indicating fixed vegetation scheme
  i_veg_vn_2b = 2
      ! Constant indicating interactive vegetation scheme

! Parameters identifying alternative models of leaf photosynthesis.
! These should have unique values.
! In all cases photosynthesis of C4 plants uses the model of Collatz et
! al., 1992, Aust.J.Plant Physiol., 19, 519-538.
INTEGER, PARAMETER ::                                                         &
  photo_collatz = 1,                                                          &
    ! C3 plants use the model of Collatz et al., 1991, Agricultural and
    ! Forest Meteorology, 54, 107-136.
  photo_farquhar = 2
    ! C3 plants use the model of Farquhar et al. ,1980, Planta, 149: 78-90.

! Parameters identifying alternative stomatal conductance models.
! These should have unique values.
INTEGER, PARAMETER ::                                                         &
  stomata_jacobs = 1,                                                         &
    ! Use the original JULES model, including the Jacobs closure
    !   - see Eqn.9 of Best et al. (2011), doi:10.5194/gmd-4-677-2011.
  stomata_medlyn = 2
    ! Use the model of Medlyn et al. (2011) - see Eqn.11,
    !   doi: 10.1111/j.1365-2486.2010.02375.x.

!-----------------------------------------------------------------------------
! Items set in namelist
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_nrun_mid_trif = .FALSE.,                                                  &
      ! Switch for starting NRUN mid-way through a TRIFFID period
  l_trif_init_accum = .FALSE.,                                                &
      ! Switch so that an NRUN will bit-compare with a CRUN when FALSE
  l_phenol = .FALSE.,                                                         &
      ! Switch for leaf phenology
  l_triffid = .FALSE.,                                                        &
      ! Switch for interactive veg model
  l_trif_eq = .FALSE.,                                                        &
      ! Switch for running TRIFFID in equilibrium mode
  l_veg_compete  = .FALSE.,                                                   &
      ! Switch for competing vegetation
      ! Setting l_triffid = .TRUE. and this as .FALSE. means that the carbon
      ! pools evolve but the PFT distribution does not change
  l_trait_phys  = .FALSE.,                                                    &
      ! .TRUE. for new trait-based PFTs (uses nmass & lma)
      ! .FALSE. for pre-new PFT configuration (nl0, sigl, and neff)
  l_ht_compete  = .FALSE.,                                                    &
      ! Switch for TRIFFID competition
      ! (T for height F for lotka)
      ! Must be true if npft > 5!
  l_landuse = .FALSE.,                                                        &
      ! Switch for landuse change that invokes wood product pools
  l_nitrogen = .FALSE.,                                                       &
      ! Switch for Nitrogen limiting NPP
  l_bvoc_emis = .FALSE.,                                                      &
      ! Switch to enable calculation of BVOC emissions
  l_o3_damage    = .FALSE.,                                                   &
      ! Switch for ozone damage
  l_prescsow = .FALSE.,                                                       &
      ! Only used if crop model is on ( ncpft > 0 )
      !   T => read in the sowing dates for each crop
      !   F => let the model determine sowing date
  l_irrig_dmd = .FALSE.,                                                      &
      !   Switch for using irrigation demand code
  l_irrig_limit = .FALSE.,                                                    &
      !   Switch for limiting irrigation supply
  l_recon = .TRUE.,                                                           &
      ! Used to switch on reconfiguration of veg fractions for TRIFFID
  l_trif_crop = .FALSE.,                                                      &
      ! switch to prevent crop and natural PFTs competing
  l_inferno = .FALSE.,                                                        &
      ! Switch used to control whether the Interactive fire scheme is used
  l_trif_fire = .FALSE.,                                                      &
      ! Switch used to control whether interactive fire is used
      !   T => if l_inferno is also true, g_burn is calculated in INFERNO
!   and passed to TRIFFID to calculate emissions and vegetation
       !   dynamics
       !   T => if l_inferno is false, interactive fire is calculated via
!   ancillary if provided, and is 0 if not provided
       !   F => g_burn is calculated via ancillary if provided, and is 0 if
!   not provided
   l_use_pft_psi = .FALSE.,                                                   &
       ! Switch used to control what parameters are used in the calculation
       ! of the soil moisture stress factor
       !   T => use psi_close and psi_open
       !   F => use sm_wilt, sm_crit and fsmc_p0
   l_spec_veg_z0 = .FALSE.,                                                   &
       ! Switch used to specify the roughness length for vegetation rather
       ! than calculate it from the vegetation canopy height
   l_limit_canhc = .FALSE.,                                                   &
       ! Switch to limit the value of vegetation canopy heat capacity to
       ! the value specified in the subroutine CANCAP

! Switches for bug fixes.
    l_leaf_n_resp_fix = .FALSE.,                                              &
        ! Switch to use correct forms for canopy-average leaf nitrogen.
        ! This affects can_rad_mod = 1, 4 and 5, not 6 (which is correct).
    l_stem_resp_fix = .FALSE.,                                                &
        ! Switch used to control whether LAI or LAI_BAL is used in stem
        ! resp calculation. Only affects non-crop PFTs with l_trait_phys=F.
    l_vegcan_soilfx = .FALSE.,                                                &
        ! Switch to modify the canopy model to allow for conduction through
        ! the soil below vegetation.
    l_gleaf_fix = .FALSE.,                                                    &
        ! Used to fix a bug accumulating g_leaf_phen_ac in standalone JULES
        !
        ! This bug occurs in standalone JULES because veg2 is called on TRIFFID
        ! timesteps and veg1 is called on phenol timesteps
        !
        ! In the UM, veg2 (where the accumulation is working) is called on
        ! TRIFFID AND phenol timesteps if l_triffid = TRUE and veg1 is called
        ! on phenol timesteps if l_triffid = FALSE, hence this bug doesn't arise
        !
        ! This means we don't bother adding it to the UM namelist transfer
        ! between PEs
    l_scale_resp_pm = .FALSE.,                                                &
        ! Switch for scaling whole plant maintenance respiraiton by the soil
        ! moisture stress factor. FALSE = Only scale leaf respiration;
        ! TRUE = scale whole plant respiration.
    l_vegdrag_surft(nsurft_max),                                              &
        ! Switch for using vegetation canopy drag scheme on each tile.
        ! Must be false for non-PFT tiles
    l_vegdrag_pft(npft_max),                                                  &
        ! Switch for using vegetation canopy drag scheme.
        ! This is what appears in the namelist, to ensure that true
        ! values can only given for pfts
    l_rsl_scalar = .FALSE.
        ! Switch for using roughness sublayer correction scheme in scalar
        ! variables. This is only valid when l_vegdrag_surft = .true.

DATA l_vegdrag_surft / nsurft_max * .FALSE. /
DATA l_vegdrag_pft / npft_max * .FALSE. /

INTEGER ::                                                                    &
  can_model = 4,                                                              &
      ! Switch for thermal vegetation
  can_rad_mod = 4,                                                            &
      ! Canopy radiation model
  ilayers = imdi,                                                             &
      ! Number of layers for canopy radiation model
  irr_crop = 0,                                                               &
      ! Switch for irrigation cropping model, default is 0.
  irrigtiles(npft_max)
      ! Indices of the tiles to be irrigated. Only used if frac_irrig_all_tiles = F.
LOGICAL ::                                                                    &
  frac_irrig_all_tiles = .TRUE.
      ! T - assign irr fraction to all tiles
      ! F - assign irr fraction to specfied tiles only

INTEGER ::                                                                    &
  ignition_method = 1,                                                        &
      ! Switch for the calculation method of INFERNO fire ignitions
      ! IGNITION_METHOD=1:Constant (1.67 per km2 per s)
      ! IGNITION_METHOD=2:Constant (Human - 1.5 per km2 per s)
      !                   Varying  (Lightning - see Pechony and Shindell,2009)
      ! IGNITION_METHOD=3:Vary Human and Lightning (Pechony and Shindell,2009)
  fsmc_shape = 0,                                                             &
      ! shape of the soil moisture stress function fsmc
      ! 0: piece-wise linear in vol. soil moisture.
      ! 1: piece-wise linear in soil potential.
  photo_model = imdi,                                                         &
      ! Chosen model of leaf photosynthesis.
  stomata_model = stomata_jacobs
      ! Stomatal conductance model.

INTEGER ::                                                                    &
  phenol_period = imdi,                                                       &
      ! Update frequency for leaf phenology (days)
  triffid_period = imdi
      ! Update frequency for TRIFFID (days)

INTEGER :: errcode   ! error code to pass to ereport.

REAL(KIND=real_jlslsm) ::                                                     &
  frac_min  = rmdi,                                                           &
      ! Minimum areal fraction for PFTs.
  frac_seed = rmdi,                                                           &
      ! "Seed" fraction for PFTs.
  pow = rmdi,                                                                 &
      ! Power in sigmoidal function.
  cd_leaf = rmdi,                                                             &
      ! Leaf level drag coefficient
  c1_usuh = rmdi,                                                             &
      ! u*/U(h) at the top of dense canopy
  c2_usuh = rmdi,                                                             &
      ! u*/U(h) above bare soil
  c3_usuh = rmdi,                                                             &
      ! Used in the exponent of equation weighting dense and sparse
      ! vegetation to get u*/U(h) in neutral condition
  stanton_leaf = rmdi
      ! Leaf-level Stanton number

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_vegetation/                                                 &
! UM only
    l_nrun_mid_trif, l_trif_init_accum,                                       &
! Shared
    l_phenol, l_triffid, l_trif_eq, l_veg_compete,                            &
    phenol_period, triffid_period, l_trait_phys, l_ht_compete,                &
    l_bvoc_emis, l_o3_damage, can_model, can_rad_mod, ilayers,                &
    frac_min, frac_seed, pow, l_landuse, l_leaf_n_resp_fix, l_stem_resp_fix,  &
    l_nitrogen, l_vegcan_soilfx, l_trif_crop, l_trif_fire,                    &
    l_inferno, ignition_method, l_vegdrag_pft, l_rsl_scalar, l_irrig_dmd,     &
    cd_leaf, c1_usuh, c2_usuh, c3_usuh, stanton_leaf, photo_model,            &
    stomata_model, l_spec_veg_z0, l_limit_canhc,                              &
! Not used in the UM yet
    l_prescsow, l_recon, l_irrig_limit, irr_crop,l_gleaf_fix,                 &
    l_scale_resp_pm, l_use_pft_psi, fsmc_shape

!-----------------------------------------------------------------------------
! Items derived from namelist inputs
!-----------------------------------------------------------------------------
INTEGER :: i_veg_vn = imdi
    ! Switch to determine version of vegetation scheme
    ! Must be one of i_veg_vn_1b or i_veg_vn_2b

LOGICAL :: l_crop = .FALSE.
    ! Indicates if the crop model is on
    ! This is derived from the value of ncpft

LOGICAL :: l_fapar_diag = .FALSE.
    ! Only true if fapar or apar is in one of the output profiles

LOGICAL :: l_fao_ref_evapotranspiration = .FALSE.
    ! Only true if fao_et0 is in one of the output profiles

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_VEGETATION_MOD'

CONTAINS


SUBROUTINE check_jules_vegetation()

USE ereport_mod, ONLY: ereport

USE jules_surface_types_mod, ONLY: npft, ncpft, nnpft, c3_grass, c4_grass

USE jules_surface_mod, ONLY: l_aggregate

USE jules_hydrology_mod, ONLY: l_top

USE jules_print_mgr, ONLY: jules_message

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_VEGETATION namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE


! Phenology or TRIFFID cannot be used with the aggregate surface scheme
IF ( l_aggregate .AND. (l_phenol .OR. l_triffid) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Phenology or TRIFFID cannot be used with the ' //             &
               'aggregated surface scheme (i.e. l_aggregate = true)')
END IF

! Check that phenol_period is specified if phenology is on
IF ( l_phenol .AND. phenol_period < 0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Phenology is on but phenol_period is not given')
END IF
! Same for triffid_period if TRIFFID is on
IF ( l_triffid .AND. triffid_period < 0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'TRIFFID is on but triffid_period is not given')
END IF

! Check options that depend on TRIFFID if it is not enabled
IF ( .NOT. l_triffid .AND. ANY( (/ l_veg_compete, l_trif_eq, l_landuse,       &
   l_ht_compete, l_nitrogen, l_trif_crop, l_trif_fire /) ) ) THEN
  errcode = 101
  WRITE(jules_message,'(A,8(1x,L1))')                                         &
     'These should be false when l_triffid = F: l_veg_compete, ' //           &
     'l_trif_eq, l_landuse, l_ht_compete, l_nitrogen, l_trif_crop, ' //       &
     'l_trif_fire = ', l_veg_compete, l_trif_eq, l_landuse, l_ht_compete,     &
     l_nitrogen, l_trif_crop, l_trif_fire
  CALL ereport( "check_jules_vegetation", errcode, jules_message )
END IF

! Always make sure that a veg version is selected
! If TRIFFID is on, select interactive veg, otherwise select fixed veg
i_veg_vn = i_veg_vn_1b
IF ( l_triffid ) THEN
  i_veg_vn = i_veg_vn_2b
END IF

! Check can_model and can_rad_mod are suitable
IF ( can_model < 1 .OR. can_model > 4 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'can_model should be in range 1 to 4')
END IF

IF ( can_model == 4 .AND. l_aggregate ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'can_model=4 cannot be used with the aggregated ' //           &
               'surface scheme')
END IF

SELECT CASE ( can_rad_mod )
CASE ( 1, 4, 5, 6 )
  ! These are valid, so nothing to do.
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'can_rad_mod should be 1, 4, 5 or 6')
END SELECT

! Check that the photosynthesis option is reasonable.
SELECT CASE ( photo_model )
CASE ( photo_collatz, photo_farquhar )
  ! These are valid, nothing more to do.
CASE DEFAULT
  errcode = 101  !  a fatal error
  CALL ereport("check_jules_vegetation", errcode,                             &
               "Invalid value given for photo_model.")
END SELECT

! The Farquhar model of photosynthesis has only been coded for certain
! values of can_rad_mod.
IF ( photo_model == photo_farquhar ) THEN
  SELECT CASE ( can_rad_mod )
  CASE ( 1, 5, 6 )
    ! These are supported, nothing more to do.
  CASE DEFAULT
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                           &
                 "Farquhar model can only be used with can_rad_mod =1, 5 " // &
                 " or 6.")
  END SELECT
END IF

SELECT CASE ( stomata_model )
CASE ( stomata_jacobs, stomata_medlyn )
  ! These are valid, so nothing to do.
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               "Invalid value for stomata_model" )
END SELECT

IF ( l_triffid .AND. ( .NOT. l_phenol ) ) THEN
  errcode = -105 ! warning
  CALL ereport("check_jules_vegetation", errcode,                             &
               "When triffid is on, l_phenol=T is recommended. You " //       &
               "have set l_phenol=F. The LAI will be set to the " //          &
               "balanced LAI.")
END IF

! Check triffid-crop options are sensible
IF ( l_trif_crop .AND. l_trif_eq ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'trif_crop and trif_eq are incompatible')
END IF

IF (l_veg_compete .AND. ( .NOT. l_ht_compete ) .AND. ( nnpft /= 5 )) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'l_ht_compete=F requires 5 natural PFTs: ' //                  &
               'BT, NT, C3, C4, SH')
END IF

! Check crop options are sensible
l_crop = ncpft > 0

IF ( l_crop .AND. l_triffid ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Crop model and triffid are incompatible')
END IF

IF ( l_aggregate .AND. l_crop ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Crop model cannot be used with the aggregated surface ' //    &
               'scheme (i.e. l_aggregate = true)')
END IF

! Irrigation demand is in the process of being implemented in the UM in a
! limited configuration. As a starting point it will only be activated all year
! round (i.e. irr_crop = 0) and only used on short vegetation (i.e. C3 and C4
! grass). There is still some work to be done even in this limited state and
! until this is completed it will be triggered off and will return a fatal
! ereport if turned on to prevent its accidental use. Work will progress in
! jules:#949.
#if defined(UM_JULES)
IF ( l_irrig_dmd ) THEN
  irr_crop = 0
  frac_irrig_all_tiles = .FALSE.
  nirrtile = 2
  irrigtiles(1) = c3_grass
  irrigtiles(2) = c4_grass
  errcode = -1
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Irrigation demand is in the process of being implemented' //  &
               ' in the UM in a limited way. It is currently turned off.')
END IF
#endif

IF ( l_irrig_limit .AND. .NOT. l_irrig_dmd ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'l_irrig_limit=T requires l_irrig_dmd=T ')
END IF

IF ( l_irrig_limit .AND. .NOT. l_top ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'l_irrig_limit=T requires l_top=T ')
END IF

IF ( l_irrig_dmd .AND. irr_crop == 2 .AND. .NOT. l_crop ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Irrigation triggered by crop DVI (irr_crop = 2) ' //          &
               'requires crop model to be active')
END IF

IF ( .NOT. l_crop ) THEN
  ! Note that in the UM, since ncpft_max = 0, l_crop is always .FALSE. and hence
  ! l_prescsow is .FALSE.
  l_prescsow = .FALSE.
END IF

! Check a suitable ignition_method was given
IF ( ignition_method < 1 .OR. ignition_method > 3 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'ignition_method must be 1, 2 or 3')
END IF

IF ( fsmc_shape == 1 .AND. .NOT. l_use_pft_psi ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'fsmc_shape=1 requires l_use_pft_psi=T ')
END IF

IF ( l_aggregate .AND. ANY(l_vegdrag_pft) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                             &
               'Vegetative drag scheme cannot be used with the ' //           &
               'aggregated surface scheme (i.e. l_aggregate = true)')
ELSE
  ! Copy the values for the given number of pfts across
  l_vegdrag_surft(1:npft) = l_vegdrag_pft(1:npft)
END IF

END SUBROUTINE check_jules_vegetation

SUBROUTINE print_nlist_jules_vegetation()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('jules_vegetation_mod',                                      &
                 'Contents of namelist jules_vegetation')

WRITE(lineBuffer,*)' l_nrun_mid_trif = ', l_nrun_mid_trif
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_init_accum = ', l_trif_init_accum
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_phenol = ',l_phenol
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_triffid = ',l_triffid
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_eq = ',l_trif_eq
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_veg_compete = ',l_veg_compete
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_ht_compete = ',l_ht_compete
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_crop = ',l_trif_crop
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_fire = ',l_trif_fire
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trait_phys = ',l_trait_phys
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_landuse = ',l_landuse
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_leaf_n_resp_fix = ',l_leaf_n_resp_fix
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_stem_resp_fix = ',l_stem_resp_fix
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_scale_resp_pm = ',l_scale_resp_pm
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_nitrogen = ',l_nitrogen
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_vegcan_soilfx = ',l_vegcan_soilfx
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' phenol_period = ',phenol_period
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' triffid_period = ',triffid_period
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_bvoc_emis = ',l_bvoc_emis
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_o3_damage = ',l_o3_damage
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_vegdrag_pft = ',l_vegdrag_pft
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_rsl_scalar = ',l_rsl_scalar
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_irrig_dmd = ',l_irrig_dmd
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' can_model = ',can_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' can_rad_mod = ',can_rad_mod
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' ilayers = ',ilayers
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' frac_min = ',frac_min
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' frac_seed = ',frac_seed
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' pow = ',pow
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' cd_leaf = ',cd_leaf
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c1_usuh = ',c1_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c2_usuh = ',c2_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c3_usuh = ',c3_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' stanton_leaf = ',stanton_leaf
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' photo_model = ',photo_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' stomata_model = ',stomata_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

END SUBROUTINE print_nlist_jules_vegetation

#if defined(UM_JULES) && !defined(LFRIC)

SUBROUTINE read_nml_jules_vegetation (unitnumber)

! Description:
!  Read the JULES_VEGETATION namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype


USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_VEGETATION'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 8
INTEGER, PARAMETER :: n_real = 8
INTEGER, PARAMETER :: n_log = 25 + npft_max

TYPE my_namelist
  SEQUENCE
  INTEGER :: phenol_period
  INTEGER :: triffid_period
  INTEGER :: can_model
  INTEGER :: can_rad_mod
  INTEGER :: ilayers
  INTEGER :: ignition_method
  INTEGER :: photo_model
  INTEGER :: stomata_model
  REAL(KIND=real_jlslsm) :: frac_min
  REAL(KIND=real_jlslsm) :: frac_seed
  REAL(KIND=real_jlslsm) :: pow
  REAL(KIND=real_jlslsm) :: cd_leaf
  REAL(KIND=real_jlslsm) :: c1_usuh
  REAL(KIND=real_jlslsm) :: c2_usuh
  REAL(KIND=real_jlslsm) :: c3_usuh
  REAL(KIND=real_jlslsm) :: stanton_leaf
  LOGICAL :: l_nrun_mid_trif
  LOGICAL :: l_trif_init_accum
  LOGICAL :: l_phenol
  LOGICAL :: l_triffid
  LOGICAL :: l_trif_eq
  LOGICAL :: l_veg_compete
  LOGICAL :: l_bvoc_emis
  LOGICAL :: l_o3_damage
  LOGICAL :: l_prescsow
  LOGICAL :: l_trait_phys
  LOGICAL :: l_ht_compete
  LOGICAL :: l_trif_crop
  LOGICAL :: l_trif_fire
  LOGICAL :: l_landuse
  LOGICAL :: l_nitrogen
  LOGICAL :: l_recon
  LOGICAL :: l_leaf_n_resp_fix
  LOGICAL :: l_stem_resp_fix
  LOGICAL :: l_scale_resp_pm
  LOGICAL :: l_vegcan_soilfx
  LOGICAL :: l_inferno
  LOGICAL :: l_vegdrag_pft(npft_max)
  LOGICAL :: l_rsl_scalar
  LOGICAL :: l_spec_veg_z0
  LOGICAL :: l_limit_canhc
  LOGICAL :: l_irrig_dmd
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,              &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_vegetation, IOSTAT = errorstatus,      &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_vegetation", iomessage)

  my_nml % phenol_period   = phenol_period
  my_nml % triffid_period  = triffid_period
  my_nml % can_model       = can_model
  my_nml % can_rad_mod     = can_rad_mod
  my_nml % ilayers         = ilayers
  my_nml % ignition_method = ignition_method
  my_nml % photo_model     = photo_model
  my_nml % stomata_model   = stomata_model
  my_nml % frac_min        = frac_min
  my_nml % frac_seed       = frac_seed
  my_nml % pow             = pow
  my_nml % cd_leaf         = cd_leaf
  my_nml % c1_usuh         = c1_usuh
  my_nml % c2_usuh         = c2_usuh
  my_nml % c3_usuh         = c3_usuh
  my_nml % stanton_leaf    = stanton_leaf
  my_nml % l_nrun_mid_trif = l_nrun_mid_trif
  my_nml % l_trif_init_accum   = l_trif_init_accum
  my_nml % l_phenol        = l_phenol
  my_nml % l_triffid       = l_triffid
  my_nml % l_trif_eq       = l_trif_eq
  my_nml % l_veg_compete   = l_veg_compete
  my_nml % l_bvoc_emis     = l_bvoc_emis
  my_nml % l_o3_damage     = l_o3_damage
  my_nml % l_prescsow      = l_prescsow
  my_nml % l_trait_phys    = l_trait_phys
  my_nml % l_ht_compete    = l_ht_compete
  my_nml % l_trif_crop     = l_trif_crop
  my_nml % l_trif_fire     = l_trif_fire
  my_nml % l_landuse       = l_landuse
  my_nml % l_nitrogen      = l_nitrogen
  my_nml % l_recon         = l_recon
  my_nml % l_leaf_n_resp_fix = l_leaf_n_resp_fix
  my_nml % l_stem_resp_fix = l_stem_resp_fix
  my_nml % l_scale_resp_pm = l_scale_resp_pm
  my_nml % l_vegcan_soilfx = l_vegcan_soilfx
  my_nml % l_inferno       = l_inferno
  my_nml % l_vegdrag_pft   = l_vegdrag_pft
  my_nml % l_rsl_scalar    = l_rsl_scalar
  my_nml % l_spec_veg_z0   = l_spec_veg_z0
  my_nml % l_limit_canhc   = l_limit_canhc
  my_nml % l_irrig_dmd     = l_irrig_dmd
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  phenol_period   = my_nml % phenol_period
  triffid_period  = my_nml % triffid_period
  can_model       = my_nml % can_model
  can_rad_mod     = my_nml % can_rad_mod
  ilayers         = my_nml % ilayers
  ignition_method = my_nml % ignition_method
  photo_model     = my_nml % photo_model
  stomata_model   = my_nml % stomata_model
  frac_min        = my_nml % frac_min
  frac_seed       = my_nml % frac_seed
  pow             = my_nml % pow
  cd_leaf         = my_nml % cd_leaf
  c1_usuh         = my_nml % c1_usuh
  c2_usuh         = my_nml % c2_usuh
  c3_usuh         = my_nml % c3_usuh
  stanton_leaf    = my_nml % stanton_leaf
  l_nrun_mid_trif = my_nml % l_nrun_mid_trif
  l_trif_init_accum = my_nml % l_trif_init_accum
  l_phenol        = my_nml % l_phenol
  l_triffid       = my_nml % l_triffid
  l_trif_eq       = my_nml % l_trif_eq
  l_veg_compete   = my_nml % l_veg_compete
  l_bvoc_emis     = my_nml % l_bvoc_emis
  l_o3_damage     = my_nml % l_o3_damage
  l_prescsow      = my_nml % l_prescsow
  l_trait_phys    = my_nml % l_trait_phys
  l_ht_compete    = my_nml % l_ht_compete
  l_trif_crop     = my_nml % l_trif_crop
  l_trif_fire     = my_nml % l_trif_fire
  l_landuse       = my_nml % l_landuse
  l_nitrogen      = my_nml % l_nitrogen
  l_recon         = my_nml % l_recon
  l_leaf_n_resp_fix = my_nml % l_leaf_n_resp_fix
  l_stem_resp_fix = my_nml % l_stem_resp_fix
  l_scale_resp_pm = my_nml % l_scale_resp_pm
  l_vegcan_soilfx = my_nml % l_vegcan_soilfx
  l_inferno       = my_nml % l_inferno
  l_vegdrag_pft   = my_nml % l_vegdrag_pft
  l_rsl_scalar    = my_nml % l_rsl_scalar
  l_spec_veg_z0   = my_nml % l_spec_veg_z0
  l_limit_canhc   = my_nml % l_limit_canhc
  l_irrig_dmd     = my_nml % l_irrig_dmd
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_vegetation
#endif

END MODULE jules_vegetation_mod
