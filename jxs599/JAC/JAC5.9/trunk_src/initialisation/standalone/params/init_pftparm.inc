#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_pftparm_jules(nml_dir,progs)

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_soil_biogeochem_mod, ONLY: l_layeredC, soil_bgc_model,              &
                                      soil_model_rothc

USE jules_surface_types_mod, ONLY: npft, nnpft

USE ancil_info, ONLY: land_pts

USE pftparm, ONLY: act_jmax, act_vcmax, deact_jmax,                           &
deact_vcmax, albsnc_max, albsnc_min, albsnf_maxu,                             &
albsnf_max, albsnf_maxl, alpha_elec, alpha, alniru,                           &
alnir, alnirl, alparu, alpar, alparl, a_wl, a_ws, b_wl,                       &
catch0, c3, dcatch_dlai, dgl_dm, dgl_dt, dqcrit, ds_jmax,                     &
ds_vcmax, dz0v_dh, z0v, emis_pft, eta_sl, fd, fsmc_of, f0,                    &
g1_stomata, glmin, g_leaf_0, infil_f, jv25_ratio, kext,                       &
kpar, lai_alb_lim, neff, nl0, nr_nl, ns_nl, nsw, nr, hw_sw,                   &
can_struct_a, gsoil_f, omegau, omega, omegal, omniru, omnir, omnirl,          &
orient, r_grow, rootd_ft, psi_close, psi_open, fsmc_p0, fsmc_mod,             &
sigl, tleaf_of, tlow, tupp, lma, nmass, vsl, vint, kn, knl, q10_leaf,         &
fl_o3_ct, dfp_dcuo, ci_st, gpp_st, ief, tef, mef, aef, ccleaf_min,            &
ccleaf_max, ccwood_min, ccwood_max, avg_ba, fire_mort, fef_co2,               &
fef_co, fef_ch4, fef_nox, fef_so2, fef_oc, fef_bc, dust_veg_scj

USE c_z0h_z0m, ONLY: z0h_z0m, z0h_z0m_classic

USE pftparm_io, ONLY: jules_pftparm,                                          &
                         fsmc_mod_io, canht_ft_io, lai_io,                    &
                         psi_close_io, psi_open_io,                           &
                         c3_io,orient_io,a_wl_io,a_ws_io,                     &
                         act_jmax_io,act_vcmax_io,                            &
                         albsnc_max_io,albsnc_min_io,albsnf_maxu_io,          &
                         albsnf_max_io,albsnf_maxl_io,                        &
                         alpha_io,alpha_elec_io,alniru_io,alnir_io,           &
                         alnirl_io,alparu_io,alpar_io,alparl_io,b_wl_io,      &
                         catch0_io,dcatch_dlai_io,deact_jmax_io,              &
                         deact_vcmax_io,dgl_dm_io,                            &
                         dgl_dt_io,dqcrit_io,ds_jmax_io,                      &
                         ds_vcmax_io,dz0v_dh_io,z0v_io,eta_sl_io,             &
                         fd_io,fsmc_of_io,fsmc_p0_io,f0_io,g1_stomata_io,     &
                         g_leaf_0_io, glmin_io,infil_f_io,jv25_ratio_io,      &
                         kext_io,kpar_io,lai_alb_lim_io,                      &
                         neff_io,nl0_io,nr_nl_io,ns_nl_io,                    &
                         nr_io,nsw_io,hw_sw_io,                               &
                         can_struct_a_io,gsoil_f_io,                          &
                         omegau_io,omega_io,omegal_io,omniru_io,              &
                         omnir_io,omnirl_io,r_grow_io,rootd_ft_io,            &
                         sigl_io,tleaf_of_io,tlow_io,tupp_io,                 &
                         emis_pft_io,z0hm_pft_io,z0hm_classic_pft_io,         &
                         dust_veg_scj_io,fl_o3_ct_io,dfp_dcuo_io,             &
                         ci_st_io,gpp_st_io,                                  &
                         ief_io,tef_io,mef_io,aef_io,                         &
                         fef_co2_io, fef_co_io, fef_ch4_io,                   &
                         fef_nox_io, fef_so2_io, fef_oc_io,                   &
                         fef_bc_io,ccleaf_min_io, ccleaf_max_io,              &
                         ccwood_min_io, ccwood_max_io, avg_ba_io,             &
                         lma_io,nmass_io,vsl_io,vint_io,kn_io,knl_io,         &
                         q10_leaf_io, fire_mort_io


USE jules_vegetation_mod, ONLY: can_rad_mod, l_crop, l_trait_phys,            &
                                 l_use_pft_psi, l_bvoc_emis, l_inferno,       &
                                 l_o3_damage, l_trif_fire, photo_acclim_model,&
                                 photo_farquhar, photo_model,                 &
                                 stomata_jacobs, stomata_medlyn,              &
                                 stomata_model, l_spec_veg_z0

USE jules_radiation_mod, ONLY: l_spec_albedo, l_albedo_obs, l_snow_albedo

USE errormessagelength_mod, ONLY: errormessagelength

!TYPE definitions
USE prognostics, ONLY: progs_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the PFT parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

!TYPES containing the data
TYPE(progs_type), INTENT(IN OUT) :: progs

! Work variables
INTEGER :: error  ! Error indicator
INTEGER :: i      ! Loop counter.
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=*), PARAMETER :: routinename='INIT_PFTPARM_JULES'

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
c3_io(:)                    = -9
orient_io(:)                = -9
fsmc_mod_io(:)              = -9

a_wl_io(:)                  = rmdi
a_ws_io(:)                  = rmdi
act_jmax_io(:)              = rmdi
act_vcmax_io(:)             = rmdi
albsnc_max_io(:)            = rmdi
albsnc_min_io(:)            = rmdi
albsnf_maxu_io(:)           = rmdi
albsnf_max_io(:)            = rmdi
albsnf_maxl_io(:)           = rmdi
alpha_elec_io(:)            = rmdi
alpha_io(:)                 = rmdi
alniru_io(:)                = rmdi
alnir_io(:)                 = rmdi
alnirl_io(:)                = rmdi
alparu_io(:)                = rmdi
alpar_io(:)                 = rmdi
alparl_io(:)                = rmdi
b_wl_io(:)                  = rmdi
catch0_io(:)                = rmdi
dcatch_dlai_io(:)           = rmdi
deact_jmax_io(:)            = rmdi
deact_vcmax_io(:)           = rmdi
dgl_dm_io(:)                = rmdi
dgl_dt_io(:)                = rmdi
dqcrit_io(:)                = rmdi
ds_jmax_io(:)               = rmdi
ds_vcmax_io(:)              = rmdi
dz0v_dh_io(:)               = rmdi
z0v_io(:)                   = rmdi
eta_sl_io(:)                = rmdi
fd_io(:)                    = rmdi
fsmc_of_io(:)               = rmdi
f0_io(:)                    = rmdi
g1_stomata_io(:)            = rmdi
g_leaf_0_io(:)              = rmdi
glmin_io(:)                 = rmdi
infil_f_io(:)               = rmdi
jv25_ratio_io(:)            = rmdi
kext_io(:)                  = rmdi
kpar_io(:)                  = rmdi
lai_alb_lim_io(:)           = rmdi
neff_io(:)                  = rmdi
nl0_io(:)                   = rmdi
nr_nl_io(:)                 = rmdi
ns_nl_io(:)                 = rmdi
nsw_io(:)                   = rmdi
nr_io(:)                    = rmdi
hw_sw_io(:)                 = rmdi
can_struct_a_io(:)          = rmdi
gsoil_f_io(:)               = rmdi
omegau_io(:)                = rmdi
omega_io(:)                 = rmdi
omegal_io(:)                = rmdi
omniru_io(:)                = rmdi
omnir_io(:)                 = rmdi
omnirl_io(:)                = rmdi
r_grow_io(:)                = rmdi
rootd_ft_io(:)              = rmdi
psi_close_io(:)             = rmdi
psi_open_io(:)              = rmdi
fsmc_p0_io(:)               = rmdi
sigl_io(:)                  = rmdi
tleaf_of_io(:)              = rmdi
tlow_io(:)                  = rmdi
tupp_io(:)                  = rmdi
emis_pft_io(:)              = rmdi
z0hm_pft_io(:)              = rmdi
z0hm_classic_pft_io(1:npft) = rmdi
dust_veg_scj_io(:)          = rmdi
fl_o3_ct_io(:)              = rmdi
dfp_dcuo_io(:)              = rmdi
ci_st_io(:)                 = rmdi
gpp_st_io(:)                = rmdi
ief_io(:)                   = rmdi
tef_io(:)                   = rmdi
mef_io(:)                   = rmdi
aef_io(:)                   = rmdi
q10_leaf_io(:)              = rmdi
lma_io(:)                   = rmdi
nmass_io(:)                 = rmdi
vsl_io(:)                   = rmdi
vint_io(:)                  = rmdi
kn_io(:)                    = rmdi
knl_io(:)                   = rmdi
fef_co2_io(:)               = rmdi
fef_co_io(:)                = rmdi
fef_ch4_io(:)               = rmdi
fef_nox_io(:)               = rmdi
fef_so2_io(:)               = rmdi
fef_oc_io(:)                = rmdi
fef_bc_io(:)                = rmdi
ccleaf_min_io(:)            = rmdi
ccleaf_max_io(:)            = rmdi
ccwood_min_io(:)            = rmdi
ccwood_max_io(:)            = rmdi
avg_ba_io(:)                = rmdi
fire_mort_io(:)             = rmdi

canht_ft_io(:) = rmdi
lai_io(:)      = rmdi

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(routinename, "Reading JULES_PFTPARM namelist...")

! Open the pft parameters namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'pft_params.nml'),          &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error opening namelist file pft_params.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_pftparm, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error reading namelist JULES_PFTPARM " //                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Error closing namelist file pft_params.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Radiation and albedo parameters.
orient(:)       = orient_io(1:npft)
albsnc_max(:)   = albsnc_max_io(1:npft)
albsnc_min(:)   = albsnc_min_io(1:npft)
albsnf_maxu(:)  = albsnf_maxu_io(1:npft)
albsnf_max(:)   = albsnf_max_io(1:npft)
albsnf_maxl(:)  = albsnf_maxl_io(1:npft)
alniru(:)       = alniru_io(1:npft)
alnir(:)        = alnir_io(1:npft)
alnirl(:)       = alnirl_io(1:npft)
alparu(:)       = alparu_io(1:npft)
alpar(:)        = alpar_io(1:npft)
alparl(:)       = alparl_io(1:npft)
kext(:)         = kext_io(1:npft)
kpar(:)         = kpar_io(1:npft)
lai_alb_lim(:)  = lai_alb_lim_io(1:npft)
can_struct_a(:) = can_struct_a_io(1:npft)
omegau(:)       = omegau_io(1:npft)
omega(:)        = omega_io(1:npft)
omegal(:)       = omegal_io(1:npft)
omniru(:)       = omniru_io(1:npft)
omnir(:)        = omnir_io(1:npft)
omnirl(:)       = omnirl_io(1:npft)

! Photosynthesis and respiration parameters.
c3(:)           = c3_io(1:npft)
act_jmax(:)     = act_jmax_io(1:npft)
act_vcmax(:)    = act_vcmax_io(1:npft)
alpha(:)        = alpha_io(1:npft)
alpha_elec(:)   = alpha_elec_io(1:npft)
deact_jmax(:)   = deact_jmax_io(1:npft)
deact_vcmax(:)  = deact_vcmax_io(1:npft)
dqcrit(:)       = dqcrit_io(1:npft)
ds_jmax(:)      = ds_jmax_io(1:npft)
ds_vcmax(:)     = ds_vcmax_io(1:npft)
fd(:)           = fd_io(1:npft)
f0(:)           = f0_io(1:npft)
g1_stomata(:)   = g1_stomata_io(1:npft)
jv25_ratio(:)   = jv25_ratio_io(1:npft)
kn(:)           = kn_io(1:npft)
knl(:)          = knl_io(1:npft)
neff(:)         = neff_io(1:npft)
nl0(:)          = nl0_io(1:npft)
nr_nl(:)        = nr_nl_io(1:npft)
ns_nl(:)        = ns_nl_io(1:npft)
r_grow(:)       = r_grow_io(1:npft)
tlow(:)         = tlow_io(1:npft)
tupp(:)         = tupp_io(1:npft)

! Trait physiology parameters
lma(:)          = lma_io(1:npft)
nmass(:)        = nmass_io(1:npft)
vsl(:)          = vsl_io(1:npft)
vint(:)         = vint_io(1:npft)
q10_leaf(:)     = q10_leaf_io(1:npft)
nr(:)           = nr_io(1:npft)
nsw(:)          = nsw_io(1:npft)
hw_sw(:)        = hw_sw_io(1:npft)

! Allometric and other parameters.
a_wl(:)         = a_wl_io(1:npft)
a_ws(:)         = a_ws_io(1:npft)
b_wl(:)         = b_wl_io(1:npft)
eta_sl(:)       = eta_sl_io(1:npft)
sigl(:)         = sigl_io(1:npft)

! Phenology parameters.
g_leaf_0(:)     = g_leaf_0_io(1:npft)
dgl_dm(:)       = dgl_dm_io(1:npft)
fsmc_of(:)      = fsmc_of_io(1:npft)
dgl_dt(:)       = dgl_dt_io(1:npft)
tleaf_of(:)     = tleaf_of_io(1:npft)

! Hydrological, thermal and other "physical" characteristics.
! Note that z0h_z0m, canht_ft and lai are read in via the pftparm namelist
! since they logically belong there, but the actual variables aren't in pftparm.
catch0(:)       = catch0_io(1:npft)
dcatch_dlai(:)  = dcatch_dlai_io(1:npft)
infil_f(:)      = infil_f_io(1:npft)
gsoil_f(:)      = gsoil_f_io(1:npft)
glmin(:)        = glmin_io(1:npft)
dz0v_dh(:)      = dz0v_dh_io(1:npft)
z0v(:)          = z0v_io(1:npft)
rootd_ft(:)     = rootd_ft_io(1:npft)
psi_close(:)    = psi_close_io(1:npft)
psi_open(:)     = psi_open_io(1:npft)
fsmc_p0(:)      = fsmc_p0_io(1:npft)
fsmc_mod(:)     = fsmc_mod_io(1:npft)
emis_pft(:)     = emis_pft_io(1:npft)
z0h_z0m(1:npft) = z0hm_pft_io(1:npft)
z0h_z0m_classic(1:npft) = z0hm_classic_pft_io(1:npft)
progs%canht_pft(:,:)   = SPREAD(canht_ft_io(1:npft), 1, land_pts)
progs%lai_pft(:,:)        = SPREAD(lai_io(1:npft), 1, land_pts)

! Ozone damage parameters.
fl_o3_ct(:)     = fl_o3_ct_io(1:npft)
dfp_dcuo(:)     = dfp_dcuo_io(1:npft)

! BVOC emission parameters.
ci_st(:)        = ci_st_io(1:npft)
gpp_st(:)       = gpp_st_io(1:npft)
ief(:)          = ief_io(1:npft)
tef(:)          = tef_io(1:npft)
mef(:)          = mef_io(1:npft)
aef(:)          = aef_io(1:npft)

! INFERNO emission parameters
fef_co2(:)      = fef_co2_io(1:npft)
fef_co(:)       = fef_co_io(1:npft)
fef_ch4(:)      = fef_ch4_io(1:npft)
fef_nox(:)      = fef_nox_io(1:npft)
fef_so2(:)      = fef_so2_io(1:npft)
fef_oc(:)       = fef_oc_io(1:npft)
fef_bc(:)       = fef_bc_io(1:npft)

! INFERNO combustion parameters
ccleaf_min(:)   = ccleaf_min_io(1:npft)
ccleaf_max(:)   = ccleaf_max_io(1:npft)
ccwood_min(:)   = ccwood_min_io(1:npft)
ccwood_max(:)   = ccwood_max_io(1:npft)
avg_ba(:)       = avg_ba_io(1:npft)
fire_mort(:)    = fire_mort_io(1:npft)

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to rmdi.
! Some configurations don't need all parameters but in some cases these are
! still tested below.
!-----------------------------------------------------------------------------
error = 0
IF ( ANY( orient(:) < 0 ) ) THEN  ! orient was initialised to < 0
  error = 1
  CALL log_error(routinename, "No value for orient")
END IF
IF ( ANY( fsmc_mod(:) < 0 ) ) THEN  ! fsmc_mod was initialised to < 0
  error = 1
  CALL log_error(routinename, "No value for fsmc_mod")
END IF
IF ( ANY( ABS( kext(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for kext")
END IF
IF ( ANY( ABS( kpar(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for kpar")
END IF
IF ( ANY( ABS( lai_alb_lim(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for lai_alb_lim")
END IF
IF ( ANY( ABS( can_struct_a(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for can_struct_a")
END IF
IF ( ANY( ABS( gsoil_f(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for gsoil_f")
END IF
IF ( ANY( c3(:) < 0 ) ) THEN  ! c3 was initialised to < 0
  error = 1
  CALL log_error(routinename, "No value for c3")
END IF
IF ( ANY( ABS( alpha(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for alpha")
END IF
IF ( ANY( ABS( fd(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for fd")
END IF
IF ( ANY( ABS( neff(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for neff")
END IF
IF ( ANY( ABS( nl0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for nl0")
END IF
IF ( ANY( ABS( nr_nl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for nr_nl")
END IF
IF ( ANY( ABS( ns_nl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for ns_nl")
END IF
IF ( ANY( ABS( r_grow(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for r_grow")
END IF
! Note that tlow and tupp are always required, for some PFTs at least.
! If using the Farquhar model for C3 plants, we still need tlow and
! tupp for C4 plants.
IF ( ANY( ABS( tlow(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for tlow")
END IF
IF ( ANY( ABS( tupp(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for tupp")
END IF

SELECT CASE ( photo_model )
CASE ( photo_farquhar )
  !---------------------------------------------------------------------------
  ! First check parameters that are always required with this model.
  !---------------------------------------------------------------------------
  ! Note that these parameter values are not used for C4 plants, but
  ! here we're still checking that they have been provided.
  IF ( ANY( ABS( act_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for act_jmax")
  END IF
  IF ( ANY( ABS( act_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for act_vcmax")
  END IF
  IF ( ANY( ABS( alpha_elec(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alpha_elec")
  END IF
  IF ( ANY( ABS( deact_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for act_jmax")
  END IF
  IF ( ANY( ABS( deact_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for act_vcmax")
  END IF
  IF ( ANY( ABS( jv25_ratio(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for jv25_ratio")
  END IF
  !---------------------------------------------------------------------------
  ! Check parameters that depend on any chosen acclimation model.
  !---------------------------------------------------------------------------
  IF ( photo_acclim_model == 0 ) THEN
    ! No acclimation.
    IF ( ANY( ABS( ds_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      error = 1
      CALL log_error("init_pftparm", "No value for ds_jmax")
    END IF
    IF ( ANY( ABS( ds_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      error = 1
      CALL log_error("init_pftparm", "No value for ds_vcmax")
    END IF
  END IF  !  photo_acclim_model == 0

END SELECT  !  photo_model

SELECT CASE ( stomata_model )
CASE ( stomata_jacobs )
  IF ( ANY( ABS( dqcrit(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dqcrit")
  END IF
  IF ( ANY( ABS( f0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for f0")
  END IF
CASE ( stomata_medlyn )
  IF ( ANY( ABS( g1_stomata(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for g1_stomata")
  END IF
END SELECT

IF ( .NOT. l_spec_albedo .AND. can_rad_mod == 1 ) THEN
  ! These don't need to be set
ELSE
  IF ( ANY( ABS( alnir(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alnir")
  END IF
  IF ( ANY( ABS( alpar(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alpar")
  END IF
  IF ( ANY( ABS( omega(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omega")
  END IF
  IF ( ANY( ABS( omnir(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omnir")
  END IF
END IF

IF ( l_albedo_obs ) THEN
  IF ( ANY( ABS( alnirl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alnirl")
  END IF
  IF ( ANY( ABS( alniru(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alniru")
  END IF
  IF ( ANY( ABS( alparl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alparl")
  END IF
  IF ( ANY( ABS( alparu(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for alparu")
  END IF
  IF ( ANY( ABS( omegal(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omegal")
  END IF
  IF ( ANY( ABS( omegau(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omegau")
  END IF
  IF ( ANY( ABS( omnirl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omnirl")
  END IF
  IF ( ANY( ABS( omniru(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for omniru")
  END IF
END IF

IF ( .NOT. l_spec_albedo ) THEN
  IF ( ANY( ABS( albsnf_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for albsnf_max")
  END IF
  IF ( l_albedo_obs ) THEN
    IF ( ANY( ABS( albsnf_maxl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      error = 1
      CALL log_error(routinename, "No value for albsnf_maxl")
    END IF
    IF ( ANY( ABS( albsnf_maxu(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      error = 1
      CALL log_error(routinename, "No value for albsnf_maxu")
    END IF
  END IF
END IF

IF ( .NOT. l_snow_albedo ) THEN
  IF ( ANY( ABS( albsnc_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for albsnc_max")
  END IF
  IF ( ANY( ABS( albsnc_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for albsnc_min")
  END IF
END IF

IF (l_trait_phys) THEN
  IF ( ANY( ABS( lma(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for lma")
  END IF
  IF ( ANY( ABS( nmass(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for nmass")
  END IF
  IF ( ANY( ABS( vsl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for vsl")
  END IF
  IF ( ANY( ABS( vint(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for vint")
  END IF
  IF ( ANY( ABS( nr(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for nr")
  END IF
  IF ( ANY( ABS( nsw(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for nsw")
  END IF
  IF ( ANY( ABS( hw_sw(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for hw_sw")
  END IF

END IF !l_trait_phys

IF ( ANY( ABS( kn(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for kn")
END IF
IF ( can_rad_mod == 6 .AND. ANY( ABS( knl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for knl")
END IF
IF ( ANY( ABS( q10_leaf(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for q10_leaf")
END IF
IF ( ANY( ABS( a_wl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for a_wl")
END IF
IF ( ANY( ABS( a_ws(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for a_ws")
END IF
IF ( ANY( ABS( b_wl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for b_wl")
END IF
IF ( ANY( ABS( eta_sl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for eta_sl")
END IF
IF ( ANY( ABS( sigl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for sigl")
END IF
IF ( ANY( ABS( g_leaf_0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for g_leaf_0")
END IF
IF ( ANY( ABS( dgl_dm(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for dgl_dm")
END IF
IF ( ANY( ABS( fsmc_of(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for fsmc_of")
END IF
IF ( ANY( ABS( dgl_dt(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for dgl_dt")
END IF
IF ( ANY( ABS( tleaf_of(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for tleaf_of")
END IF
IF ( ANY( ABS( catch0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for catch0")
END IF
IF ( ANY( ABS( dcatch_dlai(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for dcatch_dlai")
END IF
IF ( ANY( ABS( infil_f(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for infil_f")
END IF
IF ( ANY( ABS( glmin(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for glmin")
END IF
IF ( .NOT. l_spec_veg_z0) THEN
  IF ( ANY( ABS( dz0v_dh(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for dz0v_dh")
  END IF
ELSE
  IF ( ANY( ABS( z0v(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for z0v")
  END IF
END IF
IF ( ANY( ABS( rootd_ft(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for rootd_ft")
END IF

IF ( l_use_pft_psi ) THEN
  IF ( ANY( ABS( psi_close(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for psi_close")
  END IF
  IF ( ANY( ABS( psi_open(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for psi_open")
  END IF
END IF !l_use_pft_psi

IF ( l_bvoc_emis ) THEN
  IF ( ANY( ABS( ci_st(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ci_st")
  END IF
  IF ( ANY( ABS( gpp_st(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for gpp_st")
  END IF
  IF ( ANY( ABS( ief(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ief")
  END IF
  IF ( ANY( ABS( tef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for tef")
  END IF
  IF ( ANY( ABS( mef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for mef")
  END IF
  IF ( ANY( ABS( aef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for aef")
  END IF
END IF

IF ( l_inferno ) THEN
  IF ( ANY( ABS( fef_co2(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fef_co2")
  END IF
  IF ( ANY( ABS( fef_co(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fef_co")
  END IF
  IF ( ANY( ABS( fef_ch4(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fef_ch4")
  END IF
  IF ( ANY( ABS( fef_nox(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fef_nox")
  END IF
  IF ( ANY( ABS( fef_so2(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fef_so2")
  END IF
  IF ( ANY( ABS( ccleaf_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ccleaf_min")
  END IF
  IF ( ANY( ABS( ccleaf_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ccleaf_max")
  END IF
  IF ( ANY( ABS( ccwood_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ccwood_min")
  END IF
  IF ( ANY( ABS( ccwood_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for ccwood_max")
  END IF
  IF ( ANY( ABS( avg_ba(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for avg_ba")
  END IF
END IF

IF ( l_trif_fire ) THEN
  IF ( ANY( ABS( fire_mort(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fire_mort")
  END IF
END IF

IF ( l_o3_damage ) THEN
  IF ( ANY( ABS( fl_o3_ct(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for fl_o3_ct")
  END IF
  IF ( ANY( ABS( dfp_dcuo(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error(routinename, "No value for dfp_dcuo")
  END IF
END IF

IF ( ANY( ABS( fsmc_p0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for fsmc_p0")
END IF
IF ( ANY( ABS( emis_pft(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for emis_pft")
END IF
IF ( ANY( ABS( z0h_z0m(1:npft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for z0hm_pft")
END IF
IF ( ANY( ABS( z0h_z0m_classic(1:npft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for z0hm_classic_pft")
END IF
IF ( ANY( ABS( progs%canht_pft(:,:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for canht_ft")
END IF
IF ( ANY( ABS( progs%lai_pft(:,:) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error(routinename, "No value for lai")
END IF

IF ( error /= 0 )                                                             &
  CALL log_fatal(routinename,                                                 &
                 "Variable(s) missing from namelist - see earlier " //        &
                 "error message(s)")

!-----------------------------------------------------------------------------
! Check that glmin is >0.
! This ensures that wt_ext in subroutine soil_evap cannot become a NaN (which
! it would if gs=glmin and gsoil=0), or blow up, and might well be required
! elsewhere too.
!-----------------------------------------------------------------------------
IF ( ANY(glmin < 1.0e-10) )                                                   &
  CALL log_warn(routinename,                                                  &
                "Increasing one or more values of glmin - very small " //     &
                "values can cause model to blow up or NaNs")

WHERE ( glmin < 1.0e-10 )
  glmin = 1.0e-10
END WHERE

IF ( l_crop ) THEN
  IF ( ANY( ABS( a_ws(nnpft+1: npft) - 1.0 ) > EPSILON(1.0) ) ) THEN
    CALL log_fatal(routinename, "crop tiles should have a_ws=1.0")
  END IF
END IF

IF ( l_use_pft_psi ) THEN
  IF ( ANY( psi_close(1: npft) > EPSILON(1.0) ) ) THEN
    CALL log_fatal(routinename, "psi_close should be negative")
  END IF
  IF ( ANY( psi_open(1: npft) > EPSILON(1.0) ) ) THEN
    CALL log_fatal(routinename, "psi_open should be negative")
  END IF
END IF

!-----------------------------------------------------------------------------
! fsmc_mod=1 should not be allowed with a layered RothC  model until this has
! been properly evaluated. (With fsmc_mod=1, subroutine root_frac does not
! return the exponential root profile that users might expect.)
! Note that l_layeredC=T is not currently allowed with the UM.
!-----------------------------------------------------------------------------
IF ( l_layeredC .AND. ( soil_bgc_model == soil_model_rothc ) .AND.            &
     ANY( fsmc_mod(:) == 1 ) ) THEN
  CALL log_error(routinename,                                                 &
                 "fsmc_mod=1 is not allowed with l_layeredC and RothC")
END IF

RETURN

END SUBROUTINE init_pftparm_jules

#endif
