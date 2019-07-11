#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_pftparm(nml_dir)

  USE c_mdi, ONLY :                                                           &
!  imported scalar parameters
     rmdi

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE nstypes, ONLY : npft

  USE ancil_info, ONLY : land_pts

  USE pftparm

  USE c_z0h_z0m, ONLY : z0h_z0m, z0h_z0m_classic

  USE prognostics, ONLY : canht_ft, lai

  USE pftparm_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the PFT parameters
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


!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
  c3_io(:)     = -9
  orient_io(:) = -9

  a_wl_io(:)                  = RMDI
  a_ws_io(:)                  = RMDI
  albsnc_max_io(:)            = RMDI
  albsnc_min_io(:)            = RMDI
  albsnf_maxu_io(:)           = RMDI
  albsnf_max_io(:)            = RMDI
  albsnf_maxl_io(:)           = RMDI
  alpha_io(:)                 = RMDI
  alniru_io(:)                = RMDI
  alnir_io(:)                 = RMDI
  alnirl_io(:)                = RMDI
  alparu_io(:)                = RMDI
  alpar_io(:)                 = RMDI
  alparl_io(:)                = RMDI
  b_wl_io(:)                  = RMDI
  catch0_io(:)                = RMDI
  dcatch_dlai_io(:)           = RMDI
  dgl_dm_io(:)                = RMDI
  dgl_dt_io(:)                = RMDI
  dqcrit_io(:)                = RMDI
  dz0v_dh_io(:)               = RMDI
  eta_sl_io(:)                = RMDI
  fd_io(:)                    = RMDI
  fsmc_of_io(:)               = RMDI
  f0_io(:)                    = RMDI
  g_leaf_0_io(:)              = RMDI
  glmin_io(:)                 = RMDI
  infil_f_io(:)               = RMDI
  kext_io(:)                  = RMDI
  kpar_io(:)                  = RMDI
  neff_io(:)                  = RMDI
  nl0_io(:)                   = RMDI
  nr_nl_io(:)                 = RMDI
  ns_nl_io(:)                 = RMDI
  omegau_io(:)                = RMDI
  omega_io(:)                 = RMDI
  omegal_io(:)                = RMDI
  omniru_io(:)                = RMDI
  omnir_io(:)                 = RMDI
  omnirl_io(:)                = RMDI
  r_grow_io(:)                = RMDI
  rootd_ft_io(:)              = RMDI
  sigl_io(:)                  = RMDI
  tleaf_of_io(:)              = RMDI
  tlow_io(:)                  = RMDI
  tupp_io(:)                  = RMDI
  emis_pft_io(:)              = RMDI
  z0hm_pft_io(:)              = RMDI
  z0hm_classic_pft_io(1:npft) = RMDI
  dust_veg_scj_io(:)          = RMDI
  fl_o3_ct_io(:)              = RMDI
  dfp_dcuo_io(:)              = RMDI
  ief_io(:)                   = RMDI
  tef_io(:)                   = RMDI
  mef_io(:)                   = RMDI
  aef_io(:)                   = RMDI

  canht_ft_io(:) = RMDI
  lai_io(:)      = RMDI

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_pftparm", "Reading JULES_PFTPARM namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'pft_params.nml'),        &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_pftparm",                                            &
                   "Error opening namelist file pft_params.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_pftparm, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_pftparm",                                            &
                   "Error reading namelist JULES_PFTPARM " //                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_pftparm",                                            &
                   "Error closing namelist file pft_params.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
  pftname(:)      = pftname_io(1:npft)

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
  omegau(:)       = omegau_io(1:npft)
  omega(:)        = omega_io(1:npft)
  omegal(:)       = omegal_io(1:npft)
  omniru(:)       = omniru_io(1:npft)
  omnir(:)        = omnir_io(1:npft)
  omnirl(:)       = omnirl_io(1:npft)

! Photosynthesis and respiration parameters.
  c3(:)           = c3_io(1:npft)
  alpha(:)        = alpha_io(1:npft)
  dqcrit(:)       = dqcrit_io(1:npft)
  fd(:)           = fd_io(1:npft)
  f0(:)           = f0_io(1:npft)
  neff(:)         = neff_io(1:npft)
  nl0(:)          = nl0_io(1:npft)
  nr_nl(:)        = nr_nl_io(1:npft)
  ns_nl(:)        = ns_nl_io(1:npft)
  r_grow(:)       = r_grow_io(1:npft)
  tlow(:)         = tlow_io(1:npft)
  tupp(:)         = tupp_io(1:npft)

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
  glmin(:)        = glmin_io(1:npft)
  dz0v_dh(:)      = dz0v_dh_io(1:npft)
  rootd_ft(:)     = rootd_ft_io(1:npft)
  emis_pft(:)     = emis_pft_io(1:npft)
  z0h_z0m(1:npft) = z0hm_pft_io(1:npft)
  z0h_z0m_classic(1:npft) = z0hm_classic_pft_io(1:npft)
  canht_ft(:,:)   = SPREAD(canht_ft_io(1:npft), 1, land_pts)
  lai(:,:)        = SPREAD(lai_io(1:npft), 1, land_pts)

! Ozone damage parameters.
  fl_o3_ct(:)     = fl_o3_ct_io(1:npft)
  dfp_dcuo(:)     = dfp_dcuo_io(1:npft)

! BVOC emission parameters.
  ief(:)          = ief_io(1:npft)
  tef(:)          = tef_io(1:npft)
  mef(:)          = mef_io(1:npft)
  aef(:)          = aef_io(1:npft)

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to RMDI.
! Some configurations don't need all parameters but for now we insist on
! getting all parameters (and that there are not RMDI!).
!-----------------------------------------------------------------------------
  error = 0
  IF ( ANY( orient(:) < 0 ) ) THEN  ! orient was initialised to < 0
    error = 1
    CALL log_error("init_pftparm", "No value for orient")
  END IF
  IF ( ANY( ABS( albsnc_max(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for albsnc_max")
  END IF
  IF ( ANY( ABS( albsnc_min(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for albsnc_min")
  END IF
  IF ( ANY( ABS( albsnf_maxu(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for albsnf_maxu")
  END IF
  IF ( ANY( ABS( albsnf_max(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for albsnf_max")
  END IF
  IF ( ANY( ABS( albsnf_maxl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for albsnf_maxl")
  END IF
  IF ( ANY( ABS( alniru(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alniru")
  END IF
  IF ( ANY( ABS( alnir(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alnir")
  END IF
  IF ( ANY( ABS( alnirl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alnirl")
  END IF
  IF ( ANY( ABS( alparu(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alparu")
  END IF
  IF ( ANY( ABS( alpar(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alpar")
  END IF
  IF ( ANY( ABS( alparl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alparl")
  END IF
  IF ( ANY( ABS( kext(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for kext")
  END IF
  IF ( ANY( ABS( kpar(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for kpar")
  END IF
  IF ( ANY( ABS( omegau(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omegau")
  END IF
  IF ( ANY( ABS( omega(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omega")
  END IF
  IF ( ANY( ABS( omegal(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omegal")
  END IF
  IF ( ANY( ABS( omniru(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omniru")
  END IF
  IF ( ANY( ABS( omnir(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omnir")
  END IF
  IF ( ANY( ABS( omnirl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for omnirl")
  END IF
  IF ( ANY( c3(:) < 0 ) ) THEN  ! c3 was initialised to < 0
    error = 1
    CALL log_error("init_pftparm", "No value for c3")
  END IF
  IF ( ANY( ABS( alpha(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for alpha")
  END IF
  IF ( ANY( ABS( dqcrit(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dqcrit")
  END IF
  IF ( ANY( ABS( fd(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for fd")
  END IF
  IF ( ANY( ABS( f0(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for f0")
  END IF
  IF ( ANY( ABS( neff(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for neff")
  END IF
  IF ( ANY( ABS( nl0(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for nl0")
  END IF
  IF ( ANY( ABS( nr_nl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for nr_nl")
  END IF
  IF ( ANY( ABS( ns_nl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for ns_nl")
  END IF
  IF ( ANY( ABS( r_grow(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for r_grow")
  END IF
  IF ( ANY( ABS( tlow(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for tlow")
  END IF
  IF ( ANY( ABS( tupp(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for tupp")
  END IF
  IF ( ANY( ABS( a_wl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for a_wl")
  END IF
  IF ( ANY( ABS( a_ws(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for a_ws")
  END IF
  IF ( ANY( ABS( b_wl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for b_wl")
  END IF
  IF ( ANY( ABS( eta_sl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for eta_sl")
  END IF
  IF ( ANY( ABS( sigl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for sigl")
  END IF
  IF ( ANY( ABS( g_leaf_0(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for g_leaf_0")
  END IF
  IF ( ANY( ABS( dgl_dm(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dgl_dm")
  END IF
  IF ( ANY( ABS( fsmc_of(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for fsmc_of")
  END IF
  IF ( ANY( ABS( dgl_dt(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dgl_dt")
  END IF
  IF ( ANY( ABS( tleaf_of(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for tleaf_of")
  END IF
  IF ( ANY( ABS( catch0(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for catch0")
  END IF
  IF ( ANY( ABS( dcatch_dlai(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dcatch_dlai")
  END IF
  IF ( ANY( ABS( infil_f(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for infil_f")
  END IF
  IF ( ANY( ABS( glmin(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for glmin")
  END IF
  IF ( ANY( ABS( dz0v_dh(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dz0v_dh")
  END IF
  IF ( ANY( ABS( rootd_ft(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for rootd_ft")
  END IF
  IF ( ANY( ABS( emis_pft(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for emis_pft")
  END IF
  IF ( ANY( ABS( z0h_z0m(1:npft) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for z0hm_pft")
  END IF
  IF ( ANY( ABS( z0h_z0m_classic(1:npft) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for z0hm_classic_pft")
  END IF
  IF ( ANY( ABS( canht_ft(:,:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for canht_ft")
  END IF
  IF ( ANY( ABS( lai(:,:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for lai")
  END IF
  IF ( ANY( ABS( fl_o3_ct(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for fl_o3_ct")
  END IF
  IF ( ANY( ABS( dfp_dcuo(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for dfp_dcuo")
  END IF
  IF ( ANY( ABS( ief(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for ief")
  END IF
  IF ( ANY( ABS( tef(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for tef")
  END IF
  IF ( ANY( ABS( mef(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for mef")
  END IF
  IF ( ANY( ABS( aef(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_pftparm", "No value for aef")
  END IF

  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_pftparm",                                            &
                   "Variable(s) missing from namelist - see earlier " //      &
                   "error message(s)")

!-----------------------------------------------------------------------------
! Check that glmin is >0.
! This ensures that wt_ext in subroutine soil_evap cannot become a NaN (which
! it would if gs=glmin and gsoil=0), or blow up, and might well be required
! elsewhere too.
!-----------------------------------------------------------------------------
  IF ( ANY(glmin < 1.0e-10) )                                                 &
    CALL log_warn("init_pftparm",                                             &
                  "Increasing one or more values of glmin - very small " //   &
                  "values can cause model to blow up or NaNs")

  WHERE ( glmin < 1.0e-10 )
    glmin = 1.0e-10
  END WHERE

  RETURN

END SUBROUTINE init_pftparm
#endif
