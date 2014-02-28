#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/params/init_pftparm.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

SUBROUTINE init_pftparm()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE nstypes, ONLY : npft

  USE ancil_info, ONLY : land_pts

  USE pftparm

  USE c_z0h_z0m, ONLY : z0h_z0m

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
! Work variables
  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------


  CALL log_info("init_pftparm", "Reading JULES_PFTPARM namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE='pft_params.nml',                                  &
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
  c3(:)           = c3_io(1:npft)
  orient(:)       = orient_io(1:npft)
  a_wl(:)         = a_wl_io(1:npft)
  a_ws(:)         = a_ws_io(1:npft)
  albsnc_max(:)   = albsnc_max_io(1:npft)
  albsnc_min(:)   = albsnc_min_io(1:npft)
  albsnf_max(:)   = albsnf_max_io(1:npft)
  alpha(:)        = alpha_io(1:npft)
  alnir(:)        = alnir_io(1:npft)
  alpar(:)        = alpar_io(1:npft)
  b_wl(:)         = b_wl_io(1:npft)
  catch0(:)       = catch0_io(1:npft)
  dcatch_dlai(:)  = dcatch_dlai_io(1:npft)
  dgl_dm(:)       = dgl_dm_io(1:npft)
  dgl_dt(:)       = dgl_dt_io(1:npft)
  dqcrit(:)       = dqcrit_io(1:npft)
  dz0v_dh(:)      = dz0v_dh_io(1:npft)
  eta_sl(:)       = eta_sl_io(1:npft)
  fd(:)           = fd_io(1:npft)
  fsmc_of(:)      = fsmc_of_io(1:npft)
  f0(:)           = f0_io(1:npft)
  g_leaf_0(:)     = g_leaf_0_io(1:npft)
  glmin(:)        = glmin_io(1:npft)
  infil_f(:)      = infil_f_io(1:npft)
  kext(:)         = kext_io(1:npft)
  kpar(:)         = kpar_io(1:npft)
  neff(:)         = neff_io(1:npft)
  nl0(:)          = nl0_io(1:npft)
  nr_nl(:)        = nr_nl_io(1:npft)
  ns_nl(:)        = ns_nl_io(1:npft)
  omega(:)        = omega_io(1:npft)
  omnir(:)        = omnir_io(1:npft)
  r_grow(:)       = r_grow_io(1:npft)
  rootd_ft(:)     = rootd_ft_io(1:npft)
  sigl(:)         = sigl_io(1:npft)
  tleaf_of(:)     = tleaf_of_io(1:npft)
  tlow(:)         = tlow_io(1:npft)
  tupp(:)         = tupp_io(1:npft)
  emis_pft(:)     = emis_pft_io(1:npft)
  fl_o3_ct(:)     = fl_o3_ct_io(1:npft)
  dfp_dcuo(:)     = dfp_dcuo_io(1:npft)
  ief(:)          = ief_io(1:npft)
  tef(:)          = tef_io(1:npft)
  mef(:)          = mef_io(1:npft)
  aef(:)          = aef_io(1:npft)

! Variables that are read in via the pftparm namelist since they logically
! belong there, but the actual variables aren't in pftparm
  z0h_z0m(1:npft) = z0hm_pft_io(1:npft)
  canht_ft(:,:) = SPREAD(canht_ft_io(1:npft), 1, land_pts)
  lai(:,:)      = SPREAD(lai_io(1:npft), 1, land_pts)

! Check that glmin is >0.
! This ensures that wt_ext in subroutine soil_evap cannot become a NaN (which
! it would if gs=glmin and gsoil=0), or blow up, and might well be required
! elsewhere too
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
