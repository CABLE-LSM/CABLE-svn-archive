! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in pftparm data
! and initialisations

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE pftparm_io

  USE max_dimensions, ONLY:                                           &
    npft_max

  IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in pftparm that we want to initialise).
!---------------------------------------------------------------------
  INTEGER ::                                                          &
    c3_io(npft_max),                                                  &
    orient_io(npft_max)

  REAL ::                                                             &
    a_wl_io(npft_max),                                                &
    a_ws_io(npft_max),                                                &
    albsnc_max_io(npft_max),                                          &
    albsnc_min_io(npft_max),                                          &
    albsnf_maxu_io(npft_max),                                         &
    albsnf_max_io(npft_max),                                          &
    albsnf_maxl_io(npft_max),                                         &
    alpha_io(npft_max),                                               &
    alniru_io(npft_max),                                              &
    alnir_io(npft_max),                                               &
    alnirl_io(npft_max),                                              &
    alparu_io(npft_max),                                              &
    alpar_io(npft_max),                                               &
    alparl_io(npft_max),                                              &
    b_wl_io(npft_max),                                                &
    catch0_io(npft_max),                                              &
    dcatch_dlai_io(npft_max),                                         &
    dgl_dm_io(npft_max),                                              &
    dgl_dt_io(npft_max),                                              &
    dqcrit_io(npft_max),                                              &
    dz0v_dh_io(npft_max),                                             &
    eta_sl_io(npft_max),                                              &
    fd_io(npft_max),                                                  &
    fsmc_of_io(npft_max),                                             &
    f0_io(npft_max),                                                  &
    g_leaf_0_io(npft_max),                                            &
    glmin_io(npft_max),                                               &
    infil_f_io(npft_max),                                             &
    kext_io(npft_max),                                                &
    kpar_io(npft_max),                                                &
    neff_io(npft_max),                                                &
    nl0_io(npft_max),                                                 &
    nr_nl_io(npft_max),                                               &
    ns_nl_io(npft_max),                                               &
    omegau_io(npft_max),                                              &
    omega_io(npft_max),                                               &
    omegal_io(npft_max),                                              &
    omniru_io(npft_max),                                              &
    omnir_io(npft_max),                                               &
    omnirl_io(npft_max),                                              &
    r_grow_io(npft_max),                                              &
    rootd_ft_io(npft_max),                                            &
    sigl_io(npft_max),                                                &
    tleaf_of_io(npft_max),                                            &
    tlow_io(npft_max),                                                &
    tupp_io(npft_max),                                                &
    emis_pft_io(npft_max),                                            &
    z0hm_pft_io(npft_max),                                            &
    z0hm_classic_pft_io(npft_max),                                    &
    dust_veg_scj_io(npft_max),                                        &
    fl_o3_ct_io(npft_max),                                            &
    dfp_dcuo_io(npft_max),                                            &
    ief_io(npft_max),                                                 &
    tef_io(npft_max),                                                 &
    mef_io(npft_max),                                                 &
    aef_io(npft_max)

#if !defined(UM_JULES)
  CHARACTER(LEN=20) :: pftname_io(npft_max)

  REAL ::                                                             &
    canht_ft_io(npft_max),                                            &
    lai_io(npft_max)
#endif

!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
  NAMELIST /jules_pftparm/                                              &
#if !defined(UM_JULES)
                           pftname_io,canht_ft_io,lai_io,               &
#endif
                           c3_io,orient_io,a_wl_io,a_ws_io,             &
                           albsnc_max_io,albsnc_min_io,albsnf_maxu_io,  & 
                           albsnf_max_io,albsnf_maxl_io,                & 
                           alpha_io,alniru_io,alnir_io,alnirl_io,       & 
                           alparu_io,alpar_io,alparl_io,b_wl_io,        & 
                           catch0_io,dcatch_dlai_io,dgl_dm_io,          &
                           dgl_dt_io,dqcrit_io,dz0v_dh_io,eta_sl_io,    &
                           fd_io,fsmc_of_io,f0_io,g_leaf_0_io,          &
                           glmin_io,infil_f_io,kext_io,kpar_io,         &
                           neff_io,nl0_io,nr_nl_io,ns_nl_io,            & 
                           omegau_io,omega_io,omegal_io,omniru_io,      & 
                           omnir_io,omnirl_io,r_grow_io,rootd_ft_io,    & 
                           sigl_io,tleaf_of_io,tlow_io,tupp_io,         & 
                           emis_pft_io,z0hm_pft_io,z0hm_classic_pft_io, &
                           dust_veg_scj_io,fl_o3_ct_io,dfp_dcuo_io,     &
                           ief_io,tef_io,mef_io,aef_io

CONTAINS

  SUBROUTINE print_nlist_jules_pftparm()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer

    CALL jules_print('pftparm_io', &
        'Contents of namelist jules_pftparm')

#if !defined(UM_JULES)
    WRITE(lineBuffer,*)' pftname_io = ',pftname_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' canht_ft_io = ',canht_ft_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' lai_io = ',lai_io
    CALL jules_print('pftparm_io',lineBuffer)
#endif
    WRITE(lineBuffer,*)' c3_io = ',c3_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' orient_io = ',orient_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' a_wl_io = ',a_wl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' a_ws_io = ',a_ws_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' albsnc_max_io = ',albsnc_max_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' albsnc_min_io = ',albsnc_min_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' albsnf_max_io = ',albsnf_max_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' alpha_io = ',alpha_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' alnir_io = ',alnir_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' alpar_io = ',alpar_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' b_wl_io = ',b_wl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' catch0_io = ',catch0_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dcatch_dlai_io = ',dcatch_dlai_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dgl_dm_io = ',dgl_dm_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dgl_dt_io = ',dgl_dt_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dqcrit_io = ',dqcrit_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dz0v_dh_io = ',dz0v_dh_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' eta_sl_io = ',eta_sl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' fd_io = ',fd_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' fsmc_of_io = ',fsmc_of_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' f0_io = ',f0_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' g_leaf_0_io = ',g_leaf_0_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' glmin_io = ',glmin_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' infil_f_io = ',infil_f_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' kext_io = ',kext_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' kpar_io = ',kpar_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' neff_io = ',neff_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' nl0_io = ',nl0_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' nr_nl_io = ',nr_nl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' ns_nl_io = ',ns_nl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' omega_io = ',omega_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' omnir_io = ',omnir_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' r_grow_io = ',r_grow_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' rootd_ft_io = ',rootd_ft_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' sigl_io = ',sigl_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' tleaf_of_io = ',tleaf_of_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' tlow_io = ',tlow_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' tupp_io = ',tupp_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' emis_pft_io = ',emis_pft_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' z0hm_pft_io = ',z0hm_pft_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dust_veg_scj_io = ',dust_veg_scj_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' fl_o3_ct_io = ',fl_o3_ct_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' dfp_dcuo_io = ',dfp_dcuo_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' ief_io = ',ief_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' tef_io = ',tef_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' mef_io = ',mef_io
    CALL jules_print('pftparm_io',lineBuffer)
    WRITE(lineBuffer,*)' aef_io = ',aef_io
    CALL jules_print('pftparm_io',lineBuffer)

    CALL jules_print('pftparm_io', &
        '- - - - - - end of namelist - - - - - -')

  END SUBROUTINE print_nlist_jules_pftparm

END MODULE pftparm_io
