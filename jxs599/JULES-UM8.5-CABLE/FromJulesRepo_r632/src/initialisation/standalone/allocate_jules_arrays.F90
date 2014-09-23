#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE allocate_jules_arrays()

  USE switches, ONLY : l_triffid, l_phenol, l_albedo_obs, l_spec_albedo

  USE nstypes, ONLY : npft, nnvg, ntype

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE dust_parameters_mod, ONLY : ndiv

  USE ancil_info

  USE prognostics

  USE fluxes

  USE p_s_parms

  USE forcing

  USE screen

  USE pftparm

  USE nvegparm

  USE trif

  USE ozone_vars

  USE aero

  USE orog

  USE u_v_grid

  USE trifctl

  USE top_pdm

  USE coastal

  USE surf_param, ONLY : diff_frac

  USE soil_param, ONLY : dzsoil

  USE c_z0h_z0m, ONLY : z0h_z0m, z0h_z0m_classic

  USE c_elevate, ONLY : surf_hgt, z_land

  USE jules_mod, ONLY : snowdep_surf, albobs_scaling

  USE snow_param, ONLY : dzsnow, ds, cansnowtile

  USE bvoc_vars, ONLY : isoprene, isoprene_ft, terpene, terpene_ft,           &
                        methanol, methanol_ft, acetone, acetone_ft

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates the model arrays using sizes determined during initialisation
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: error, error_sum  ! Error indicators


!-----------------------------------------------------------------------------

! Initialise.
  error_sum = 0

! Runoff components.
  ALLOCATE( SUB_SURF_ROFF(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SURF_ROFF(land_pts), STAT=error )
  error_sum = error_sum + error

! Index variables
  ALLOCATE( TILE_INDEX(land_pts,ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SOIL_INDEX(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LICE_INDEX(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FRAC(land_pts,ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ICE_FRACT(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ICE_FRACT_NCAT(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Z1_UV(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Z1_TQ(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SSI_INDEX(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SEA_INDEX(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SICE_INDEX(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SICE_PTS_NCAT(nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FSSI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Screen variables
  ALLOCATE( Q1P5M(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Q1P5M_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( T1P5M(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( T1P5M_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( U10M(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( V10M(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Plant and soil parameters
  ALLOCATE( ALBSOIL(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ALBOBS_SW(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ALBOBS_VIS(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ALBOBS_NIR(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CATCH(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CATCH_SNOW(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( COSZ(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DIFF_FRAC(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( B(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DZSOIL(sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SATHH(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( HCON(land_pts,0:sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SMVCCL(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SMVCST(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SMVCWT(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( HCAP(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SATCON(land_pts,0:sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( INFIL_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Z0_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( STHU(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( STHF(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SOIL_CLAY(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Surface type variables.
  ALLOCATE( z0h_z0m(ntype), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0h_z0m_classic(ntype), STAT=error )
  error_sum = error_sum + error

! Veg surface type variables.
  ALLOCATE( nvgName(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( pftName(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnc_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnc_min(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_maxu(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_max(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_maxl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alpha(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alniru(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alnir(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alnirl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alparu(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alpar(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( alparl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( a_wl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( a_ws(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( b_wl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( c3(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dcatch_dlai(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dgl_dm(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dgl_dt(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dqcrit(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dz0v_dh(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emis_pft(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( eta_sl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fd(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc_of(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( f0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( glmin(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( g_leaf_0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( infil_f(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( kext(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( kpar(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( neff(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nl0(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( nr_nl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ns_nl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omegau(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omega(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omegal(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omniru(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omnir(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( omnirl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( orient(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( r_grow(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( rootd_ft(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( sigl(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tleaf_of(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tlow(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tupp(npft), STAT=error )
  error_sum = error_sum + error

! Ozone damage parameters
  ALLOCATE( fl_o3_ct(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( dfp_dcuo(npft), STAT=error )
  error_sum = error_sum + error

! BVOC emission parameters
  ALLOCATE( ief(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( tef(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( mef(npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( aef(npft), STAT=error )
  error_sum = error_sum + error

! Non-veg surface type variables.
  ALLOCATE( albsnc_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvgu(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( albsnf_nvgl(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( catch_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( emis_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( gs_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( infil_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( z0_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ch_nvg(nnvg), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( vf_nvg(nnvg), STAT=error )
  error_sum = error_sum + error

! Height above mean grid-box
  ALLOCATE( surf_hgt(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error

! Land height
  ALLOCATE( Z_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! TRIFFID variables - only needed if TRIFFID and/or phenology is selected.
  IF ( l_triffid .OR. l_phenol ) THEN
    ALLOCATE( crop(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_area(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_grow(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_root(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( g_wood(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lai_max(npft), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( lai_min(npft), STAT=error )
    error_sum = error_sum + error
  ENDIF

! Forcing variables
  ALLOCATE( QW_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TL_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( U_0(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( V_0(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( U_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( V_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( PSTAR(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LS_RAIN(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CON_RAIN(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LS_SNOW(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CON_SNOW(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SW_DOWN(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LW_DOWN(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( diff_rad(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Ozone variables
  ALLOCATE( o3(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( flux_o3_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fo3_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) o3(:) = 0.0

! BVOC diagnostics
  ALLOCATE( isoprene(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( isoprene_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( terpene(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( terpene_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( methanol(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( methanol_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( acetone(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( acetone_ft(land_pts,npft), STAT=error )
  error_sum = error_sum + error

! Prognostics variables
  ALLOCATE( LAI(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CANHT_FT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SMCL(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( T_SOIL(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RGRAIN(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RHO_SNOW_GRND(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOW_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SOOT(t_i_length*t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSTAR_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CANOPY(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CANOPY_GB(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CS(land_pts,DIM_CS1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Z0MSEA(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( GS(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( GC(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SMC(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DI_NCAT(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( K_SICE(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOW_GRND(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOW_MASS(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOW_MASS_SEA_NCAT(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( NSNOW(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOWDEPTH(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOWDEP_SURF(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  IF ( l_albedo_obs .AND. l_spec_albedo) THEN  
    ALLOCATE( albobs_scaling(land_pts,ntype,2), STAT=error )  
    error_sum = error_sum + error
  ELSE IF ( l_albedo_obs .AND. (.NOT.l_spec_albedo) ) THEN  
    ALLOCATE( albobs_scaling(land_pts,ntype,1), STAT=error )  
    error_sum = error_sum + error
  END IF  
! Snow layer variables.
! Some compilers seem to be unhappy passing these if they have not been
! allocated, others don't mind. So now allocating even if nsmax=0.
  ALLOCATE( DZSNOW(nsmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DS(land_pts,ntiles,nsmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RGRAINL(land_pts,ntiles,nsmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SICE(land_pts,ntiles,nsmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SLIQ(land_pts,ntiles,nsmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSNOW(land_pts,ntiles,nsmax), STAT=error )
  error_sum = error_sum + error

! (Forcing) fluxes
  ALLOCATE( ALB_TILE(land_pts,ntiles,4) )
  error_sum = error_sum + error
  ALLOCATE( TSTAR(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( E_SEA(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FQW_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( fsmc(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FTL_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FTL_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LE_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( H_SEA(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUX_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUY_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FQW_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FQW_ICE(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FTL_ICE(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ECAN(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ESOIL_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SEA_ICE_HTF(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SURF_HT_FLUX(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_htf_tile(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SICE_MLT_HTF(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOMLT_SURF_HTF(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LAND_ALBEDO(t_i_length,t_j_length,4) )
  error_sum = error_sum + error
  ALLOCATE( LATENT_HEAT(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( EI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( EI_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ECAN_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ESOIL(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( EXT(land_pts,sm_levels), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOWMELT(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( MELT_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( HF_SNOW_MELT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RADNET_TILE(land_pts,ntiles) )
  error_sum = error_sum + error
  ALLOCATE( SW_TILE(land_pts,ntiles) )
  error_sum = error_sum + error
  ALLOCATE( EMIS_TILE(land_pts,ntiles) )
  error_sum = error_sum + error
  ALLOCATE( SNOW_MELT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SNOMLT_SUB_HTF(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TOT_TFALL(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( surf_ht_store(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( anthrop_heat(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ALB_SICE(t_i_length*t_j_length,nice_use,4), STAT=error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    alb_tile(:,:,:) = 0.0
    ftl_tile(:,:) = 0.0
    le_tile(:,:) = 0.0
    fqw_tile(:,:) = 0.0
    esoil_tile(:,:) = 273.15
    ecan_tile(:,:) = 0.0
    ei_tile(:,:) = 0.0
    melt_tile(:,:) = 0.0
    radnet_tile(:,:) = 0.0
    sw_tile(:,:) = 0.0
    emis_tile(:,:) = 0.0
    snow_grnd(:,:) = 0.0
    surf_htf_tile(:,:) = 0.0
    surf_ht_store(:,:) = 0.0
    anthrop_heat(:,:) = 0.0
    alb_sice(:,:,:) = 0.0
  ENDIF

! Aerosol variables
  ALLOCATE( CO2_3D(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RHO_CD_MODV1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RHO_ARESIST(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ARESIST(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESIST_B(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RHO_ARESIST_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ARESIST_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESIST_B_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( R_B_DUST(t_i_length,t_j_length,NDIV), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CD_STD_DUST(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( U_S_STD_TILE(land_pts,ntiles), STAT=error )
  error_sum = error_sum + error

! Orographic roughness variables
  ALLOCATE( SIL_OROG_LAND(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( HO2R2_OROG(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( H_BLEND_OROG(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( Z0M_EFF(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Grid-change variables
  ALLOCATE( U_0_P(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( V_0_P(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( U_1_P(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( V_1_P(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DTRDZ_CHARNEY_GRID_1(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error

! Triffid variables
  ALLOCATE( G_LEAF_ACC(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( NPP_FT_ACC(land_pts_TRIF,npft_trif), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_W_FT_ACC(land_pts_TRIF,npft_trif), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_S_ACC(land_pts_TRIF,DIM_CS1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( G_LEAF_PHEN_ACC(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( GPP(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( NPP(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_P(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( G_LEAF(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( G_LEAF_PHEN(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( GPP_FT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( NPP_FT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_P_FT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_S(land_pts,DIM_CS1), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_W_FT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LAI_PHEN(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( C_VEG(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( CV(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( G_LEAF_DAY(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( G_LEAF_DR_OUT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LIT_C(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( LIT_C_MN(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( NPP_DR_OUT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_W_DR_OUT(land_pts,npft), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( RESP_S_DR_OUT(land_pts,5), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FRAC_AGR(land_pts), STAT=error )
  error_sum = error_sum + error

!       TOPMODEL and PDM variables
  ALLOCATE( A_FSAT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( A_FWET(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( C_FSAT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( C_FWET(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DRAIN(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( DUN_ROFF(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FCH4_WETL(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FEXP(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FSAT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( FWETL(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( GAMTOT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( QBASE(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( QBASE_ZW(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( STHZW(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TI_MEAN(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TI_SIG(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( ZW(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( INLANDOUT_ATM(land_pts), STAT=error )
  error_sum = error_sum + error

! Inititialise sthzw in all cases, although it is only required if l_top=T
  IF ( error_sum == 0 ) sthzw(:) = 0.0

! Coastal tiling variables
  ALLOCATE( TSTAR_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSTAR_SEA(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSTAR_SICE(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSTAR_SICE_NCAT(t_i_length,t_j_length,nice_use), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TSTAR_SSI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUX_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUX_SSI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUY_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( TAUY_SSI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( VSHR_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( VSHR_SSI(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SURF_HT_FLUX_LAND(t_i_length,t_j_length), STAT=error )
  error_sum = error_sum + error
  ALLOCATE( SURF_HT_FLUX_SICE(t_i_length,t_j_length,nice), STAT=error )
  error_sum = error_sum + error

! ANCIL variables.
  ALLOCATE( tile_pts(ntype), STAT=error )
  error_sum = error_sum + error

! Snow variables.
  ALLOCATE( canSnowTile(ntiles), STAT=error )
  error_sum = error_sum + error

! Check for error.
  IF ( error_sum /= 0 )                                                       &
    CALL log_fatal("allocate_jules_arrays",                                   &
                   "Error allocating JULES model arrays")

  RETURN

END SUBROUTINE allocate_jules_arrays
#endif
