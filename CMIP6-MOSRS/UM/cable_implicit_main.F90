!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose:
!
! Called from: JULES: surf_couple_ pathway
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_implicit_main_mod
  
contains

SUBROUTINE cable_implicit_main(                                                &        
              cycleno,                                                         & 
              row_length,rows, land_pts, ntiles, npft, sm_levels,              &
        dim_cs1, dim_cs2,                  & 
        Fland,                      &
!forcing: precip
        ls_rain_cable, conv_rain_cable,                                        &
        ls_snow_cable, conv_snow_cable, & !
!forcing: temp and humidity 
        tl_1, qw_1,&
!forcing: increments to temp and humidity 
        dtl1_1, dqw1_1,  &
!prog: canopy water storage
        canopy_gb, & ! aggregate over tiles, per land_pt 
        canopy,    & ! per tile
!prog: UM soil quantities, non-tiled  aggregate over CABLE tiles per land_pt
        T_SOIL, &
        smcl, &
        STHF,&
        STHU,&
        snow_surft, & 
        ftl_1,                   &
        ftl_surft, &
        fqw_1, &
        fqw_surft,  &
        tstar_surft, &
        surf_ht_flux_land,         &
        surf_htf_surft,     &
        ecan_surft, esoil_surft,                 &
        ei_surft, radnet_surft,                  &
        gs, &
        gs_surft, &
        t1p5m_surft, &
        q1p5m_surft, &
        melt_surft, &
        NPP, NPP_FT, GPP, GPP_FT,                                             &
        RESP_S, RESP_S_TOT, &   !RESP_S_TILE, !Kathy intro-ed as diag         &
        RESP_P, RESP_P_FT, G_LEAF,                                            &
        snow_depth                                                            &
)
  
  USE cable_implicit_driv_mod, ONLY :cable_implicit_driver

  USE cable_common_module, ONLY : knode_gl,        & ! processor number
                                  ktau_gl,         & ! number
                                  kwidth_gl          ! width in S 
  
  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable,           &
                                  soil_froz_frac_cable, snow_dpth_cable,       & 
                                  snow_mass_cable, snow_temp_cable,            &
                                  snow_rho_cable, snow_avg_rho_cable,          &   
                                  snow_age_cable, snow_flg_cable
  
  USE cable_gather_um_data_decs, ONLY : smvcst_cable  
  USE atmos_physics2_alloc_mod, ONLY : resp_s_tile


  USE atm_fields_real_mod, ONLY : C_pool_casa, N_pool_casa, P_pool_casa,       &
                                  SOIL_ORDER_casa, N_DEP_casa, N_FIX_casa,     &
                                  P_DUST_casa, P_weath_casa, LAI_casa,         &
                                  PHENPHASE_casa, NPP_PFT_ACC, RSP_W_PFT_ACC

 USE model_time_mod, ONLY: &
    target_end_stepim, i_day, i_day_number
     
  implicit none
 
  !--- IN ARGS FROM sf_impl2_cable, passed from surf_couple_implicit() down ----
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 
  
  REAL,  DIMENSION(land_pts) :: & 
      fland
   
  REAL, DIMENSION(row_length,rows) :: &
    ls_rain_cable,    &!forcing: precip:rain , large scale
    ls_snow_cable,    &!forcing: precip: snow, large scale
    conv_rain_cable,  &!forcing: precip:rain , convective 
    conv_snow_cable    !forcing: precip:snow, convective 

  REAL, DIMENSION(row_length,rows) :: &
   tl_1, qw_1,       &  !forcing: temp and humidity 
   dtl1_1, dqw1_1       !forcing: increments to temp and humidity 

  REAL ::                                                                       &
    canopy_gb(land_pts),         & !prog: canopy water store aggregate over tiles 
    canopy(land_pts, ntiles)       !prog:  per tile

!prog: UM soil quantities, non-tiled  aggregate over CABLE tiles per land_pt
real, dimension(land_pts,sm_levels) ::                           &
  T_SOIL, &
  smcl, &
  STHF,&
  STHU

real, dimension( land_pts, ntiles ) ::                           &
  snow_surft, &     ! snow ammount on tile nee:snow_tile 
  ftl_surft, &
  fqw_surft,  &
  tstar_surft, &
  surf_htf_surft, &
  ecan_surft, esoil_surft,                 &
  ei_surft, radnet_surft, &
  gs_surft, &
  t1p5m_surft, &
  q1p5m_surft, &
  melt_surft, &
  snow_depth
  
real, dimension( land_pts ) ::                           &
  gs
  
real, dimension(row_length,rows) ::      &
  ftl_1, &
  fqw_1, &
  surf_ht_flux_land

REAL ::                                                     &
   gpp(land_pts),                                                    &
                                ! IN Gross primary productivity
                                !    (kg C/m2/s).
   gpp_ft(land_pts,ntiles),                                            &
                                ! IN Gross primary productivity
                                !    on PFTs (kg C/m2/s).
   npp(land_pts),                                                    &
                                ! IN Net primary productivity
                                !    (kg C/m2/s).
   npp_ft(land_pts,ntiles),                                            &
                                ! IN Net primary productivity
                                !    on PFTs (kg C/m2/s).
   resp_p(land_pts),                                                 &
                                ! IN Plant respiration (kg C/m2/s).
   resp_p_ft(land_pts,ntiles),                                         &
                                  !IN Plant respiration on PFTs
                                !     (kg C/m2/s).
   resp_s(land_pts,dim_cs1),                                         &
                                   ! IN Soil respiration (kg C/m2/s).
   resp_s_tot(dim_cs2),                                               &
                                   ! IN Total soil resp'n (kg C/m2/s).
   g_leaf(land_pts,ntiles)
                                ! IN Leaf turnover rate (/360days).
 


  !--- End IN ARGS  -----------------------------------------------------------

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_implicit_main"
  logical, save :: first_call = .true.
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable
  
  !--- End header -------------------------------------------------------------
  
  if(knode_gl==0) then
    write (6, *) "CABLE_LSM:Start Subr: ", subr_name
    write (6, *) "@ knode_gl: ", knode_gl 
    write (6, *) "@ timestep: ", ktau_gl         
    write (6, *) "@ cycleno: ", cycleno         
    write (6, *) "============================="
  endif
     
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- Progress log and IN args @ timestep X,Y,Z                  -------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  
  isnow_flg_cable = int(snow_flg_cable)

 call cable_implicit_driver( i_day_number,                             &
cycleno,                                                         & 
row_length,rows, land_pts, ntiles, npft, sm_levels,              &
dim_cs1, dim_cs2,                                                & 
Fland,                                                           &
LS_RAIN_cable, CONV_RAIN_cable, LS_SNOW_cable, CONV_SNOW_cable,                           &
DTL1_1,DQW1_1, &
T_SOIL, &
soil_temp_cable, &
SMCL,        &
soil_moist_cable, &
real(kwidth_gl), &
SMVCST_cable,&
STHF, &
soil_froz_frac_cable, &
STHU,&
!STHU_TILE, &
snow_surft, &
snow_avg_rho_cable,       &
isnow_flg_cable, &
snow_dpth_cable, &
snow_mass_cable,      &
snow_rho_cable, &
snow_temp_cable,&
! this is diag only. initialize to huge(-1.)
!snow_cond, &
FTL_1, FTL_surft, FQW_1, FQW_surft,    &
TSTAR_surft, &
SURF_HT_FLUX_LAND, ECAN_surft, ESOIL_surft,    &
EI_surft, RADNET_surft, &
!TOT_ALB, &
SNOW_AGE_cable,   &
CANOPY, GS, &
gs_surft, &
T1P5M_surft, Q1P5M_surft,     &
CANOPY_GB, MELT_surft, &
NPP, NPP_FT, GPP, GPP_FT, RESP_S,   &
RESP_S_TOT, &
!Kathy intro'ed as a diag
RESP_S_TILE, &
RESP_P, RESP_P_FT,  &
G_LEAF, & 
TL_1, QW_1, &
SURF_HTF_surft, &
                              C_pool_casa, N_pool_casa, P_pool_casa,           &
                              LAI_casa, PHENPHASE_casa,            &
                              NPP_PFT_ACC, RSP_W_PFT_ACC )
 
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  snow_flg_cable = real(isnow_flg_cable)
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- OUT args @ timestep X,Y,Z                                  -------------
  !----------------------------------------------------------------------------

  !jhan: call checks as required by namelis      
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
    
  first_call = .false.        

return

End subroutine cable_implicit_main
  
End module cable_implicit_main_mod


