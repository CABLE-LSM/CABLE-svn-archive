! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL-------------------------------------------------

! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.

! Documentation : UM Documentation Paper 25

! Subroutine Interface:
MODULE hydrol_cable_mod
CONTAINS
SUBROUTINE hydrol_cable (                                         &
                   lice_pts,lice_index,soil_pts,soil_index,       &
                   nsnow,                                         &
                   npnts,nshyd,b,can_cpy,con_rain,                &
                   e_canopy,ext,hcap,hcon,ls_rain,                &
                   con_rainfrac, ls_rainfrac,                     &
                   satcon,sathh,snowdepth,                        &
                   surf_ht_flux,timestep,                         &
                   v_sat,v_wilt,                                  &
                   can_wcnt,stf_hf_snow_melt,                     &
                   stf_sub_surf_roff,smcl,sthf,sthu,tsoil,        &
                   can_wcnt_gb,hf_snow_melt,smc,                  &
                   snow_melt,                                     &
                   sub_surf_roff,surf_roff,tot_tfall,             &
! add new inland basin variable
                   inlandout_atm,l_inland,                        &
! Additional variables for MOSES II
                   ntiles,tile_pts,tile_index,                    &
                   infil_tile,                                    &
                   melt_tile,tile_frac,                           &
! Additional variables required for large-scale hydrology:
                   l_top,l_pdm,fexp,gamtot,ti_mean,ti_sig,cs,     &
                   dun_roff,drain,fsat,fwetl,qbase,qbase_zw,      &
                   zw,sthzw,a_fsat,c_fsat,a_fwet,c_fwet,          &
                   fch4_wetl,dim_cs1,l_soil_sat_down,l_triffid)

USE jules_hydrology_mod, ONLY :                                   &
 ti_max                                                           &
,zw_max

USE jules_soil_mod, ONLY :                                        &
 dzsoil       !  Thicknesses of the soil layers (m)

USE crop_vars_mod, ONLY : sthu_irr, frac_irr, ext_irr
             
USE jules_vegetation_mod, ONLY : l_irrig_dmd 

USE c_densty, ONLY :                                              &
   rho_water  !  density of pure water (kg/m3)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                            &
 lice_pts                                                         &
                     ! IN Number of land ice points.
,npnts                                                            &
                     ! IN Number of gridpoints.
,nshyd                                                            &
                     ! IN Number of soil moisture levels.
,soil_pts                                                         &
                     ! IN Number of soil points.
,ntiles                                                           &
                     ! IN Number of tiles
,dim_cs1             ! IN Number of soil carbon pools

REAL, INTENT(IN) ::                                               &
 timestep            ! IN Model timestep (s).

LOGICAL , INTENT(IN) ::                                           &
 stf_hf_snow_melt                                                 &
                     ! IN Stash flag for snowmelt heat flux.
!cxyz STF_HF_SNOW_MELT is not used in this version.
,stf_sub_surf_roff                                                &
                     ! IN Stash flag for sub-surface runoff.
,l_top                                                            &
             ! IN Flag for TOPMODEL-based hydrology.
,l_pdm                                                            &
             ! IN Flag for PDM hydrology.
,l_soil_sat_down                                                  &
             ! IN Switch controlling direction of movement of
             !    soil moisture in excess of saturation
,l_triffid   ! IN Switch to use TRIFFID.


!   Array arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                            &
 lice_index(npnts)                                                &
                     ! IN Array of land ice points.
,soil_index(npnts)                                                &
                     ! IN Array of soil points.
,nsnow(npnts,ntiles) ! IN Number of snow layers

REAL, INTENT(IN) ::                                               &
 b(npnts,nshyd)                                                   &
                     ! IN Clapp-Hornberger exponent.
,can_cpy(npnts,ntiles)                                            &
                      !IN Canopy/surface capacity of
!                          !    land tiles (kg/m2).
,con_rain(npnts)                                                  &
                     ! IN Convective rain (kg/m2/s).
,e_canopy(npnts,ntiles)                                           &
!                          ! IN Canopy evaporation from
!                          !    land tiles (kg/m2/s).
,ext(npnts,nshyd)                                                 &
                     ! IN Extraction of water from each soil
!                          !    layer (kg/m2/s).
,hcap(npnts,nshyd)                                                &
                     ! IN Soil heat capacity (J/K/m3).
,hcon(npnts,0:nshyd)                                              &
                     ! IN Soil thermal conductivity (W/m/K).
,ls_rain(npnts)                                                   &
                     ! IN Large-scale rain (kg/m2/s).
,con_rainfrac(npnts)                                              &
                     ! IN Convective rain fraction
,ls_rainfrac(npnts)                                               &
                     ! IN large scale rain fraction

,satcon(npnts,0:nshyd)                                            &
                     ! IN Saturated hydraulic conductivity
!                          !    (kg/m2/s).
,sathh(npnts,nshyd)                                               &
                     ! IN Saturated soil water pressure (m).
,snow_melt(npnts)                                                 &
                     ! IN Snowmelt (kg/m2/s).
,snowdepth(npnts,ntiles)                                          &
                     ! Snow depth (on ground) (m)
,surf_ht_flux(npnts)                                              &
                     ! IN Net downward surface heat flux (W/m2)
,v_sat(npnts,nshyd)                                               &
                     ! IN Volumetric soil moisture
!                          !    concentration at saturation
!                          !    (m3 H2O/m3 soil).
,v_wilt(npnts,nshyd)                                              &
                     ! IN Volumetric soil moisture
!                          !    concentration below which
!                          !    stomata close (m3 H2O/m3 soil).
,fexp(npnts)                                                      &
                     ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
,gamtot(npnts)                                                    &
                     ! IN Integrated complete Gamma function.
,ti_mean(npnts)                                                   &
                     ! IN Mean topographic index.
,ti_sig(npnts)                                                    &
                     ! IN Standard dev. of topographic index.
,cs(npnts,dim_cs1)                                                &
                     ! IN Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !    are DPM, RPM, biomass and humus.
,a_fsat(npnts)                                                    &
                     ! IN Fitting parameter for Fsat in LSH model
,c_fsat(npnts)                                                    &
                     ! IN Fitting parameter for Fsat in LSH model
,a_fwet(npnts)                                                    &
                     ! IN Fitting parameter for Fwet in LSH model
,c_fwet(npnts)
                     ! IN Fitting parameter for Fwet in LSH model


!   Array arguments with intent(INOUT) :

REAL, INTENT(INOUT) ::                                            &
 can_wcnt(npnts,ntiles)                                           &
!                          ! INOUT Canopy water content for
!                          !       land tiles (kg/m2).
,smcl(npnts,nshyd)                                                &
                     ! INOUT Soil moisture content of each
!                          !       layer (kg/m2).
,sthf(npnts,nshyd)                                                &
                     ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,sthu(npnts,nshyd)                                                &
                     ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
,tsoil(npnts,nshyd)                                               &
                     ! INOUT Sub-surface temperatures (K).
,fsat(npnts)                                                      &
                      ! INOUT Surface saturation fraction.
,fwetl(npnts)                                                     &
                      ! INOUT Wetland fraction.
,zw(npnts)                                                        &
                      ! INOUT Water table depth (m).
,sthzw(npnts)         ! INOUT soil moist fract. in deep-zw layer.

! Arguments which are required for compatibility but not used.
! Shouldnt use OUT.
REAL              ::                                              &
 hf_snow_melt(npnts)
                      ! OUT Gridbox snowmelt heat flux (W/m2).


!   Array arguments with intent(OUT) :
REAL, INTENT(OUT) ::                                              &
 can_wcnt_gb(npnts)                                               &
                      ! OUT Gridbox canopy water content (kg/m2).
,smc(npnts)                                                       &
                      ! OUT Available soil moisture in a layer at
                      !     the surface (kg/m2)
,sub_surf_roff(npnts)                                             &
                      ! OUT Sub-surface runoff (kg/m2/s).
,surf_roff(npnts)                                                 &
                      ! OUT Surface runoff (kg/m2/s).
,tot_tfall(npnts)                                                 &
                      ! OUT Total throughfall (kg/m2/s).
,dun_roff(npnts)                                                  &
                      ! OUT Dunne part of sfc runoff (kg/m2/s).
,qbase(npnts)                                                     &
                      ! OUT Base flow (kg/m2/s).
,qbase_zw(npnts)                                                  &
                      ! OUT Base flow from ZW layer (kg/m2/s).
,drain(npnts)                                                     &
                      ! OUT Drainage out of nshyd'th level (kg/m2/s).
,fch4_wetl(npnts)     ! OUT Scaled wetland methane flux.
                      !     (10^-9 kg C/m2/s).


! Additional variables for MOSES II
INTEGER, INTENT(IN) ::                                            &
 tile_pts(ntiles)                                                 &
                     ! IN Number of tile points.
,tile_index(npnts,ntiles)
!                          ! IN Index of tile points.

REAL, INTENT(IN) ::                                               &
 infil_tile(npnts,ntiles)                                         &
!                          ! IN Maximum surface infiltration
,melt_tile(npnts,ntiles)                                          &
!                          ! IN Snowmelt on tiles (kg/m2/s).
,tile_frac(npnts,ntiles)                                          &
                     ! IN Tile fractions.

! Declare variable for inland basin outflow
,inlandout_atm(npnts)            ! IN TRIP INLAND BASIN
!                       OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm
LOGICAL, INTENT(IN) ::                                            &
 l_inland                   ! IN True if re-routing inland
                            !   basin flow to soil moisture

! Local scalars:
INTEGER                                                           &
 i,j                                                              &
                      ! WORK Loop counters.
,n                    ! WORK Tile loop counter.

! Local arrays:

REAL                                                              &
 dsmc_dt(npnts)                                                   &
                      ! WORK Rate of change of soil moisture
!                           !      due to water falling onto the
!                           !      surface after surface runoff
!                           !      (kg/m2/s).
,w_flux(npnts,0:nshyd)                                            &
                      ! WORK Fluxes of water between layers
!                           !      (kg/m2/s).
,ksz(npnts,0:nshyd)                                               &
                      ! WORK Saturated hydraulic
!                           !      conductivity in layer (kg/m2/s).
,qbase_l(npnts,nshyd+1)                                           &
!                           ! WORK Base flow from each level (kg/m2/s).
,top_crit(npnts)                                                  &
                      ! WORK Critical TI when ZW <=0.0
,zdepth(0:nshyd)                                                  &
                      ! WORK Lower soil layer boundary depth (m).
,tsoil_d(npnts)                                                   &
                      ! WORK Soil temperature in the top metre
,wutot(npnts)         ! WORK Ratio of unfrozen to total soil
!                                             !    moisture at ZW.

! WORK variables required for irrigation code
REAL                                                              &
w_flux_irr(npnts,0:nshyd)                                         &
                      ! WORK The fluxes of water between layers
!                          !     in irrigated fraction (kg/m2/s).
,w_flux_nir(npnts,0:nshyd)                                        &
                      ! WORK The fluxes of water between layers
!                          !     in non-irrigated fraction (kg/m2/s).
,smcl_irr(npnts,nshyd)                                            &
                      ! WORK Total soil moisture contents of each
!                          !       layer in irrigated fraction (kg/m2).
,smcl_nir(npnts,nshyd)                                            &
                      ! WORK Total soil moisture contents of each
!                          !       layer in non-irrigated fraction (kg/m2).
,sthu_nir(npnts,nshyd)                                            &
                      ! WORK Unfrozen soil moisture content of
!                          !    each layer as a fraction of
!                          !    saturation in irrigated fraction.
,ext_nir(npnts,nshyd)                                             &
                      ! WORK Extraction of water from each soil
!                          !    layer in non-irrigated fraction (kg/m2/s).
,smclsat(npnts,nshyd)                                             &
                      ! WORK The saturation moisture content of
!                           !     each layer (kg/m2).
,smclzw(npnts)                                                    &             
                      ! WORK moisture content in deep layer.
,smclsatzw(npnts)                                                 
                      ! WORK moisture content in deep layer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header--------------------------------------------------------
IF (lhook) CALL dr_hook('HYDROL_CABLE',zhook_in,zhook_handle)

! Initialise w_flux variables that are used in irrigation code
  w_flux(:,:) = 0.0
  w_flux_irr(:,:) = 0.0 ! to prevent random values reported over areas 
                      ! that are not included as soil points (i.e. ice points)

!----------------------------------------------------------------------
! Set up variables required for LSH scheme:
!----------------------------------------------------------------------
zdepth(:)=0.0

DO n=1,nshyd
   zdepth(n)=zdepth(n-1)+dzsoil(n)
END DO

!-----------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth:
! Initialiase base flow to zero:
!-----------------------------------------------------------------------

DO n=0,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    i=soil_index(j)
    ksz(i,n)=satcon(i,n)
  END DO
END DO

DO n=1,nshyd
!CDIR NODEP
  DO j=1,soil_pts
    qbase_l(soil_index(j),n)=0.0
  END DO
END DO

DO i=1,npnts
  qbase(i)=0.0
  qbase_zw(i)=0.0
  wutot(i)=0.0
  drain(i)=0.0
END DO


IF(l_inland)THEN

 DO i=1,npnts

! Add inland basin outflow to change in soil moisture store

   dsmc_dt(i)=dsmc_dt(i)+inlandout_atm(i)

 END DO
END IF


!-----------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------
! DEPENDS ON: soilmc
CALL soilmc ( npnts,nshyd,soil_pts,soil_index,                    &
              dzsoil,sthu,v_sat,v_wilt,smc )

!-----------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------

DO i=1,npnts
  fch4_wetl(i)=0.0
END DO

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('HYDROL_CABLE',zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol_cable
END MODULE hydrol_cable_mod
