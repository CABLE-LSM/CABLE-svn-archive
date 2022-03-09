MODULE cable_land_sf_explicit_mod

CONTAINS

SUBROUTINE cable_land_sf_explicit(                                             &
    !RETURNED per tile fields as seen by atmosphere
    ftl_surft, fqw_surft, tstar_surft,
    u_s_std_surft, 
! I think these need to be aggreg. over tiles for land_pt
cd_surft, ch_surft,               &
!this appears to be a local variable in jules_land_sf_explicit, formerly in sf_expl. pumped  into fcdch                      
radnet_surft, 
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
    fraca, 
    !Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
    resfs,                                                      &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for
    !snow.
    resft,                                                      &
    ! Roughness lengths for heat and moisture/ momentum  (m)
    fluxes%z0h_surft, fluxes%z0m_surft,           &
!this appears to be a local variable in jules_land_sf_explicit, formerly in sf_expl. pumped  into fcdch           
recip_l_MO_surft, 
    epot_surft, 
!This might as well be local to CABLE
!l_tile_pts, 
    !Surface friction velocity (m/s) (grid box)
    u_s, 
  !CALL cable_explicit_main(                                                    &
      !mype, timestep_number, endstep,!abandoned at 5.7 comparrison as not actually neeeded - useful ONLY 
      !timestep, !this is available via standalone module ONLY. rightly so it is not here as it is not needed
!IN Model fundamentals & dimensions 
      cycleno, numcycles curr_day_number,                                      &
      land_pts, nsurft, npft, sm_levels,                                       & !consistency - half the time we use nsl from CABLE
      latitude, longitude,
      land_index, surft_pts, surft_index,                     &
      !IN Parametrs 
      tile_frac,  !Dont fush this please  - maybe call tile_pts
      !Re-examine consequence of setting to soilt=1      
      psparms%bexp_soilt(:,1,:), 
      psparms%hcon_soilt(:,1,:), 
      psparms%satcon_soilt(:,1,:), 
      psparms%sathh_soilt(:,1,:), 
      psparms%smvcst_soilt(:,1,:), 
      psparms%smvcwt_soilt(:,1,:), 
      psparms%smvccl_soilt(:,1,:), 
      psparms%sthu_soilt(:,1,:), !Unfrozen soil fraction of moisture content !why in parms? 
      psparms%albsoil_soilt,   
      psparms%cosz_ij       !cosine Zenith angle changes with dt? param??
      co2_mmr, !USEd 
      !IN prognostics
      progs%snow_surft, !verify where we need snow depth on explicit call. init only?  
      progs%canopy_surft, !water content retained in canopy ? 
      progs%canht_pft, progs%lai_pft, !canopy height/ LAI                                &
      !IN met forcing
      forcing%lw_down_ij ,                                               &
      fluxes%sw_surft, &    ! Surface net shortwave on tiles (W/m2) 
      forcing%ls_rain, forcing%ls_snow, !large scale precip 
      forcing%qw_1_ij,  !  Total water content (Kg/Kg)  !this appears to be what we need now
      forcing%tl_1_ij, ! Ice/liquid water temperature (k) !this appears to be what we need now
      forcing%pstar_ij, ! Surface pressure (Pascals)
      ainfo%z1_tq_ij, !  height of temperature data
      ainfo%z1_uv_ij, !  height of wind data
      coast%Fland, !this seems to be the one we want but why is it in coastal%        &
      coast%vshr_land_ij !! wind shear surface-to-atm (m per s)!why coastal%
!INOUT: CABLE TYPES containing field data (IN OUT)
    progs_cbl, work_cbl, pars_io_cbl )


! In general CABLE utilizes a required subset of tbe JULES types, however;
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep
USE params_io_mod_cbl,  ONLY: params_io_type      ! and veg/soil parameters

implicit NONE

REAL :: _surft(land_pts,nsurft) 

!RETURNED per tile fields as seen by atmosphere
REAL :: ftl_surft(land_pts,nsurft)    ! Surface FTL for land tiles
REAL :: fqw_surft(land_pts,nsurft)    ! Surface FTL for land tiles
REAL :: tstar_surft(land_pts,nsurft)  ! Surface temperature (K) 
REAL :: u_s_std_surft _surft(land_pts,nsurft) ! Surface friction velocity (m/s) 
                                              ! (standard value) -> aerosols
REAL :: cd_surft(land_pts,nsurft)     ! Drag coefficient
REAL :: ch_surft(land_pts,nsurft)     ! Transfer coeff for heat/moisture
REAL :: radnet_surft(land_pts,nsurft) ! Net SW ( appears to be local )

REAL :: fraca(land_pts,nsurft)        ! Frac of moist. flux with aerodynamic 
                                      ! resistance for !snow-free land tiles
  
!Combined soil, stomatal and aerodynamic resistance factor for fraction
!(1-FRACA) {-!!THIS IS NOT CONSISTENT WITH ABOVE COMMENT FOR fraca } of snow-free land tiles.
REAL :: resfs(land_pts,nsurft) 

!Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for snow
REAL :: resft(land_pts,nsurft) 
  
  
    
REAL :: epot_surft(land_pts,nsurft) 


!CABLE TYPES containing field data
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs
TYPE(work_vars_type), INTENT(OUT)        :: work
TYPE(params_io_type), INTENT(IN)         :: pars

!local Vars


!This might as well be local to CABLE
!l_tile_pts, 

ftl_surft(:,:)  = 0.0
fqw_surft(:,:)  = 0.0
epot_surft(:,:) = 0.0
z0h_surft(:,:) = 0.0
z0m_surft(:,:) = 0.0


DO n = 1,ntype
  DO l = 1, land_pts
    tile_frac(l,n) = frac(l,n)
  END DO
END DO

WRITE(6,*) "Currently CABLE explicit@6.3 is not implemented"
RETURN

END SUBROUTINE cable_land_sf_explicit

END MODULE cable_land_sf_explicit_mod


