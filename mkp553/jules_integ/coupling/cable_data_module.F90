! this module is to live in the JULES src code to be compiled with the UM
! data/fields/vars in this module are required when calling CABLE from UM
! this module serves as a mechanism of transporting all those vars from the
! top_level in the UM to the calling points. It also serves as a library
! defining what vars are required for CABLE and what it can give back.
! It is further intended that this module serve as a template for the future 
! data atructure to be adopted by CABLE

!jhan***************************************************************************
   !this needs to be cleaned up, types renamed appropriately etc.
   !impl_params (eg)are vars pointed to by cable%typed BUT lazily are passed just
   !before the call anyway. it has already been shown that the method works SO 
   !for the sake of getting the model working - put them here and clan up later 
   !by redistributing to higher level calls
!jhan***************************************************************************
   
module cable_data_mod
   !use define_dimensions, only : nrb 
   implicit none
  
   ! remove hard wired values
   integer, parameter :: ms= 6, & ! soil levels must be same as UM
                         msn = 3   


!------------------------------------------------------------------------------

   type model_params

      INTEGER ::                                                      &
         endstep,            & !
         row_length, rows,  & !
         sm_levels,          & !
         land_pts,          & !
         ntiles,          & !
         npft,          & !
         timestep_number!,    & !
         !mype

      !REAL, POINTER ::                                                      &
      INTEGER ::                                                      &
         timestep_width
      
      REAL, DIMENSION(:), POINTER ::                                           &
         dzsoil
      
      REAL, DIMENSION(:,:), POINTER ::                                         &
         latitude, &
         longitude!, &

   end type model_params
  
!------------------------------------------------------------------------------
   ! these are prognostic initializations that are generally not calculated
   ! dynamically per timestep
   type prognostic_params
   
      REAL, DIMENSION(:,:), POINTER ::                                         &
         tile_frac!, &
   
   end type prognostic_params
   
!------------------------------------------------------------------------------

   ! cable prognostic_vars
   type cable_vars

      INTEGER, DIMENSION(:,:), POINTER :: &
         snow_flg3l
      
      REAL, DIMENSION(:,:), POINTER :: &
         snow_rho1l,    & !
         snow_age

      REAL, DIMENSION(:,:,:), POINTER :: &
         tsoil_tile, &
         smcl_tile, &
         sthf_tile, &
         snow_depth3l, &
         snow_mass3l, &
         snow_tmp3l, &
         snow_rho3l, &  
         sthu_tile !jhan: c nwe kill this

      REAL, DIMENSION(:,:,:), POINTER :: &
         snow_cond ! no init from file

   end type cable_vars
 
!------------------------------------------------------------------------------

   ! forcing vars 
   type forcing_vars
      
      real, dimension(:,:,:), pointer :: & 
      ShortWave => NULL()!, & ! => surf_down_sw

      !LongWave, &
      !AirTemper, &
      !SurfPressure, &
      !Humidity, &
      !WindSpeed, &
      !Precip 
   end type forcing_vars

!------------------------------------------------------------------------------
   ! TYPEd vars pased onto cable after UM version being pointed to 
   type UM_params

      INTEGER :: dim_cs1, dim_cs2

      INTEGER, DIMENSION(:), POINTER ::                                        &
         land_index,       &
         tile_pts

      INTEGER, DIMENSION(:,:), POINTER ::                                      &
         tile_index
 
      REAL, DIMENSION(:,:), POINTER :: &
         bexp, & !
         hcon, & !
         satcon, & !
         sathh, & !
         smvcst, & !
         smvcwt, & !
         smvccl!, & !
         !albsoil, &
         !CANOPY_GB
      
      REAL, DIMENSION(:), POINTER :: &
!         bexp, & !
!         hcon, & !
!         satcon, & !
!         sathh, & !
!         smvcst, & !
!         smvcwt, & !
!         smvccl, & !
         albsoil, &
         CANOPY_GB, &
         GS
      
      REAL :: co2_mmr

      REAL, DIMENSION(:,:), POINTER :: &
         sthu, &
         smcl, &
         sthf, &
         tot_alb

      !REAL, DIMENSION(:,:,:), POINTER :: &
      !   land_alb
      
      REAL, DIMENSION(:,:), POINTER :: &
         snow_tile, &
         vshr_land, &
         sin_theta_latitude, &
         pstar, &
         canht_ft, & !
         lai_ft,   & !
         land_alb,   & !
         canopy

      REAL, DIMENSION(:), POINTER :: &
        cos_zenith_angle

      REAL, DIMENSION(:,:), POINTER :: &
        lw_down, &
        ls_rain, &
        ls_snow 

      REAL, DIMENSION(:,:), POINTER :: &
         tl_1, &
         qw_1

       real, dimension(:,:), pointer :: & 
         Z1_TQ, &
         Z1_UV, &
         U_S, &
         conv_rain, & 
         conv_snow
        
      real, dimension(:,:), pointer :: & 
         FTL_TILE, &
         fqw_tile, &
         tstar_tile, &
         U_S_STD_TILE, &
         CD_TILE, &
         CH_TILE, &
         RADNET_TILE, &
         FRACA, &
         rESFS, &
         RESFT, &
         Z0H_TILE, &
         Z0M_TILE, &
         RECIP_L_MO_TILE, &
         EPOT_TILE  
   
      real, dimension(:), POINTER :: & 
         FLAND(:)

    real, dimension(:,:), pointer :: & 
      !(LAND_PTS,NTILES)                                       &
      NPP_FT, &
      GPP_FT,&
      RESP_S_TILE,     &
      RESP_P_FT, &
      G_LEAF
      
   real, dimension(:), pointer :: & 
      !land_pts
      NPP,&
      GPP,&
      RESP_P
      
   real, dimension(:), pointer :: & 
      !(DIM_CS2)                                           &
      RESP_S_TOT
     
   real, dimension(:,:), pointer :: & 
     !(LAND_PTS,DIM_CS1)                                      &
     RESP_S

   Real, dimension(:,:,:), pointer ::                       &
         alb_tile

   Real, dimension(:,:), pointer ::                       &
         land_albedo

   end type UM_params
  
!------------------------------------------------------------------------------
   
   type impl_params 
   
   real, dimension(:,:), pointer :: & 
      !(ROW_LENGTH,ROWS),                                         &  
      dtl_1, &
      dqw_1, &
      FTL_1,&
      FQW_1,  &
      SURF_HT_FLUX_LAND
      
   real, dimension(:,:), pointer :: & 
      !(LAND_PTS,SM_LEVELS)                                       &
      T_SOIL
      
   real, dimension(:,:), pointer :: & 
      !(LAND_PTS,NTILES)                                       &
      ECAN_TILE,&
      ESOIL_TILE,&
      EI_TILE,&
      T1P5M_TILE, &
      Q1P5M_TILE, &
      MELT_TILE

   end type impl_params

!------------------------------------------------------------------------------

   type hyd_type
      real, dimension(:), pointer :: &
         sub_surf_roff, &
         surf_roff, &
         tot_tfall, &
         LYING_SNOW
   end type hyd_type


!------------------------------------------------------------------------------

   type tmp_pvars
   
      LOGICAL, DIMENSION(:,:), POINTER ::                                      &
         L_tile_pts

      Real ::                                                                  &
         Epsilon,                                                              &
         c_virtual,                                                            &
         D_T,                                                                  &
         DS_RATIO,                                                             & 
         LH 
       
      REAL :: rho_water

   end type tmp_pvars

!------------------------------------------------------------------------------

   type model 
      type (UM_params) :: um 
      type (tmp_pvars) :: tmp 
      type (model_params) :: mp
      type (prognostic_params) :: ppar
      type (cable_vars) :: cable 
      type (forcing_vars) :: forcing 
      type (impl_params) :: im 
      type (hyd_type) :: hyd

      !integer, allocatable :: gridcell(:) 
      !real, allocatable :: lat(:), lon(:)
   end type model 

!------------------------------------------------------------------------------
      
   !instantiate types
   !type (model_constants),save, target :: const
   type (model),save :: cable 
   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

CONTAINS


!=============================================================================
! Routine for allocating any explicitly allocated pointers during
! initialisation
!=============================================================================
  SUBROUTINE cable_allocate_ptrs(row_length, rows, land_pts, ntiles, sm_levels)

    USE ereport_mod, ONLY : ereport

    INTEGER, INTENT(IN) ::                                                    &
      row_length, rows, land_pts, ntiles, sm_levels

    INTEGER :: i
    REAL :: smcl_layers(6), tsoil_layers(6)

!-----------------------------------------------------------------------------

    IF( sm_levels > 6 ) THEN
      CALL ereport("cable_allocate_ptrs", 101,                                &
                   "Only sm_levels <= 6 is currently supported with CABLE")
    END IF

    ALLOCATE( cable% tmp% l_tile_pts(land_pts,ntiles))
    cable% tmp% l_tile_pts = .FALSE.

    ALLOCATE(cable% um% sin_theta_latitude(row_length, rows))

    ALLOCATE(cable% um% tile_pts(ntiles))
    ALLOCATE(cable% um% tile_index(land_pts, ntiles))
    ALLOCATE(cable% forcing% ShortWave(row_length, rows, 4))

    ALLOCATE(cable% cable% snow_rho1l(land_pts, ntiles))
    cable% cable% snow_rho1l(:,:) = 140.0

    ALLOCATE(cable% cable% snow_age(land_pts, ntiles))
    cable% cable% snow_age(:,:) = 0.0

    ALLOCATE(cable% cable% snow_flg3l(land_pts, ntiles))
    cable% cable% snow_flg3l(:,:) = 0

    ALLOCATE(cable% cable% snow_rho3l(land_pts, ntiles, 3))
    cable% cable% snow_rho3l(:,:,:) = 140.0

    ALLOCATE(cable% cable% snow_depth3l(land_pts, ntiles, 3))
    cable% cable% snow_depth3l(:,:,:) = 0.0

    ALLOCATE(cable% cable% snow_tmp3l(land_pts, ntiles, 3))
    cable% cable% snow_tmp3l(:,:,:) = 273.1

    ALLOCATE(cable% cable% snow_mass3l(land_pts, ntiles, 3))
    cable% cable% snow_mass3l(:,:,:) = 0.0

    ALLOCATE(cable% cable% snow_cond(land_pts, ntiles, 3))
    cable% cable% snow_cond = -huge(1.0)

    ALLOCATE(cable% cable% smcl_tile(land_pts, ntiles, sm_levels))
    smcl_layers(:) = (/ 7.475, 19.78, 52.72, 141.5, 386.5, 1111.0 /)
    DO i = 1,sm_levels
      cable% cable% smcl_tile(:,:,i) = smcl_layers(i)
    END DO

    ALLOCATE(cable% cable% sthf_tile(land_pts, ntiles, sm_levels))
    cable% cable% sthf_tile(:,:,:) = 0.0

    ALLOCATE(cable% cable% tsoil_tile(land_pts, ntiles, sm_levels))
    tsoil_layers(:) = (/ 277.7, 277.8, 278.0, 278.6, 280.1, 283.6 /)
    DO i = 1,sm_levels
      cable% cable% tsoil_tile(:,:,i) = tsoil_layers(i)
    END DO

    ALLOCATE(cable% cable% sthu_tile(land_pts, ntiles, sm_levels))
    cable% cable% sthu_tile = -huge(1.0)

    cable% tmp% Epsilon = 0.62198
    cable% tmp% c_virtual =  1.0 / cable% tmp% Epsilon - 1.0

  END SUBROUTINE cable_allocate_ptrs


!=============================================================================
! Routine for setting the pointers and variables required by cable_rad_driver,
! except sw_down
!
! NOTE: cable_rad_driver is not called on the first timestep - JULES values
!       are used instead - because the CABLE albedo scheme requires sw_down
!       which is not known until after the call to tile_albedo
!=============================================================================
  SUBROUTINE cable_radiation_setup(                                           &
             timestep_number, cosz, snow_tile, albsoil, land_albedo, alb_tile)

    INTEGER, INTENT(IN) :: timestep_number

    REAL, INTENT(IN), TARGET ::                                               &
      cosz(:), snow_tile(:,:), albsoil(:)

    REAL, INTENT(INOUT), TARGET ::                                            &
      land_albedo(:,:), alb_tile(:,:,:)

!-----------------------------------------------------------------------------

    cable% mp% timestep_number  = timestep_number

    cable% um% cos_zenith_angle => cosz
    cable% um% snow_tile        => snow_tile
    cable% um% albsoil          => albsoil
    cable% um% land_albedo      => land_albedo
    cable% um% alb_tile         => alb_tile

  END SUBROUTINE cable_radiation_setup


!=============================================================================
! Routine for setting up variables and pointers required for the call to
! cable_explicit_driver
!=============================================================================
  SUBROUTINE cable_explicit_setup(                                            &
    row_length, rows, land_pts, ntiles, npft, sm_levels, dim_cs1, dim_cs2,    &
    timestep, timestep_number,                                                &
    land_index, tile_pts, tile_index,                                         &
    co2_mmr, rho_water,                                                       &
    latitude, longitude, dzsoil, bexp, hcon, satcon, sathh, smvcst, smvcwt,   &
    smvccl, albsoil, fland, cos_zenith_angle, sw_down_4band, lw_down,         &
    ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, snow_tile,  &
    tile_frac, canht_ft, lai_ft, canopy, sthu,                                &
    ftl_tile, fqw_tile, tstar_tile, z0h_tile, z0m_tile, cd_tile, ch_tile,     &
    u_s_std_tile, u_s, radnet_tile, resfs, resft, fraca, recip_l_mo_tile,     &
    epot_tile, gs, npp, npp_ft, gpp, gpp_ft, resp_s, resp_s_tot, resp_p,      &
    resp_p_ft, g_leaf )

    INTEGER, INTENT(IN) ::                                                    &
      row_length, rows, land_pts, ntiles, npft, sm_levels, dim_cs1, dim_cs2,  &
      timestep, timestep_number,                                              &
      tile_pts(:), tile_index(:,:)

    INTEGER, INTENT(IN), TARGET ::                                            &
      land_index(:)

    REAL, INTENT(IN) :: co2_mmr, rho_water

    REAL, INTENT(INOUT), TARGET ::                                            &
      latitude(:,:), longitude(:,:), dzsoil(:), bexp(:,:), hcon(:,:),         &
      satcon(:,:), sathh(:,:), smvcst(:,:), smvcwt(:,:), smvccl(:,:),         &
      albsoil(:), fland(:), cos_zenith_angle(:), sw_down_4band(:,:,:),        &
      lw_down(:,:), ls_rain(:,:), ls_snow(:,:), tl_1(:,:), qw_1(:,:),         &
      vshr_land(:,:), pstar(:,:), z1_tq(:,:), z1_uv(:,:), snow_tile(:,:),     &
      tile_frac(:,:), canht_ft(:,:), lai_ft(:,:), canopy(:,:), sthu(:,:)

    REAL, INTENT(INOUT), TARGET ::                                            &
      ftl_tile(:,:), fqw_tile(:,:), tstar_tile(:,:), z0h_tile(:,:),           &
      z0m_tile(:,:), cd_tile(:,:), ch_tile(:,:), u_s_std_tile(:,:), u_s(:,:), &
      radnet_tile(:,:), resfs(:,:), resft(:,:), fraca(:,:),                   &
      recip_l_mo_tile(:,:), epot_tile(:,:), gs(:), npp(:), npp_ft(:,:),       &
      gpp(:), gpp_ft(:,:), resp_s(:,:), resp_s_tot(:), resp_p(:),             &
      resp_p_ft(:,:), g_leaf(:,:)

!-----------------------------------------------------------------------------

    cable% mp% row_length      = row_length
    cable% mp% rows            = rows
    cable% mp% land_pts        = land_pts
    cable% mp% ntiles          = ntiles
    cable% mp% npft            = npft
    cable% mp% sm_levels       = sm_levels

    cable% um% dim_cs1 = dim_cs1
    cable% um% dim_cs2 = dim_cs2

    cable% mp% timestep_width  = timestep
    cable% mp% timestep_number = timestep_number

    cable% tmp% rho_water = rho_water
    cable% um% co2_mmr    = co2_mmr

    cable% mp% dzsoil    => dzsoil
    cable% mp% latitude  => latitude
    cable% mp% longitude => longitude

    cable% um% sin_theta_latitude = SIN( latitude )

    cable% um% tile_pts       = tile_pts
    cable% um% tile_index     = tile_index
    cable% forcing% ShortWave = sw_down_4band

    cable% um% land_index       => land_index
    cable% ppar% tile_frac      => tile_frac
    cable% um% bexp             => bexp
    cable% um% hcon             => hcon
    cable% um% satcon           => satcon
    cable% um% sathh            => sathh
    cable% um% smvcst           => smvcst
    cable% um% smvcwt           => smvcwt
    cable% um% smvccl           => smvccl
    cable% um% albsoil          => albsoil
    cable% um% snow_tile        => snow_tile
    cable% um% lw_down          => lw_down
    cable% um% cos_zenith_angle => cos_zenith_angle
    cable% um% ls_rain          => ls_rain
    cable% um% ls_snow          => ls_snow
    cable% um% tl_1             => tl_1
    cable% um% qw_1             => qw_1
    cable% um% vshr_land        => vshr_land
    cable% um% pstar            => pstar
    cable% um% z1_tq            => z1_tq
    cable% um% z1_uv            => z1_uv
    cable% um% canopy           => canopy
    cable% um% fland            => fland
    cable% um% sthu             => sthu
    cable% um% canht_ft         => canht_ft
    cable% um% lai_ft           => lai_ft
    cable% um% ftl_tile         => ftl_tile
    cable% um% fqw_tile         => fqw_tile
    cable% um% tstar_tile       => tstar_tile
    cable% um% u_s              => u_s
    cable% um% u_s_std_tile     => u_s_std_tile
    cable% um% cd_tile          => cd_tile
    cable% um% ch_tile          => ch_tile
    cable% um% radnet_tile      => radnet_tile
    cable% um% fraca            => fraca
    cable% um% resfs            => resfs
    cable% um% resft            => resft
    cable% um% z0h_tile         => z0h_tile
    cable% um% z0m_tile         => z0m_tile
    cable% um% recip_l_mo_tile  => recip_l_mo_tile
    cable% um% epot_tile        => epot_tile
    cable% um% gs               => gs
    cable% um% npp              => npp
    cable% um% npp_ft           => npp_ft
    cable% um% gpp              => gpp
    cable% um% gpp_ft           => gpp_ft
    cable% um% resp_s           => resp_s
    cable% um% resp_s_tot       => resp_s_tot
    cable% um% resp_p           => resp_p
    cable% um% resp_p_ft        => resp_p_ft
    cable% um% g_leaf           => g_leaf

  END SUBROUTINE cable_explicit_setup


!=============================================================================
! Routine for setting up variables and pointers required for the call to
! cable_implicit_driver
!=============================================================================
  SUBROUTINE cable_implicit_setup(                                            &
    timestep_width,                                                           &
    ls_rain, conv_rain, ls_snow, conv_snow, dtl_1, dqw_1, ftl_1, ftl_tile,    &
    fqw_1, fqw_tile, tstar_tile, surf_ht_flux_land, ecan_tile, esoil_tile,    &
    ei_tile, radnet_tile, t1p5m_tile, q1p5m_tile, fland, melt_tile )

    INTEGER, INTENT(IN) :: timestep_width

    REAL, INTENT(IN), TARGET ::                                               &
      dtl_1(:,:), dqw_1(:,:), fland(:)

    REAL, INTENT(INOUT), TARGET ::                                            &
      ls_rain(:,:), conv_rain(:,:), ls_snow(:,:), conv_snow(:,:), ftl_1(:,:), &
      ftl_tile(:,:), fqw_1(:,:), fqw_tile(:,:), tstar_tile(:,:),              &
      surf_ht_flux_land(:,:), ecan_tile(:,:), esoil_tile(:,:), ei_tile(:,:),  &
      radnet_tile(:,:), t1p5m_tile(:,:), q1p5m_tile(:,:), melt_tile(:,:)

!-----------------------------------------------------------------------------

    cable% mp% timestep_width = timestep_width

    cable% um% ls_rain           => ls_rain
    cable% um% conv_rain         => conv_rain
    cable% um% ls_snow           => ls_snow
    cable% um% conv_snow         => conv_snow
    cable% im% dtl_1             => dtl_1
    cable% im% dqw_1             => dqw_1
    cable% im% ftl_1             => ftl_1
    cable% um% ftl_tile          => ftl_tile
    cable% im% fqw_1             => fqw_1
    cable% um% fqw_tile          => fqw_tile
    cable% um% tstar_tile        => tstar_tile
    cable% im% surf_ht_flux_land => surf_ht_flux_land
    cable% im% ecan_tile         => ecan_tile
    cable% im% esoil_tile        => esoil_tile
    cable% im% ei_tile           => ei_tile
    cable% um% radnet_tile       => radnet_tile
    cable% im% t1p5m_tile        => t1p5m_tile
    cable% im% q1p5m_tile        => q1p5m_tile
    cable% um% fland             => fland
    cable% im% melt_tile         => melt_tile

  END SUBROUTINE cable_implicit_setup


!=============================================================================
! Routine for setting up variables and pointers required for the call to
! cable_extra_driver
!=============================================================================
  SUBROUTINE cable_extra_setup( smvcst, t_soil, smcl, sthf, sthu, snow_tile,  &
                                canopy, canopy_gb, lying_snow, surf_roff,     &
                                sub_surf_roff, tot_tfall )

    REAL, INTENT(INOUT), TARGET ::                                            &
      smvcst(:,:), t_soil(:,:), smcl(:,:), sthf(:,:), sthu(:,:),              &
      snow_tile(:,:), canopy(:,:), canopy_gb(:), lying_snow(:), surf_roff(:), &
      sub_surf_roff(:), tot_tfall(:)

!------------------------------------------------------------------------------

    cable% um% smvcst         => smvcst
    cable% im% t_soil         => t_soil
    cable% um% smcl           => smcl
    cable% um% sthf           => sthf
    cable% um% sthu           => sthu
    cable% um% snow_tile      => snow_tile
    cable% um% canopy         => canopy
    cable% um% canopy_gb      => canopy_gb
    cable% hyd% lying_snow    => lying_snow
    cable% hyd% surf_roff     => surf_roff
    cable% hyd% sub_surf_roff => sub_surf_roff
    cable% hyd% tot_tfall     => tot_tfall

  END SUBROUTINE cable_extra_setup

END MODULE cable_data_mod
