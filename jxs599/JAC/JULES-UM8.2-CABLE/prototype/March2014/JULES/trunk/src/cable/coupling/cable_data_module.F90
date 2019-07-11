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
  
   ! public variables. ALL vars above "contains" are deliberately public
    
   ! public subroutines 
   public :: cable_control

   ! remove hard wired values
   integer, parameter :: ms= 6, & ! soil levels must be same as UM
                         msn = 3   

   !! introduced tiled prognostics that have to go through UM's I/O 
   !real, dimension(:), pointer, save ::           &
   !   TSOIL_TILE,    & !
   !   SMCL_TILE,     & !
   !   STHF_TILE,     & !
   !   SNOW_DEPTH3L,  & !
   !   SNOW_MASS3L,   & !
   !   SNOW_TMP3L,    & !
   !   SNOW_RHO3L,    & !
   !   SNOW_RHO1L,    & !
   !   SNOW_AGE,      & !
   !   SNOW_FLG3L       !
   ! introduced tiled prognostics that have to go through JULES's I/O 
   real, dimension(:,:,:), pointer, save ::           &
      SOIL_Temp_CABLE, & !
      SMCL_CABLE,     & !
      STHF_CABLE        !

   real, dimension(:,:,:), pointer, save ::           &
      SNOW_DEPTH_CABLE,  & !
      SNOW_MASS_CABLE,   & !
      SNOW_TEMP_CABLE,    & !
      SNOW_RHO_CABLE       !

   real, dimension(:,:), pointer, save ::           &
      SNOW_RHO1L_CABLE,    & !
      SNOW_AGE_CABLE         !

   real, dimension(:,:), pointer, save ::           &
      SNOW_FLG3L_CABLE       !

!------------------------------------------------------------------------------

   type model_params

      INTEGER, POINTER ::                                                      &
         endstep,            & !       
         row_length, rows,  & !
         sm_levels,          & !
         land_pts,          & !
         ntiles,          & !
         npft,          & !
         timestep_number!,    & !
         !mype

      !REAL, POINTER ::                                                      &
      integer, POINTER ::                                                      &
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

      real, DIMENSION(:,:), POINTER :: &
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
      ShortWave!, & ! => surf_down_sw

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

      LOGICAL, POINTER ::                                                      &
         L_cable

      INTEGER, POINTER ::                                                      &
         dim_cs1, dim_cs2     !

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
      
      REAL, POINTER ::                                              &
         co2_mmr

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
        !cos_zenith_angle, &
        ls_rain, &
        ls_snow 

      REAL, DIMENSION(:,:), POINTER :: &
         tl_1, &
         qw_1

       LOGICAL, dimension(:,:), pointer :: & 
         LAND_MASK

       real, dimension(:,:), pointer :: & 
         SW_DOWN, &  
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
         alb_tile,      & !
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
      real, dimension(:), pointer :: &                                                               &
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
       
      real, POINTER :: & 
         rho_water

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
   

contains


!===============================================================================

!subroutine set_endstep_umodel(fendstep)
!   integer, target :: fendstep
!
!      cable% mp% endstep            => fendstep
!
!end subroutine set_endstep_umodel

!===============================================================================

 ! call cable_control( L_cable, a_step, &!mype, 
 !        timestep_len, row_length,                                             &
 !        rows, land_pts, ntiles, sm_levels, dim_cs1, dim_cs2,                  & 
 !        !sin_theta_latitude, cos_theta_longitude,                              &
 !        land_index, b, hcon, satcon, SATHH, smvcst, smvcwt,                   &
 !        smvccl, albsoil, lw_down, cosz, ls_rain, ls_snow, pstar,              &
 !        CO2_MMR, sthu, smcl, sthf, GS, canopy, land_albedo )


SUBROUTINE cable_control( L_cable, a_step, timestep_len, row_length,     &
             rows, land_pts, ntiles, sm_levels, dim_cs1, dim_cs2,              &
             latitude, longitude,                                              &
             land_index, b, hcon, satcon, sathh, smvcst, smvcwt, smvccl,       &
             albsoil, lw_down, cosz, ls_rain, ls_snow, pstar, CO2_MMR,         &
             sthu, smcl, sthf, GS, canopy_gb , land_albedo )

   LOGICAL, target :: L_CABLE
   
   !fudged at present as NA in JULES
   INTEGER, target :: endstep
                                                &
   INTEGER, target ::                                              &
      !mype,             & !
      a_step,  & !
      !endstep,          & !
      row_length, rows, &
      land_pts,      &
      ntiles,        &
      sm_levels,        &
      dim_cs1, dim_cs2

   integer, target::                                              &
      timestep_len

   REAL, target::                                              &
      !timestep_len, &
      co2_mmr
   
   INTEGER, DIMENSION(:), target::                                              &
      land_index
   
   REAL, DIMENSION(:,:), target ::                                         &
      pstar, &
      land_albedo 
 
   !REAL, DIMENSION(:,:,:), target ::                                         &
   !   land_albedo 
 
   REAL, DIMENSION(:), TARGET :: &
     albsoil, &
     cosz,    &
     canopy_gb, &
     GS 

   !REAL, DIMENSION(:), TARGET :: &
   REAL, DIMENSION(:,:), TARGET :: &
      !albsoil, &
      b, &    !
      hcon, &    !
      satcon, &
      sathh, &
      smvcst, &
      smvcwt, &
      smvccl!, &
      !canopy_gb
 
   REAL, DIMENSION(:,:), TARGET:: &
     lw_down, &
     !cosz, &
     ls_rain, &
     ls_snow, & 
     sthu, &
     smcl, &
     sthf
                           
   REAL, DIMENSION(:,:), TARGET :: &
      latitude, longitude

   !---------------------------------------------------------------------------
   !local vars
   
   LOGICAL :: first_atmstep_call
   
   integer :: i,j, k,n
   
   !---------------------------------------------------------------------------

   if( a_step == 1) first_atmstep_call = .TRUE. 
   
   if( first_atmstep_call ) then 
      !allocate( latitude(row_length,rows) )
      allocate( cable% um% sin_theta_latitude(row_length,rows) )
      !allocate( longitude(row_length,rows) )
      allocate( cable% cable% SNOW_COND(land_pts,NTILES,3))
      allocate( cable% cable% STHU_TILE(land_pts,NTILES,sm_levels))
      allocate( cable% tmp% L_TILE_PTS(land_pts,NTILES))
      !this can be deleted once rm from cable_explicit_driver (call/recieve )
      allocate( cable% um% SW_DOWN(row_length,rows)           )
      allocate( cable% forcing% ShortWave(row_length,rows,4)    )
      allocate( cable% um% TILE_INDEX(land_pts,NTiles) ) 
      allocate( cable% um% TILE_PTS(NTiles) )  
      allocate( cable% um% TOT_ALB(land_pts,ntiles)        )
   endif

   !call cable_parse_isnow( land_pts, ntiles, snow_flg3l,                       &
   !        tsoil_tile, smcl_tile, sthf_tile,                                   &
   !        snow_depth3l, snow_mass3l, snow_tmp3l,                              & 
   !        snow_rho3l, snow_rho1l, snow_age ) 
  
      cable% um% L_cable            => L_cable
      !cable% mp% mype               => mype
      cable% mp% endstep            => endstep     
      cable% mp% timestep_number    => a_step
      cable% mp% timestep_width     => timestep_len
      cable% mp% row_length         => row_length 
      cable% mp% rows               => rows        
      cable% mp% land_pts           => land_pts
      cable% mp% ntiles             => ntiles
      cable% mp% sm_levels          => sm_levels
      
      cable% um% dim_cs1           => dim_cs1
      cable% um% dim_cs2           => dim_cs2

      cable% cable% tsoil_tile       => soil_temp_CABLE
      cable% cable% smcl_tile        => smcl_CABLE
      cable% cable% sthf_tile        => sthf_CABLE
      cable% cable% snow_depth3l     => snow_depth_CABLE
      cable% cable% snow_mass3l      => snow_mass_CABLE
      cable% cable% snow_tmp3l       => snow_temp_CABLE
      cable% cable% snow_rho3l       => snow_rho_CABLE
      cable% cable% snow_flg3l       => snow_flg3l_CABLE
      cable% cable% snow_rho1l       => snow_rho1l_CABLE
      cable% cable% snow_age         => snow_age_CABLE

      !jhan: re-implement sin_theta_lat by computing here      
      !latitude  = asin( sin_theta_latitude )
      !longitude = acos( cos_theta_longitude )
      cable% mp% latitude           => latitude
      cable% mp% longitude          => longitude


      cable% um% sin_theta_latitude = sin( cable% mp% latitude )

      cable% um% land_index         => land_index 

      cable% um% sthu             => sthu        
      
      cable% um% sthf             => sthf        
      cable% um% smcl             => smcl        
      
      cable% um% land_alb         => land_albedo   
          

      cable% um% bexp     => b 
      cable% um% hcon     => hcon
      cable% um% satcon   => satcon
      cable% um% sathh    => sathh
      cable% um% smvcst   => smvcst  
      cable% um% smvcwt   => smvcwt  
      cable% um% smvccl   => smvccl  
      cable% um% albsoil  => albsoil 
      
      cable% um% pstar => pstar 

      cable% um% lw_down   => lw_down
      cable% um% cos_zenith_angle   => cosz
      cable% um% ls_rain     => ls_rain 
      cable% um% ls_snow      => ls_snow
      cable% um% co2_mmr      => co2_mmr

      cable% um% gs => gs
      cable% um% canopy_gb => canopy_gb

      cable% cable% snow_cond = -huge(1.)
      cable% cable% sthu_tile = -huge(1.)
      cable% tmp% l_tile_pts = .false.
      cable% tmp% Epsilon = 0.62198 
      cable% tmp% c_virtual =  1. / cable% tmp% Epsilon - 1. 

END SUBROUTINE cable_control
 
!===============================================================================
 
SUBROUTINE cable_control2( npft, tile_frac, snow_tile, vshr_land, canopy,      &
              canht_ft, lai_ft, conv_rain, conv_snow, NPP,NPP_FT,              &
              GPP, GPP_FT, RESP_S, rESP_S_TOT, RESP_S_TILE, RESP_P,            &
              RESP_P_FT, G_LEAF, Radnet_TILE, Lying_snow, surf_roff,           &
              sub_surf_roff, tot_tfall )

   INTEGER, target ::                                              &
      npft 

   REAL, DIMENSION(:,:), TARGET:: &
      tile_frac, &
      snow_tile, &
      vshr_land
 
   REAL, DIMENSION(:,:), TARGET ::                                         &
     canopy, &
     canht_ft, &
     lai_ft 
 
   REAL, DIMENSION(:,:), TARGET ::                                         &
     !(row_length, rows)                                     &
     conv_rain, &
     conv_snow
   
   real, dimension(:), target :: & 
      GPP, & ! Gross primary productivity (kg C/m2/s).
      NPP, & ! Net primary productivity
      RESP_P ! Plant respiration (kg C/m2/s).
   
   real, dimension(:,:), target :: & 
      GPP_FT, & !     on PFTs (kg C/m2/s).
      NPP_FT, & ! Net primary productivity (kg C/m2/s).
      G_LEAF, & ! Leaf turnover rate (/360days).
      RESP_P_FT, & !  Plant respiration on PFTs (kg C/m2/s).
      RESP_S_TILE  ! Soil respiration on tiles (kg C/m2/s).
              
   real, dimension(:,:), target :: & 
      RESP_S ! Soil respiration (kg C/m2/s).
   
   real, dimension(:), target :: & 
      RESP_S_TOT ! OUT total RESP_S over pools

   real, dimension(:,:), target :: & 
      RADNET_TILE
 
   real, dimension(:), target :: &                                                               &
      sub_surf_roff, &
      surf_roff, &
      tot_tfall, &
      LYING_SNOW
   
      cable% mp% npft            => npft
      cable% ppar% tile_frac     => tile_frac
      cable% um% snow_tile       => snow_tile
      cable% um% canopy => canopy
      cable% um% canht_ft  => canht_ft
      cable% um% lai_ft  => lai_ft
      cable% um% vshr_land       => vshr_land 
      cable% um% conv_rain       => conv_rain
      cable% um% conv_snow       => conv_snow
      cable% um% NPP => NPP
      cable% um% NPP_FT => NPP_FT
      cable% um% GPP => GPP
      cable% um% GPP_FT => GPP_FT
      cable% um% RESP_S => RESP_S
      cable% um% RESP_S_TOT => RESP_S_TOT
      cable% um% RESP_S_TILE => RESP_S_TILE
      cable% um% RESP_P => RESP_P
      cable% um% RESP_P_FT => RESP_P_FT
      cable% um% G_LEAF => G_LEAF
      cable% um% RADNET_TILE => RADNET_TILE
      cable% hyd% sub_surf_roff  => sub_surf_roff
      cable% hyd% surf_roff      => surf_roff
      cable% hyd% tot_tfall      => tot_tfall
      cable% hyd% LYING_SNOW     => LYING_SNOW

END SUBROUTINE cable_control2


!===============================================================================


!SUBROUTINE cable_bdy_layr( TL, qw )  
SUBROUTINE cable_control3( TL, qw )  

      REAL, DIMENSION(:,:), TARGET:: &
         tl, &
         qw

      cable% um% tl_1 => tl 
      cable% um% qw_1 => qw 

END SUBROUTINE cable_control3


!===============================================================================


SUBROUTINE cable_control5( alb_tile, land_albedo,         &
                  TILE_PTS, TILE_INDEX )        

   INTEGER, DIMENSION(:) ::                                        &
      tile_pts

   INTEGER, DIMENSION(:,:) ::                                      &
      tile_index

   Real, dimension(:,:,:), target ::                       &
      alb_tile
   
   Real, dimension(cable% mp% rows, cable% mp% row_length,4), target ::                       &
      land_albedo


   cable% um% alb_tile => alb_tile
   cable% um% land_albedo => land_albedo
   cable% um% TILE_PTS = TILE_PTS
   cable% um% TILE_INDEX = TILE_INDEX

END SUBROUTINE cable_control5 


!===============================================================================

!jhan: this is a very temp HACK - for offline SW is split in cable radiation
!module by spitter. online it recieves th SW calculated by the UM rad scheme 
SUBROUTINE cable_control4( sw_down )
   
   Real, dimension(:,:) :: sw_down
   !Real, dimension( cable% mp% row_length, cable% mp% rows, 4) :: surf_down_sw

     !jhan: offline receives total SW and splits (CABLE uses subr spitter)
     cable% forcing% ShortWave(:,:,1)    = sw_down(:,:) 
     cable% forcing% ShortWave(:,:,2)    = 0. 
     cable% forcing% ShortWave(:,:,3)    = 0. 
     cable% forcing% ShortWave(:,:,4)    = 0. 
     !cable% forcing% ShortWave(:,:,1)    = sw_down(:,:) / 4.
     !cable% forcing% ShortWave(:,:,2)    = sw_down(:,:) / 4.
     !cable% forcing% ShortWave(:,:,3)    = sw_down(:,:) / 4.
     !cable% forcing% ShortWave(:,:,4)    = sw_down(:,:) / 4.
     

END SUBROUTINE cable_control4



!===============================================================================

SUBROUTINE cable_control6( z1_tq, z1_uv, Fland, dzsoil, LAND_MASK, FTL_TILE, &
             FQW_TILE, TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE, FRACA, &
             rESFS, RESFT, Z0H_TILE, Z0M_TILE, RECIP_L_MO_TILE, EPOT_TILE )

   real, dimension(:,:), target :: & 
      Z1_TQ, &
      Z1_UV, &
      U_S 
     
   real, dimension(:,:), target :: & 
      !FTL_TILE(land_pts,NTILES)                                        &
      FTL_TILE, &
      fqw_tile, &
      tstar_tile, &
      U_S_STD_TILE, &
      CD_TILE, &
      CH_TILE, &
      FRACA, &
      rESFS, &
      RESFT, &
      Z0H_TILE, &
      Z0M_TILE, &
      RECIP_L_MO_TILE, &
      EPOT_TILE  

   real, target :: & 
      rho_water

   real, dimension(:), target :: & 
      !FLAND(land_pts)
      FLAND(:)

   logical, dimension(:,:), target :: & 
      LAND_MASK

   real, dimension(ms), target :: dzsoil
   !real, target :: dzsoil
 
      cable% um% Z1_TQ => Z1_TQ
      cable% um% Z1_UV => Z1_UV
      cable% um% U_S => U_S 
      
     
      cable% um% FTL_TILE => FTL_TILE
      cable% um% fqw_tile => fqw_tile
      cable% um% tstar_tile => tstar_tile
      cable% um% U_S_STD_TILE => U_S_STD_TILE
      cable% um% CD_TILE => CD_TILE
      cable% um% CH_TILE => CH_TILE
      cable% um% FRACA => FRACA
      cable% um% rESFS => rESFS
      cable% um% RESFT => RESFT
      cable% um% Z0H_TILE => Z0H_TILE
      cable% um% Z0M_TILE => Z0M_TILE
      cable% um% RECIP_L_MO_TILE => RECIP_L_MO_TILE
      cable% um% EPOT_TILE => EPOT_TILE  
      
      cable% tmp% rho_water => rho_water
      
      cable% um% FLAND => FLAND
     
      cable% um% LAND_MASK => LAND_MASK
     
      cable% mp% dzsoil => dzsoil
       

END SUBROUTINE cable_control6


!===============================================================================


SUBROUTINE cable_control7(                      &
                     dtl_1, &
                     dqw_1, &
                     T_SOIL, &
                      FTL_1,&
                      FQW_1,  &
                     SURF_HT_FLUX_LAND, &
                     ECAN_TILE,&
                     ESOIL_TILE,&
                     EI_TILE,&
                     T1P5M_TILE, &
                     Q1P5M_TILE, &
                     MELT_TILE &
                  )
   
   real, dimension(:,:), target :: & 
      !(ROW_LENGTH,ROWS),                                         &  
      dtl_1, &
      dqw_1, &
      FTL_1,&
      FQW_1,  &
      SURF_HT_FLUX_LAND
      
   real, dimension(:,:), target :: & 
      !(LAND_PTS,SM_LEVELS)                                       &
      T_SOIL
      
   real, dimension(:,:), target :: & 
      !(LAND_PTS,NTILES)                                       &
      ECAN_TILE,&
      ESOIL_TILE,&
      EI_TILE,&
      T1P5M_TILE, &
      Q1P5M_TILE, &
      MELT_TILE
      
      cable% im% dtl_1 => dtl_1
      cable% im% dqw_1 => dqw_1
      cable% im% T_SOIL => T_SOIL
      cable% im% FTL_1 => FTL_1
      cable% im% FQW_1 => FQW_1
      cable% im% SURF_HT_FLUX_LAND => SURF_HT_FLUX_LAND
      cable% im% ECAN_TILE => ECAN_TILE
      cable% im% ESOIL_TILE => ESOIL_TILE
      cable% im% EI_TILE => EI_TILE
      cable% im% T1P5M_TILE => T1P5M_TILE
      cable% im% Q1P5M_TILE => Q1P5M_TILE
      cable% im% MELT_TILE => MELT_TILE
 

END SUBROUTINE cable_control7

!===============================================================================
!
!   subroutine cable_point_isnow(isnow_flg3l, &
!                          ftsoil_tile, fsmcl_tile, fsthf_tile,     & !
!                          fsnow_depth3l, fsnow_mass3l, fsnow_tmp3l,    & !
!                          fsnow_rho3l, fsnow_rho1l, fsnow_age )
!    
!      integer, DIMENSION(:,:), target :: isnow_flg3L
!      
!      REAL, DIMENSION(:,:), target :: &
!         fsnow_rho1l,    & !
!         fsnow_age
!
!      REAL, DIMENSION(:,:,:), target :: &
!         ftsoil_tile, &
!         fsmcl_tile, &
!         fsthf_tile, &
!         fsnow_depth3l, &
!         fsnow_mass3l, &
!         fsnow_tmp3l, &
!         fsnow_rho3l
!
!      cable% cable% snow_flg3l      => isnow_flg3l
!
!      cable% cable% tsoil_tile       => ftsoil_tile
!      cable% cable% smcl_tile        => fsmcl_tile
!      cable% cable% sthf_tile        => fsthf_tile
!      cable% cable% snow_depth3l     => fsnow_depth3l
!      cable% cable% snow_mass3l      => fsnow_mass3l
!      cable% cable% snow_tmp3l       => fsnow_tmp3l
!      cable% cable% snow_rho3l       => fsnow_rho3l
!      cable% cable% snow_rho1l       => fsnow_rho1l
!      cable% cable% snow_age         => fsnow_age
!
!   end subroutine cable_point_isnow
!
!!===============================================================================
!! curiously does not work with locally declared j pointers
!subroutine cable_set_atm_pointers( SI, NITEMS,NSECTS,N_INTERNAL_MODEL, &
!                                   Sect_No,im_index, & 
!                                   jTSOIL_TILE, jSMCL_TILE, jSTHF_TILE,            &
!                                   jSNOW_DEPTH3L, jSNOW_MASS3L, jSNOW_TMP3L,       &
!                                   jSNOW_RHO3L, jSNOW_RHO1L, jSNOW_AGE, jSNOW_FLG3L )
!
!   implicit none 
!   ! Address of item from generating plug compatible routine (often
!   ! workspace) !UM: include/declaration/typsts.h
!   INTEGER :: NITEMS,NSECTS,N_INTERNAL_MODEL
!   INTEGER :: SI(  NITEMS,0:NSECTS,N_INTERNAL_MODEL)
!   integer :: Sect_No, im_index
!
!   INTEGER :: JTSOIL_TILE(ms)  ! Tiled soil temperature
!   INTEGER :: JSMCL_TILE(ms)   ! Tiled soil moisture content in layers
!   INTEGER :: JSTHF_TILE(ms)   ! Tiled frozen soil moisture fraction
!   INTEGER :: JSNOW_DEPTH3L(msn)        ! Tiled snow depth
!   INTEGER :: JSNOW_MASS3L(msn)         ! Tiled snow mass
!   INTEGER :: JSNOW_TMP3L(msn)          ! Tiled snow temperature
!   INTEGER :: JSNOW_RHO3L(msn)          ! Tiled snow density
!   INTEGER :: JSNOW_RHO1L             ! Tiled mean snow density
!   INTEGER :: JSNOW_AGE               ! Tiled snow age
!   INTEGER :: JSNOW_FLG3L             ! Flag for use of 3 level snow scheme
!
!   !--- allow for 6 layers
!   JTSOIL_TILE(1) = SI(301,Sect_No,im_index)
!   JTSOIL_TILE(2) = SI(302,Sect_No,im_index)
!   JTSOIL_TILE(3) = SI(303,Sect_No,im_index)
!   JTSOIL_TILE(4) = SI(304,Sect_No,im_index)
!   JTSOIL_TILE(5) = SI(305,Sect_No,im_index)
!   JTSOIL_TILE(6) = SI(306,Sect_No,im_index)
!
!   JSMCL_TILE(1) = SI(307,Sect_No,im_index)
!   JSMCL_TILE(2) = SI(308,Sect_No,im_index)
!   JSMCL_TILE(3) = SI(309,Sect_No,im_index)
!   JSMCL_TILE(4) = SI(310,Sect_No,im_index)
!   JSMCL_TILE(5) = SI(311,Sect_No,im_index)
!   JSMCL_TILE(6) = SI(312,Sect_No,im_index)
!   
!   JSTHF_TILE(1) = SI(313,Sect_No,im_index)
!   JSTHF_TILE(2) = SI(314,Sect_No,im_index)
!   JSTHF_TILE(3) = SI(315,Sect_No,im_index)
!   JSTHF_TILE(4) = SI(316,Sect_No,im_index)
!   JSTHF_TILE(5) = SI(317,Sect_No,im_index)
!   JSTHF_TILE(6) = SI(318,Sect_No,im_index)
!
!   JSNOW_TMP3L(1) = SI(323,Sect_No,im_index)
!   JSNOW_TMP3L(2) = SI(324,Sect_No,im_index)
!   JSNOW_TMP3L(3) = SI(325,Sect_No,im_index)
!   JSNOW_RHO3L(1) = SI(326,Sect_No,im_index)
!   JSNOW_RHO3L(2) = SI(327,Sect_No,im_index)
!   JSNOW_RHO3L(3) = SI(328,Sect_No,im_index)
!   JSNOW_RHO1L = SI(329,Sect_No,im_index)
!   JSNOW_AGE = SI(330,Sect_No,im_index)
!   JSNOW_FLG3l = SI(331,Sect_No,im_index)
!
!   JSNOW_DEPTH3L(1) = SI(332,Sect_No,im_index)
!   JSNOW_DEPTH3L(2) = SI(333,Sect_No,im_index)
!   JSNOW_DEPTH3L(3) = SI(334,Sect_No,im_index)
!   JSNOW_MASS3L(1) = SI(335,Sect_No,im_index)
!   JSNOW_MASS3L(2) = SI(336,Sect_No,im_index)
!   JSNOW_MASS3L(3) = SI(337,Sect_No,im_index)
!
!end subroutine cable_set_atm_pointers
!
!!===============================================================================
!
!
!
!subroutine cable_set_atm_fields( D1, LEN_TOT, land_pts,no_halo,sm_levels,ntiles, &
!                                 jTSOIL_TILE, jSMCL_TILE, jSTHF_TILE,            &
!                                 jSNOW_DEPTH3L, jSNOW_MASS3L, jSNOW_TMP3L,       &
!                                 jSNOW_RHO3L, jSNOW_RHO1L, jSNOW_AGE, jSNOW_FLG3L )
!
!   use field_length_mod , only : field_length
!   implicit none 
!
!   REAL,    TARGET, INTENT(IN) :: D1(LEN_TOT)
!   integer, intent(in) :: LEN_TOT, land_pts,no_halo,sm_levels,ntiles
!      INTEGER :: jTSOIL_TILE(ms)  ! Tiled soil temperature
!      INTEGER :: jSMCL_TILE(ms)   ! Tiled soil moisture content in layers
!      INTEGER :: jSTHF_TILE(ms)   ! Tiled frozen soil moisture fraction
!      INTEGER :: jSNOW_DEPTH3L(msn)        ! Tiled snow depth
!      INTEGER :: jSNOW_MASS3L(msn)         ! Tiled snow mass
!      INTEGER :: jSNOW_TMP3L(msn)          ! Tiled snow temperature
!      INTEGER :: jSNOW_RHO3L(msn)          ! Tiled snow density
!      INTEGER :: jSNOW_RHO1L             ! Tiled mean snow density
!      INTEGER :: jSNOW_AGE               ! Tiled snow age
!      INTEGER :: jSNOW_FLG3L             ! Flag for use of 3 level snow scheme
!
!     TSOIL_TILE  => D1(jTSOIL_TILE(1):jTSOIL_TILE(1)+  &
!                         field_length(land_pts,no_halo,sm_levels*ntiles))
!     SMCL_TILE   => D1(jSMCL_TILE(1):jSMCL_TILE(1)+  &
!                          field_length(land_pts,no_halo,sm_levels*ntiles))
!     STHF_TILE   => D1(jSTHF_TILE(1):jSTHF_TILE(1)+  &
!                   field_length(land_pts,no_halo,sm_levels*ntiles))
!      ! MRD - should be a parameter for number of snow levels here rather than 3
!     SNOW_DEPTH3L=> D1(jSNOW_DEPTH3L(1):jSNOW_DEPTH3L(1)+  &
!                          field_length(land_pts,no_halo,3*ntiles))
!     SNOW_MASS3L => D1(jSNOW_MASS3L(1):jSNOW_MASS3L(1)+  &
!                          field_length(land_pts,no_halo,3*ntiles))
!     SNOW_TMP3L  => D1(jSNOW_TMP3L(1):jSNOW_TMP3L(1)+  &
!                          field_length(land_pts,no_halo,3*ntiles))
!     SNOW_RHO3L  => D1(jSNOW_RHO3L(1):jSNOW_RHO3L(1)+  &
!                          field_length(land_pts,no_halo,3*ntiles))
!     SNOW_RHO1L  => D1(jSNOW_RHO1L:jSNOW_RHO1L+  &
!                          field_length(land_pts,no_halo,ntiles))
!     SNOW_AGE    => D1(jSNOW_AGE:jSNOW_AGE+  &
!                          field_length(land_pts,no_halo,ntiles))
!     SNOW_FLG3L  => D1(jSNOW_FLG3L:jSNOW_FLG3L+  &
!                          field_length(land_pts,no_halo,ntiles))
!
!end subroutine cable_set_atm_fields
!

!===============================================================================


end module cable_data_mod 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!subroutine cable_parse_isnow(land_pts, ntiles, isnow_flg3l, &
!                          tsoil_tile, smcl_tile, sthf_tile,     & !
!                          snow_depth3l, snow_mass3l, snow_tmp3l,    & !
!                          snow_rho3l, snow_rho1l, snow_age ) 
!  
!   use cable_data_mod, ONLY: cable_point_isnow, msn, ms
!
!   integer, DIMENSION(land_pts, ntiles), target :: isnow_flg3L
!   real :: TSOIL_TILE(land_pts,ntiles,ms)
!   real :: SMCL_TILE(land_pts,ntiles,ms)
!   real :: STHF_TILE(land_pts,ntiles,ms)
!   real :: SNOW_DEPTH3L(land_pts,ntiles,msn)
!   real :: SNOW_MASS3L(land_pts,ntiles,msn)
!   real :: SNOW_TMP3L(land_pts,ntiles,msn)
!   real :: SNOW_RHO3L(land_pts,ntiles,msn)
!   real :: SNOW_RHO1L(land_pts,ntiles)
!   real :: SNOW_AGE(land_pts,ntiles)
!   real :: SNOW_FLG3L(land_pts,ntiles)
!
!   call cable_point_isnow(isnow_flg3l, &
!                          tsoil_tile, smcl_tile, sthf_tile,     & !
!                          snow_depth3l, snow_mass3l, snow_tmp3l,    & !
!                          snow_rho3l, snow_rho1l, snow_age ) 
!  
!end subroutine cable_parse_isnow































