!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
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
! Purpose: Initialize and update CABLE variables from UM forcing, calls to 
!          memory allocation and initialization subroutines
!
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: No significant change from v1.8
!
!
! ==============================================================================

MODULE cable_um_init_mod

   IMPLICIT NONE   

CONTAINS

SUBROUTINE explicit_call_initialization(                                       & 
            ! cable% type yet to be officailly implemented
            row_length, & ! -> cable%row_length 
            rows,       & ! -> cable%rows
            land_pts,   & ! -> cable%land_pts
            ntiles,     & ! -> cable%ntiles
            npft,       & ! -> cable%npft 
            sm_levels,  & ! -> cable%ms
            timestep,   & ! -> cable%timestep_width
            latitude,   & ! -> cable%latitude
            longitude,  & ! -> cable%longitude
            land_index, & ! -- necessary for packing 
            tile_frac,  & ! -> cable%tile_frac
            tile_pts,   & ! -> cable%
            tile_index, & ! -- necessary for packing
            dzsoil,     & ! -> soil%zse                        
            ! soil properties from UM/JULES
            bexp,       & ! -> soil%bch
            hcon,       & ! ~> soil%cnsd
            satcon,     & ! ~> soil%hyds
            sathh,      & ! -> soil%sucs
            smvcst,     & ! -> soil%ssat
            smvcwt,     & ! -> soil%swilt
            smvccl,     & ! -> soil%sfc
            albsoil,    & ! -> soil%albsoil
            ! canopy properties from UM/JULES
            canht_ft,   & ! ~> veg%hc
            lai_ft,     & ! ~> veg%lai
            ! forcing from JULES
            sw_down,    & ! ~> met%fsd
            lw_down,    & ! -> met%fld 
            ls_rain,    & ! ~> met%precip
            ls_snow,    & ! ~> met%precip_sn
            tl_1,       & ! -> met%tk
            qw_1,       & ! -> met%qv
            vshr_land,  & ! -> met%ua
            pstar,      & ! ~> met%pmb
            z1_tq,      & ! -> rough%za_tq
            z1_uv,      & ! -> rough%za_uv
            canopy_tile,& ! -> canopy%cansto
            Fland,      & ! -> ssnow%fland
            CO2_MMR,    & ! ~> met%ca
         
            !jhan:adapted from JULES var cosz, done elsewhere, move to here 
            ! and make switchable
            cos_zenith_angle,    & ! ->met%coszen
         
            ! snow properties from UM/JULES
            snow_tile,     & ! -> ssnow%snowd
         
            ! snow properties from CABLE vars 
            snage_tile,    & ! -> ssnow%snage
            snow_rho1l,    & ! -> ssnow%ssdnn
            isnow_flg3l,   & ! -> ssnow%isflag
            snow_rho3l,    & ! -> ssnow%ssdn
            snow_depth3l,  & ! -> ssnow%sdepth
            snow_tmp3l,    & ! -> ssnow%tggsn
            snow_mass3l,   & ! -> ssnow%smass
            snow_cond,     & ! -> ssnow%sconds
         
         
            !soil properties from CABLE vars 
            smcl_tile,     & ! ~> soil%wb
            sthf_tile,     & ! ~> soil%wbice
            tsoil_tile     & ! -> ssnow%tgg
   )                         

   USE cable_um_init_subrs_mod     ! where most subrs called from here reside
   
   USE cable_um_tech_mod,   ONLY :                                             &
      alloc_um_interface_types,  & ! mem. allocation subr (um1, kblum%) 
      dealloc_vegin_soilin,      & ! mem. allocation subr (vegin%,soilin%)
      um1,                       & ! um1% type UM basics 4 convenience
      kblum_veg                    ! kblum_veg% reset UM veg vars 4 CABLE use

   USE cable_common_module, ONLY :                                             &
      cable_user,       & ! cable_user% type inherits user definition
                          ! via namelist (cable.nml) 
      get_type_parameters ! veg and soil parameters READ subroutine  

   USE cable_def_types_mod, ONLY : mp ! number of points CABLE works on


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 
   !___UM dimensions, array indexes, flags
   INTEGER, INTENT(IN) ::                                                      &
      row_length, rows, & ! UM resolution
      land_pts,         & ! number of land_pts
      ntiles,           & ! number of tiles
      npft,             & ! number of Plant Functional Types
      sm_levels           ! number of soil layers

   INTEGER, INTENT(IN), DIMENSION(land_pts) ::                                 &
      land_index  ! index of land point 
   
   INTEGER, INTENT(IN), DIMENSION(ntiles) ::                                   &
      tile_pts    ! number of land_pts per tile type
  
   INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles) ::                         &
      tile_index, &  ! index of tile 
      isnow_flg3l    ! flag for 3-layer snow 

   !___UM parameters 
   INTEGER, INTENT(IN) :: timestep
   REAL, INTENT(IN), DIMENSION(sm_levels) ::                                   &
      dzsoil

   !___UM soil/snow/radiation/met vars
   REAL, INTENT(IN), DIMENSION(land_pts) ::                                    &
      bexp, &     !
      hcon, &     !  
      satcon, &   !
      sathh,  &   !
      smvcst, &   !
      smvcwt, &   !
      smvccl, &   !
      albsoil,&   !
      fland       !
       
   REAL, INTENT(INOUT), DIMENSION(row_length,rows) ::                          &
      sw_down, &        !
      cos_zenith_angle  !

   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
      latitude, longitude, &
      lw_down, &  !
      ls_rain, &  !
      ls_snow, &  !   
      tl_1,    &  !
      qw_1,    &  !
      vshr_land,& !
      pstar,   &  !
      z1_tq,   &  !
      z1_uv       !

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         &
      snow_tile   !
   
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            & 
      tile_frac, &   !   
      snow_rho1l,&   !
      snage_tile     !

   REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
      canht_ft, & !
      lai_ft      !

   REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile !

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond   !

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l, &     !
      snow_depth3l, &   ! 
      snow_mass3l,  &   ! 
      snow_tmp3l

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles, sm_levels) ::                 &
      sthf_tile, &   !
      smcl_tile, &   !
      tsoil_tile     !

   REAL, INTENT(IN) :: co2_mmr

   LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::                  &
      L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
  
   REAL, DIMENSION(row_length,rows) ::                             & 
      sin_theta_latitude

   !------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM cable_explicit_driver() -------------------------
   !------------------------------------------------------------------------- 

   !___ local decs
   !___defs 1st call to CABLE in this run. necesary in UM & coupled.
   LOGICAL, SAVE :: first_call = .TRUE. 
   INTEGER :: i,j

   INTEGER, DIMENSION(land_pts, ntiles) ::                                     &
     tile_index_mp   !tile # on each land point, ntile for packing to CABLE

   !--- logn, vegparmnew can be set thru cable.nml
   INTEGER :: logn=6       ! 6=write to std out
   LOGICAL :: vegparmnew=.true.   ! true=read std veg params false=CASA file 
         

   
   IF( .NOT. ALLOCATED(L_tile_pts) )                                           & 
      ALLOCATE( L_tile_pts(land_pts, ntiles) )
     
   !jhan: move this into interface_UM_data
   sin_theta_latitude  = sin(latitude)

      !---------------------------------------------------------------------!
      !--- code to create type um1% conaining UM basic vars describing    --! 
      !--- dimensions etc, which are used frequently throughout interfaces. !
      !--- and then assign every time explicit is called (shouldn't need to)!
      !---------------------------------------------------------------------!
      IF(first_call) THEN
         CALL alloc_um_interface_types( row_length, rows, land_pts,            &
                                        ntiles, sm_levels )
      ENDIF       
      
      CALL assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,        &
                                    npft, sm_levels, timestep,            &
                                    latitude,      &
                                    longitude, land_index, tile_frac,          &
                                    tile_pts, tile_index &
                                  )




      !---------------------------------------------------------------------!
      !--- CABLE vars are initialized/updated from passed UM vars       ----!
      !---------------------------------------------------------------------!
      
      !---def. vector length for cable(mp) & logical l_tile_pts
      !--- IF the tile is "active"
      IF ( first_call ) THEN
      
         L_TILE_PTS = .FALSE.
         mp = SUM(TILE_PTS)
         
         CALL alloc_cable_types()
         
         DO i=1,land_pts
            DO j=1,ntiles
               
               IF( TILE_FRAC(i,j) .GT. 0.0 ) THEN 
                     L_TILE_PTS(i,j) = .TRUE.
                  !jhan:can set veg%iveg from  here ?
                  tile_index_mp(i,j) = j 
               ENDIF
            
            ENDDO
         ENDDO
         
         um1%L_TILE_PTS = L_TILE_PTS 
      
      ENDIF
         
      !jhan: turn this off until implementation finalised
      !--- initialize latitude/longitude & mapping IF required
      !if ( first_call ) & 
      !   call initialize_maps(latitude,longitude, tile_index_mp)



      !--- read in soil (and veg) parameters 
      IF(first_call)                                                        & 
         CALL  get_type_parameters(logn,vegparmnew)

      !--- initialize veg   
      CALL initialize_veg( canht_ft, lai_ft ) 
 
      !--- initialize soil
      CALL initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt,      &
                            smvccl, albsoil,                                &
                            dzsoil ) 
        
      !--- initialize roughness
      CALL initialize_roughness( z1_tq, z1_uv, kblum_veg%htveg ) 
       
      !--- initialize soilsnow
      CALL initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile,   &
                                snow_tile, snow_rho1l, snage_tile,          &
                                isnow_flg3l, snow_rho3l, snow_cond,         &
                                snow_depth3l, snow_mass3l, snow_tmp3l,      &
                                fland, sin_theta_latitude ) 

      !--- initialize canopy   
      CALL initialize_canopy(CANOPY_TILE)
 
      !--- initialize radiation & met forcing
      CALL initialize_radiation( sw_down, lw_down, cos_zenith_angle,        &
                                 ls_rain, &
                                 ls_snow, tl_1, qw_1, vshr_land, pstar,     &
                                 co2_mmr ) 
                   
 
      IF( first_call ) THEN
         CALL init_bgc_vars() 
         CALL init_sumflux_zero() 
         CALL dealloc_vegin_soilin()
         first_call = .FALSE. 
      ENDIF      
      
END SUBROUTINE explicit_call_initialization
 
!============================================================================
!============================================================================
!============================================================================


SUBROUTINE implicit_call_initialization(                                       & 
            row_length, rows, & ! grid resolution
            ls_rain,    & ! ~>  met%precip
            ls_snow,    & ! ~>  met%precip_sn
            conv_rain,  & ! ~~> met%precip
            conv_snow,  & ! ~~> met%precip_sn
            dtl_1,      & !~~> met%precip_sn
            dqw_1       & !~~> met%precip_sn
   )
    
   USE cable_def_types_mod, ONLY : mp
   USE cable_um_tech_mod, ONLY : um1, met
   USE cable_um_init_subrs_mod, ONLY : um2cable_rr
        
   INTEGER, INTENT(IN) ::                                                      &
      row_length, rows

   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
      ls_rain,    &
      ls_snow     

   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                         &
      CONV_RAIN, & ! IN Convective rain
      CONV_SNOW   ! IN Convective snow
   
   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                         &
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   REAL, DIMENSION(:),  ALLOCATABLE, SAVE ::                                   & 
      dtlc, & 
      dqwc, &
      conv_rain_prevstep, &
      conv_snow_prevstep
     
   LOGICAL, SAVE :: first_call = .TRUE. 

      IF(first_call) THEN
         ALLOCATE( &
                   dtlc(mp),               & 
                   dqwc(mp),               &
                   conv_rain_prevstep(mp), &
                   conv_snow_prevstep(mp)  &
                 )  
         first_call = .FALSE. 
      ENDIF

      dtlc = 0. ; dqwc = 0.
      
      ! should be no change from explicit call offline
      CALL um2cable_rr( (LS_RAIN+CONV_RAIN)*um1%TIMESTEP, met%precip)
      CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
      CALL um2cable_rr( dtl_1, dtlc)
      CALL um2cable_rr( dqw_1, dqwc)
      
      ! in JULES this increment is meaningless
      !jhan
      dtlc = 0. ; dqwc = 0.
      
      !--- conv_rain(snow)_prevstep are added to precip. in explicit call
      CALL um2cable_rr( (CONV_RAIN)*um1%TIMESTEP, conv_rain_prevstep)
      CALL um2cable_rr( (CONV_snow)*um1%TIMESTEP, conv_snow_prevstep)
      
      met%precip   =  met%precip + met%precip_sn
      met%tk = met%tk + dtlc
      met%qv = met%qv + dqwc
      met%tvair = met%tk
      met%tvrad = met%tk
   
END SUBROUTINE implicit_call_initialization                                  
!
!============================================================================
!============================================================================
!============================================================================

SUBROUTINE assign_um_basics_to_um1( row_length, rows, land_pts, ntiles,     &
                                    npft, sm_levels, timestep, latitude,    &
                                    longitude, land_index, tile_frac,       &
                                    tile_pts, tile_index                    &
                                  )
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_user

   INTEGER, INTENT(IN) :: row_length, rows, land_pts, ntiles, npft, sm_levels
   INTEGER, INTENT(IN) :: timestep 
   REAL, INTENT(IN), DIMENSION(row_length,rows) :: latitude, longitude 
   INTEGER, INTENT(IN), DIMENSION(land_pts)  :: land_index 
   INTEGER, INTENT(IN), DIMENSION(ntiles)  :: tile_pts 
   INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_index
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles)  :: tile_frac 
     
      um1%row_length = row_length
      um1%rows = rows
      um1%land_pts = land_pts
      um1%ntiles = ntiles   
      um1%npft = npft   
      um1%sm_levels = sm_levels
      um1%timestep = timestep
      um1%latitude = latitude
      um1%longitude = longitude
      um1%land_index = land_index
      um1%tile_frac = tile_frac
      um1%tile_pts = tile_pts
      um1%tile_index = tile_index

END SUBROUTINE assign_um_basics_to_um1

!============================================================================
!============================================================================
!============================================================================

END MODULE cable_um_init_mod
