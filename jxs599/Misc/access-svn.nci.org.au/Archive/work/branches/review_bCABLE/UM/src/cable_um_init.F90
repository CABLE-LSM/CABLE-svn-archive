
module cable_um_init_mod

    public :: interface_UM_data
    private  
   contains

   !===========================================================================!
   !=== called from cable_um_explicit_driver to initialize/update CABLE from ==! 
   !=== UM forcings etc.                                                    ===!   
   !=== comprises CALLS to mem. allocation & initialization subrs.          ===!
   !===========================================================================!

   subroutine interface_UM_data( row_length, rows, land_pts, ntiles,    & 
               npft, sm_levels, itimestep, latitude, longitude, land_index,      &
               tile_frac, tile_pts, tile_index, bexp, hcon, satcon, sathh,       &
               smvcst, smvcwt, smvccl, albsoil, snow_tile, snow_rho1l,           &
               snage_tile, isnow_flg3l, snow_rho3l, snow_cond, snow_depth3l,     &
               snow_tmp3l, snow_mass3l, sw_down, lw_down, cos_zenith_angle,      &
               surf_down_sw, ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar,     &
               z1_tq, z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,          &
               CO2_MMR, sthu_tile, smcl_tile, sthf_tile, sthu, tsoil_tile,       &
               canht_ft, lai_ft, sin_theta_latitude, dzsoil )                         

      use cable_um_init_subrs          ! where most subrs called from here reside
      
      use cable_um_tech_mod, only : &
         alloc_um_interface_types,  &  ! mem. allocation subr (um1, kblum%) 
         dealloc_vegin_soilin,      &  ! mem. allocation subr (vegin%,soilin%)
         um1,      &                   ! um1% type UM basics 4 convenience
         kblum_veg                     ! kblum_veg% reset UM veg properties 4 CABLE

      use cable_common_module, only : &
         cable_user, &  ! cable_user% type has user def from cable.nml 
         get_type_parameters

      use define_dimensions, only : &
         mp       ! number of points CABLE works on

      use cable_diag_module, only : &
         cable_stat ! basic DIAG to report progress
      
      implicit none            
      
      !-------------------------------------------------------------------------- 
      !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
      !-------------------------------------------------------------------------- 

      !___UM dimensions, array indexes, flags
      integer, intent(in) :: &
         row_length, rows, &  ! UM resolution
         land_pts, &          ! number of land_pts
         ntiles, &            ! number of tiles
         npft, &              ! number of Plant Functional Types
         sm_levels            ! number of soil layers

      integer, intent(in), dimension(land_pts) :: &
         land_index  ! index of land point 
      
      integer, intent(in), dimension(ntiles) :: &
         tile_pts    ! number of land_pts per tile type
     
      integer, intent(in), dimension(land_pts, ntiles) :: &
         tile_index, &  ! index of tile 
         isnow_flg3l    ! flag for 3-layer snow 

      !___UM parameters 
      integer, intent(in) :: itimestep
      real, intent(in) :: rho_water 
      real, intent(in), dimension(sm_levels) :: &
         dzsoil

      !___UM soil/snow/radiation/met vars
      real, intent(in), dimension(land_pts) :: &
         bexp, &     !
         hcon, &     !  
         satcon, &   !
         sathh,  &   !
         smvcst, &   !
         smvcwt, &   !
         smvccl, &   !
         albsoil,&   !
         fland       !
          
      real, intent(inout), dimension(row_length,rows) :: &
         sw_down, &        !
         cos_zenith_angle  !

      real, intent(in), dimension(row_length,rows) :: &
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

      real, intent(inout), dimension(land_pts, ntiles) :: &
         snow_tile   !
      
      real, intent(in), dimension(land_pts, ntiles) :: & 
         tile_frac, &   !   
         snow_rho1l,&   !
         snage_tile     !

      real, intent(in), dimension(row_length, rows, 4) :: &
         surf_down_sw 

      real, intent(in), dimension(land_pts, npft) :: &
         canht_ft, & !
         lai_ft      !

      real, intent(in),dimension(land_pts, ntiles) :: &
         canopy_tile !

      real, intent(inout), dimension(land_pts, ntiles,3) :: &
         snow_cond   !

      real, intent(in), dimension(land_pts, ntiles,3) :: &
         snow_rho3l, &     !
         snow_depth3l, &   ! 
         snow_mass3l,  &   ! 
         snow_tmp3l

      real, intent(in), dimension(land_pts, sm_levels) :: &
         sthu  !

      real, intent(in), dimension(land_pts, ntiles, sm_levels) :: &
         sthu_tile, &   !
         sthf_tile, &   !
         smcl_tile, &   !
         tsoil_tile     !

      real, intent(in) :: co2_mmr

      logical, intent(inout),dimension(land_pts, ntiles) :: &
         L_tile_pts  ! true IF vegetation (tile) fraction is greater than 0
     
      real, intent(in), dimension(row_length,rows) :: & 
         sin_theta_latitude

      !-------------------------------------------------------------------------- 
      !--- end INPUT ARGS FROM cable_explicit_driver() ---------------------------
      !-------------------------------------------------------------------------- 

      !___ local decs
      !___defs 1st call to CABLE in this run. necesary in UM & coupled.
      logical, save :: first_call = .true. 
      integer :: i,j

      integer, dimension(land_pts, ntiles) :: &
        tile_index_mp   !tile # on each land point, ntile for packing to CABLE

      !--- logn, vegparmnew can be set thru cable.nml
      integer :: logn=6       !--- 6=write to std out
      logical :: vegparmnew=.true.   !--- true=read std veg params false=CASA file 
            


         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('interface_UM_data')
        
         !---------------------------------------------------------------------!
         !--- code to create type um1% conaining UM basic vars describing    --! 
         !--- dimensions etc, which are used frequently throughout interfaces. !
         !--- and then assign every time explicit is called (shouldn't need to)!
         !---------------------------------------------------------------------!
         if(first_call) then
            call alloc_um_interface_types( row_length, rows, land_pts, &
                        ntiles, sm_levels )
         endif       
         
         call assign_um_basics_to_um1( row_length, rows, land_pts, ntiles, &
                                       npft, sm_levels, itimestep, latitude, &
                                       longitude, land_index, tile_frac, tile_pts, &
                                       tile_index, l_tile_pts, rho_water  )



         !---------------------------------------------------------------------!
         !--- CABLE vars are initialized/updated from passed UM vars       ----!
         !---------------------------------------------------------------------!
         
         !---def. vector length for cable(mp) & logical l_tile_pts
         !--- IF the tile is "active"
         IF ( first_call ) THEN
            um1%L_TILE_PTS = .FALSE.
            mp = sum(um1%TILE_PTS)
            call alloc_cable_types()
            do i=1,land_pts
               do j=1,ntiles
                  if( um1%TILE_FRAC(i,j) .gt. 0.0 ) then 
                        um1%L_TILE_PTS(i,j) = .TRUE.
                     !jhan:can set veg%iveg from  here ?
                     tile_index_mp(i,j) = j 
                  endif
               enddo
            enddo
         endif
            
         !jhan: turn this off until properly implemented here   
         !--- initialize latitude/longitude & mapping IF required
         !if ( first_call ) & 
         !   call initialize_maps(latitude,longitude, tile_index_mp)



         !--- read in soil (and veg) parameters 
         if (first_call) & 
            call  get_type_parameters(logn,vegparmnew)

         !--- initialize veg   
         call initialize_veg( canht_ft, lai_ft ) 
    
         !--- initialize soil
         call initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt, &
                  smvccl, albsoil, tsoil_tile, sthu, sthu_tile, dzsoil ) 
           
         !--- initialize roughness
         call initialize_roughness( z1_tq, z1_uv, kblum_veg%htveg ) 
          
         !--- initialize soilsnow
         call initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile, &
                  snow_tile, snow_rho1l, snage_tile, isnow_flg3l, snow_rho3l,   &
                  snow_cond, snow_depth3l, snow_mass3l, snow_tmp3l, fland, & 
                  sin_theta_latitude ) 
 
         !--- initialize canopy   
         call initialize_canopy(CANOPY_TILE)
    
         !--- initialize radiation & met forcing
         call initialize_radiation( sw_down, lw_down, cos_zenith_angle, &
                  surf_down_sw, sin_theta_latitude, ls_rain, ls_snow,   &
                  tl_1, qw_1, vshr_land, pstar, co2_mmr ) 
                      
    
         if (first_call) then
            call init_bgc_vars() 
            call init_sumflux_zero() 
            call dealloc_vegin_soilin()
            first_call = .false. 
         endif      
         
         
          
         return
      end subroutine interface_UM_data
                                   
   !============================================================================
   !============================================================================
   !============================================================================

   subroutine assign_um_basics_to_um1( row_length, rows, land_pts,    &
         ntiles, npft, sm_levels, timestep, latitude, longitude, land_index,     &
         tile_frac, tile_pts, tile_index, l_tile_pts, rho_water  )
      use cable_diag_module, only : cable_stat
      use cable_um_tech_mod, only : um1
      use cable_common_module, only : cable_user
      implicit none

      integer, intent(in) :: row_length, rows, land_pts, ntiles, npft, sm_levels
      integer, intent(in) :: timestep 
      real, intent(in), dimension(row_length,rows) :: latitude, longitude 
      real,intent(in):: rho_water
      integer, intent(in), dimension(land_pts)  :: land_index 
      integer, intent(in), dimension(ntiles)  :: tile_pts 
      integer, intent(in), dimension(land_pts, ntiles)  :: tile_index
      real, intent(in), dimension(land_pts, ntiles)  :: tile_frac 
      logical, intent(in), dimension(land_pts,ntiles)  :: l_tile_pts 
        
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('assign_um_basics_to_um1')
         
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
         um1%rho_water= rho_water 

        return 
     end subroutine assign_um_basics_to_um1

   !============================================================================
   !============================================================================
   !============================================================================

end module cable_um_init_mod
