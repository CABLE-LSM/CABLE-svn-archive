#include "include/cable_directives.h"

module cable_um_init_mod

    public :: interface_UM_data
      
   contains

   !===========================================================================!
   !=== called from cable_um_explicit_driver to initialize/update CABLE from ==! 
   !=== UM forcings etc.                                                    ===!   
   !===---------------------------------------------------------------------===!
   !=== comprises CALLS to 2 subrs.                                         ===!
   !=== 1. setup_um_interface_types - technical allocation/assignment of    ===!
   !===    temporary types, grouping UM vars                                ===!  
   !=== 2. intialize_cable_vars - initializes CABLE from UM vars            ===!  
   !===========================================================================!

   subroutine interface_UM_data(    &
#                    include "include/args/um_expl_pack.h"
                  &  )  

      !___"types" of similar UM land surface vars "typed" for convenience. 
      use cable_um_tech_mod, only : alloc_um_interface_types,alloc_vegin_soilin, &
                                    dealloc_vegin_soilin, um1, kblum_veg      
      use cable_um_init_subrs
      use cable_common_module
      use define_dimensions, only : mp
      use cable_diag_module, only : cable_stat
      implicit none            
      
      !___ recieved args
#     include "include/decs/um_expl_pack.h"
    
      !___ local decs
      !--- set in cable_directives.h (from read files pref. [mvtype,mstype]) 
      integer,parameter :: nsoilin= JHNSOIL 
      integer,parameter :: nvegin = JHNTILES 
      integer :: timestep
      !___defs 1st call to CABLE in this run. necesary in UM & coupled.
      logical, save :: first_call = .true. 
      integer :: i,j

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('interface_UM_data')
         timestep = itimestep   
        
         !---------------------------------------------------------------------!
         !--- code to create type um1% conaining UM basic vars describing    --! 
         !--- dimensions etc, which are used frequently throughout interfaces. !
         !--- also allocate mem for setting veg/soil from external files    ---!
         !---------------------------------------------------------------------!
               
         if(first_call) then
            call alloc_um_interface_types( row_length, rows, land_pts, &
                        ntiles, sm_levels )
            call alloc_vegin_soilin(nvegin,nsoilin)
         endif       
         
         call assign_um_fundamentals_to_um1( row_length, rows, land_pts, ntiles, &
                                       npft, sm_levels, timestep, latitude, &
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
                  if( um1%TILE_FRAC(i,j) .gt. 0.0 ) & 
                        um1%L_TILE_PTS(i,j) = .TRUE.
               enddo
            enddo
         endif

         !--- initialize soil
         call initialize_soil( & 
#                    include "include/args/um_soil.h"
                  &  )
           
         !--- initialize veg   
         call initialize_veg( canht_ft, lai_ft ) 
    
         !--- initialize roughness
         call initialize_roughness( & 
#                       include "include/args/um_met.h"
                        &, kblum_veg%htveg ) 
          
         !--- initialize soilsnow
         call initialize_soilsnow(   &
#                    include "include/args/um_soil.h"
                     , &
#                    include "include/args/um_snow.h"                  
                   & , sin_theta_latitude ) 
 
         !--- initialize canopy   
         call initialize_canopy(CANOPY_TILE)
    
         !--- initialize radiation & met forcing
         call initialize_radiation( &
#                       include "include/args/um_rad.h"
                        , &
#                       include "include/args/um_met.h"
                      & )    
    
         if (first_call) then
            call init_bgc_vars() 
            call init_sumflux_zero() 
            call dealloc_vegin_soilin(nvegin,nsoilin)
            first_call = .false. 
         endif      
          
         return
      end subroutine interface_UM_data
                                   
   !============================================================================
   !============================================================================
   !============================================================================

   subroutine assign_um_fundamentals_to_um1( row_length, rows, land_pts,    &
         ntiles, npft, sm_levels, timestep, latitude, longitude, land_index,     &
         tile_frac, tile_pts, tile_index, l_tile_pts, rho_water  )
      use cable_diag_module, only : cable_stat
      use cable_um_tech_mod, only : um1
      use cable_common_module, only : cable_runtime, cable_user
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
            call cable_stat('assign_um_fundamentals_to_um1')
         
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
         !um1%l_tile_pts= l_tile_pts
         um1%rho_water= rho_water 

        return 
     end subroutine assign_um_fundamentals_to_um1

   !============================================================================
   !============================================================================
   !============================================================================

end module cable_um_init_mod
