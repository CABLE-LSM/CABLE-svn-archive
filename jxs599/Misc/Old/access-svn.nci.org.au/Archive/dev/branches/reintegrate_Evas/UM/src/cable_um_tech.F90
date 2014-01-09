
   module cable_um_tech_mod
      use define_types
      implicit none

      TYPE (air_type), save             :: air
      TYPE (bgc_pool_type), save        :: bgc
      TYPE (met_type), save             :: met
      TYPE (balances_type), save        :: bal
      TYPE (radiation_type), save       :: rad
      TYPE (roughness_type), save       :: rough
      TYPE (soil_parameter_type), save  :: soil       ! soil parameters
      TYPE (soil_snow_type), save       :: ssoil
      TYPE (sum_flux_type), save        :: sum_flux
      TYPE (veg_parameter_type), save   :: veg        ! vegetation parameters
      TYPE (canopy_type), save          :: canopy
 
      type derived_rad_bands    
         real, allocatable ::                         &                                   
            SW_DOWN_DIR (:,:) & ! Surface downward SW direct radiation (W/m2).
            ,SW_DOWN_DIF(:,:) & ! Surface downward SW diffuse radiation (W/m2).
            ,SW_DOWN_VIS(:,:) & ! Surface downward VIS radiation (W/m2).
            ,SW_DOWN_NIR(:,:) & ! Surface downward NIR radiation (W/m2).
            ,FBEAM(:,:,:)      ! Surface downward SW radiation (W/m2).
      end type derived_rad_bands
      
      type um_dimensions 
         integer :: row_length, rows, land_pts, ntiles, npft, &
                     sm_levels, timestep 
         integer, allocatable, dimension(:) :: tile_pts, land_index
         integer, allocatable, dimension(:,:) :: tile_index
         real :: rho_water
         real,allocatable, dimension(:,:) :: tile_frac
         real,allocatable, dimension(:,:) :: latitude, longitude
         logical,allocatable, dimension(:,:) :: l_tile_pts
      endtype um_dimensions 

      type derived_veg_pars
         integer, dimension(:,:), pointer       &      
                           :: ivegt(:,:),    & ! vegetation  types
                              isoilm(:,:)      ! soil types
         real, dimension(:,:), pointer           &  
                            :: htveg(:,:),    &
                               laift(:,:)       !, hruffmax(:.:)
      end type derived_veg_pars

      interface check_nmlvar 
         module procedure check_chvar, check_intvar
      end interface check_nmlvar 
 
      type (derived_rad_bands), save  :: kblum_rad    
      type (derived_veg_pars), save  :: kblum_veg    
      type (um_dimensions), save :: um1
      
      real,allocatable, dimension(:) :: conv_rain_prevstep, conv_snow_prevstep 

   contains

   !========================================================================
   !========================================================================
   !========================================================================

      subroutine cable_um_runtime_vars(runtime_vars_file) 
         use io_variables, only: filename
         use cable_diag_module, only : cable_stat
         use cable_common_module, only : cable_runtime, cable_user, &
                                          cable_user, knode_gl
         implicit none
         character(len=*), intent(in) :: runtime_vars_file
         integer :: funit=88
         !--- namelist for CABLE runtime vars, files, switches 
         NAMELIST/CABLE/filename,cable_user
         
            if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
               call cable_stat('cable_um_runtime_vars')
            
            !--- assume namelist exists. no iostatus check 
            OPEN(unit=funit,FILE= runtime_vars_file)
               READ(funit,NML=CABLE)
               if (knode_gl==0) then
                  print *, '  '; print *, 'CABLE_log:' 
                  print *, '  Opened file - '
                  print *, '  ', trim(runtime_vars_file)
                  print *, '  for reading runtime vars.' 
                  print *, 'End CABLE_log:'; print *, '  '
              endif
            CLOSE(funit)
                         
            !--- check value of variable 
            call check_nmlvar('filename%veg', filename%veg)
            call check_nmlvar('filename%soil', filename%soil)
            call check_nmlvar('cable_user%DIAG_SOIL_RESP', cable_user%DIAG_SOIL_RESP)
            call check_nmlvar('cable_user%LEAF_RESPIRATION', cable_user%LEAF_RESPIRATION)
            call check_nmlvar('cable_user%FWSOIL_SWITCH', cable_user%FWSOIL_SWITCH)
            call check_nmlvar('cable_user%RUN_DIAG_LEVEL', cable_user%RUN_DIAG_LEVEL)

         return  
      end subroutine cable_um_runtime_vars

      !jhan: also add real, logical, int interfaces
      subroutine check_chvar(this_var, val_var)
         use cable_common_module, only : knode_gl
         implicit none
         character(len=*), intent(in) :: this_var, val_var 
            if (knode_gl==0) then
               print *, '  '; print *, 'CABLE_log:' 
               print *, '   run time variable - '
               print *, '  ', trim(this_var) 
               print *, '   defined as - '
               print *, '  ', trim(val_var) 
               print *, 'End CABLE_log:'; print *, '  '
            endif
     
         return
      end subroutine check_chvar
 
      subroutine check_intvar(this_var, val_var)
         use cable_common_module, only : knode_gl
         implicit none
         character(len=*), intent(in) :: this_var
         integer, intent(in) :: val_var 

            if (knode_gl==0) then
               print *, '  '; print *, 'CABLE_log:' 
               print *, '   run time variable - '
               print *, '  ', trim(this_var) 
               print *, '   defined as - '
               print *, '  ', val_var
               print *, 'End CABLE_log:'; print *, '  '
            endif
     
         return
      end subroutine check_intvar
    
      !========================================================================= 
      !=========================================================================
      !========================================================================= 
 
   subroutine alloc_um_interface_types( row_length, rows, land_pts, ntiles, sm_levels)
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
      integer,intent(in) :: row_length, rows, land_pts, ntiles, sm_levels   

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('alloc_um_interface_types')
         
         allocate( um1%land_index(land_pts) )
         allocate( um1%tile_pts(ntiles) )
         allocate( um1%tile_frac(land_pts, ntiles) )
         allocate( um1%tile_index(land_pts, ntiles) )
         allocate( um1%latitude(row_length, rows) )
         allocate( um1%longitude(row_length, rows) )
         allocate( um1%l_tile_pts(land_pts, ntiles) ) 
        !-------------------------------------------------------
         allocate ( kblum_rad%sw_down_dir(row_length,rows) )
         allocate ( kblum_rad%sw_down_dif(row_length,rows) )
         allocate ( kblum_rad%sw_down_vis(row_length,rows) )
         allocate ( kblum_rad%sw_down_nir(row_length,rows) )
         allocate ( kblum_rad%fbeam(row_length,rows,3) )
         allocate( kblum_veg%htveg(land_pts,ntiles) )
         allocate( kblum_veg%laift(land_pts,ntiles) )
         allocate( kblum_veg%ivegt(land_pts,ntiles) )
         allocate( kblum_veg%isoilm(land_pts,ntiles) ) 
         !-------------------------------------------------------
         
      return
   end subroutine alloc_um_interface_types 

!========================================================================
!========================================================================
!========================================================================

   subroutine dealloc_vegin_soilin()
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user, vegin, soilin
      implicit none
         
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('dealloc_vegin_soilin')
         
         deallocate(vegin%canst1)
         deallocate(vegin%dleaf)
         deallocate(vegin%vcmax)
         deallocate(vegin%ejmax)
         deallocate(vegin%hc)
         deallocate(vegin%xfang)
         deallocate(vegin%rp20)
         deallocate(vegin%rpcoef)
         deallocate(vegin% rs20)
         deallocate(vegin%shelrb)
         deallocate(vegin%vegcf)
         deallocate(vegin%frac4)
         deallocate(vegin%reflin)
         deallocate(vegin%taulin)
         deallocate(vegin%xalbnir)
         deallocate(vegin%extkn)
         deallocate(vegin%froot)
         deallocate(vegin%tminvj)
         deallocate(vegin%tmaxvj)
         deallocate(vegin%vbeta)
         deallocate(vegin%cplant)
         deallocate(vegin%csoil)
         deallocate(vegin%ratecp)
         deallocate(vegin%ratecs)
        
         deallocate(soilin%silt)
         deallocate(soilin%clay)
         deallocate(soilin%sand)
         deallocate(soilin%swilt)
         deallocate(soilin%sfc)
         deallocate(soilin%ssat)
         deallocate(soilin%bch)
         deallocate(soilin%hyds)
         deallocate(soilin%sucs)
         deallocate(soilin%rhosoil)
         deallocate(soilin%css)

      return

   end subroutine dealloc_vegin_soilin


   !========================================================================
   !========================================================================
   !========================================================================


   end module cable_um_tech_mod




