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
! Called from: JULES: surf_couple_radiation, tile_albedo_cable
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_rad_main_mod
  
contains

SUBROUTINE cable_rad_main( mype, timestep_number, cycleno, numcycles,          &
              row_length, rows, land_pts, ntiles,                              &
              sm_levels, tile_frac, surf_down_sw, cosine_zenith_angle, snow_tile,  &
              soil_alb, land_albedo, alb_surft, land_alb )
          !cable% snow_temp, 
          !cable% snow_avg_rho,                  &
          !cable% soil_temp, 
          !cable% snow_flg,&

  USE cable_rad_driv_mod, ONLY : cable_rad_driver
  USE atm_fields_real_mod, ONLY : snow_temp_cable, snow_avg_rho_cable, soil_temp_cable,&      
                                snow_flg_cable 
  USE cable_decs_mod, ONLY : sw_down_diag

  implicit none
  !--- IN ARGS FROM surf_couple_radiation() ------------------------------------
  
  integer :: mype, timestep_number, cycleno, numcycles
  integer :: row_length, rows, land_pts, ntiles, sm_levels ! grid
  real :: surf_down_sw(row_length,rows,4) ! 4-band ShortWave forcing
  real :: tile_frac(land_pts,ntiles)      ! surface type fraction 
  real :: cosine_zenith_angle(row_length,rows)        ! cosine_zenith_angle          
  real :: snow_tile(land_pts,ntiles)     ! formerly snow tile        
  real :: soil_alb(land_pts)              ! Snow-free, bare soil albedo: 
  real :: land_albedo(row_length,rows,4) 
  real :: alb_surft(Land_pts,ntiles,4)
  real :: land_alb(row_length,rows)         ! Mean land_albedo

  !--- End IN ARGS  -----------------------------------------------------------

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_rad_main"
  logical, save :: first_call = .true.
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable
  
  !--- End header -------------------------------------------------------------
  
  if(mype==0) then
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "@ timestep: ", timestep_number 
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
  if( .NOT. allocated(sw_down_diag) )allocate( sw_down_diag(row_length, rows, 4) ) 
  if(mype==0) then
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "calling _rad_driver: "
  endif
     
  call cable_rad_driver( cycleno, row_length, rows, land_pts, ntiles,          &
              sm_levels, tile_frac, & 
              snow_temp_cable, snow_avg_rho_cable, soil_temp_cable,            &
                                isnow_flg_cable,                               & 
              surf_down_sw, cosine_zenith_angle, snow_tile,                    &
              soil_alb, land_albedo, alb_surft, land_alb )
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  snow_flg_cable = real(isnow_flg_cable)
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- OUT args @ timestep X,Y,Z                                  -------------
  !----------------------------------------------------------------------------

  !jhan: call checks as required by namelis      
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  first_call = .false.        

  if(mype==0) then
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "End . @ timestep: ", timestep_number 
  endif
     
return

End subroutine cable_rad_main
  
End module cable_rad_main_mod











































