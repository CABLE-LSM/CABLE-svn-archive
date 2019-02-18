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

module cable_hyd_main_mod
  
contains

SUBROUTINE cable_hyd_main( land_pts, ntiles, lying_snow, SNOW_surft, SURF_ROFF, SUB_SURF_ROFF,  &
                             TOT_TFALL )
  
  USE cable_common_module, ONLY : knode_gl,        & ! processor number
                                  ktau_gl,         & ! number
                                  kwidth_gl          ! width in S 
  
  USE cable_hyd_driv_mod, ONLY : cable_hyd_driver

  implicit none
 
  !--- IN ARGS FROM sf_exch_cable, passed from surf_couple_hyd() down ----
  
  integer :: land_pts, ntiles
  
  real :: snow_surft(land_pts,ntiles)
  
  real, dimension(land_pts) ::                                                 &
    lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
    sub_surf_roff,  & ! OUT Sub-surface runoff (kg/m2/s).
    surf_roff,      & ! OUT Surface runoff (kg/m2/s).
    tot_tfall         ! OUT Total throughfall (kg/m2/s).

  !--- End IN ARGS  -----------------------------------------------------------

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_hyd_main"
  logical, save :: first_call = .true.
  
  !--- End header -------------------------------------------------------------
  
  if(knode_gl==0) then
    write (6, *) "CABLE_LSM:Subr: ", subr_name,  "@ timestep: ",ktau_gl 
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
  
  CALL cable_hyd_driver( land_pts, ntiles, lying_snow, SNOW_surft, SURF_ROFF,   &
                         SUB_SURF_ROFF, TOT_TFALL )
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
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

return

End subroutine cable_hyd_main
  
End module cable_hyd_main_mod











































