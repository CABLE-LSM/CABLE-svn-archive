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

SUBROUTINE cable_hyd_main( land_pts, ntiles, lying_snow, SNOW_surft, SURF_ROFF,&
                           SUB_SURF_ROFF, TOT_TFALL )
  
  !subrs called 
  USE cable_hyd_driv_mod, ONLY : cable_hyd_driver
  
  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime

  implicit none
 
  !___ re-decl input args

  integer :: land_pts, ntiles
  
  real :: snow_surft(land_pts,ntiles)
  
  real, dimension(land_pts) ::                                                 &
    lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
    sub_surf_roff,  & ! OUT Sub-surface runoff (kg/m2/s).
    surf_roff,      & ! OUT Surface runoff (kg/m2/s).
    tot_tfall         ! OUT Total throughfall (kg/m2/s).

  !___ local vars
  logical, save :: first_call = .true.
  
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_hyd_main"
 
  !-------- Unique subroutine body -----------
  
  !--- initialize cable_runtime% switches 
  cable_runtime%um =          .TRUE.
  cable_runtime%um_hydrology =.TRUE.
  
  CALL cable_hyd_driver( land_pts, ntiles, lying_snow, SNOW_surft, SURF_ROFF,   &
                         SUB_SURF_ROFF, TOT_TFALL )
  
  cable_runtime%um_hydrology =.FALSE.
  
  !-------- End Unique subroutine body -----------
  
  first_call = .false.        

return

End subroutine cable_hyd_main
  
End module cable_hyd_main_mod











































