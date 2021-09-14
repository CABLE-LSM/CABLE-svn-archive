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
! Purpose: Reads vegetation and soil parameter files, fills vegin, soilin
!          NB. Most soil parameters overwritten by spatially explicit datasets
!          input as ancillary file (for ACCESS) or surface data file (for offline)
!          Module enables accessibility of variables throughout CABLE
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: v2.0 vegin%dleaf now calculated from leaf length and width
!          Parameter files were read elsewhere in v1.8 (init_subrs)
!
!
! ==============================================================================

MODULE cable_common_module
   IMPLICIT NONE 

   !---allows reference to "gl"obal timestep in run (from atm_step)
   !---total number of timesteps, and processing node 
   INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl
   
   ! set from environment variable $HOME
   CHARACTER(LEN=200) ::                                                       & 
      myhome

   !---Lestevens Sept2012
   !---CASACNP switches and cycle index
   LOGICAL, SAVE :: l_casacnp,l_laiFeedbk,l_vcmaxFeedbk
   LOGICAL :: l_luc = .FALSE.
   LOGICAL :: l_thinforest = .FALSE.
   
   !---CABLE runtime switches def in this type
   TYPE kbl_internal_switches
      LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
            um_radiation = .FALSE.
      LOGICAL :: offline = .FALSE., mk3l = .FALSE.
   END TYPE kbl_internal_switches 

   TYPE(kbl_internal_switches), SAVE :: cable_runtime

   !---CABLE runtime switches def in this type
   TYPE kbl_user_switches
      !jhan: this is redundant now we all use filename%veg?
      CHARACTER(LEN=200) ::                                                    &
         VEG_PARS_FILE  ! 
      
      CHARACTER(LEN=20) ::                                                     &
         FWSOIL_SWITCH     !
      
      CHARACTER(LEN=5) ::                                                      &
         RUN_DIAG_LEVEL  !
      
      CHARACTER(LEN=3) ::                                                      &
         SSNOW_POTEV,      & !
         DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical) 
         LEAF_RESPIRATION    ! either ON or OFF (jhan:Make Logical) 

      LOGICAL ::                                                               &
         INITIALIZE_MAPPING = .FALSE., & ! 
         CONSISTENCY_CHECK = .FALSE.,  & !
         CASA_DUMP_READ = .FALSE.,     & !
         CASA_DUMP_WRITE = .FALSE.,    & !
         CABLE_RUNTIME_COUPLED = .FALSE., & !
         ! L.Stevens - Test Switches
         L_NEW_ROUGHNESS_SOIL  = .FALSE., & !
         L_NEW_RUNOFF_SPEED    = .FALSE., & !
         L_NEW_REDUCE_SOILEVP  = .FALSE.!


   END TYPE kbl_user_switches

   TYPE(kbl_user_switches), SAVE :: cable_user

   ! external files read/written by CABLE
   TYPE filenames_type

   CHARACTER(LEN=99) ::                                                        &
      met,        & ! name of file for CABLE input
      out,        & ! name of file for CABLE output
      log,        & ! name of file for execution log
      restart_in, & ! name of restart file to read
      restart_out,& ! name of restart file to read
      LAI,        & ! name of file for default LAI
      type,       & ! file for default veg/soil type
      veg,        & ! file for vegetation parameters
      soil,       & ! name of file for soil parameters
      inits,      & ! name of file for initialisations
      soilIGBP      ! name of file for IGBP soil map

   END TYPE filenames_type

   TYPE(filenames_type) :: filename

   ! hydraulic_redistribution switch _soilsnow module
   LOGICAL ::                                                                  &
      redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
   
   ! hydraulic_redistribution parameters _soilsnow module
   REAL :: wiltParam=0.5, satuParam=0.8


   ! soil parameters read from file(filename%soil def. in cable.nml)
   ! & veg parameters read from file(filename%veg def. in cable.nml)
!   !---parameters, tolerances, etc. could be set in _directives.h
!jhan:cable.nml   real, parameter :: RAD_TOLS = 1.0e-2

!jhan:temporary measure. improve hiding
!   real, dimension(:,:), pointer,save :: c1, rhoch
      
END MODULE cable_common_module

