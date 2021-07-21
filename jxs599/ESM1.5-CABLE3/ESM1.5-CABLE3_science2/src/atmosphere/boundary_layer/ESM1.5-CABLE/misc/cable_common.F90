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
! Purpose: *USED to* Read vegetation and soil parameter files, fills vegin, soilin
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
USE cable_runtime_opts_mod ,ONLY: cable_user
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


   CHARACTER(LEN=70), DIMENSION(:), POINTER ::                                 &
      veg_desc,   & ! decriptions of veg type
      soil_desc     ! decriptns of soil type 

!   !---parameters, tolerances, etc. could be set in _directives.h
!jhan:cable.nml   real, parameter :: RAD_TOLS = 1.0e-2

!jhan:temporary measure. improve hiding
!   real, dimension(:,:), pointer,save :: c1, rhoch
      
  REAL, SAVE ::        &!should be able to change parameters!!!
       max_glacier_snowd=1100.0,&
       snow_ccnsw = 2.0, &
                                !jh!an:clobber - effectively force single layer snow
                                !snmin = 100.0,      & ! for 1-layer;
       snmin = 1.,          & ! for 3-layer;
       max_ssdn = 750.0,    & !
       max_sconds = 2.51,   & !
       frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)

CONTAINS


! get svn revision number and status
SUBROUTINE report_version_no( logn )
   INTEGER, INTENT(IN) :: logn
   ! set from environment variable $HOME
   CHARACTER(LEN=200) ::                                                       & 
      myhome,       & ! $HOME (POSIX) environment/shell variable
      fcablerev,    & ! recorded svn revision number at build time
      icable_status   ! recorded svn STATUS at build time (ONLY 200 chars of it)

   
   INTEGER :: icable_rev, ioerror
    
   CALL getenv("HOME", myhome) 
   fcablerev = TRIM(myhome)//TRIM("/.cable_rev")
   
   OPEN(440,FILE=TRIM(fcablerev),STATUS='old',ACTION='READ',IOSTAT=ioerror)

      IF(ioerror==0) then 
         ! get svn revision number (see WRITE comments)
         READ(440,*) icable_rev
      ELSE 
         icable_rev=0 !default initialization
         PRINT *, "We'll keep running but the generated revision number "     
         PRINT *, " in the log & file will be meaningless."     
      ENDIF
      
      
      WRITE(logn,*) ''
      WRITE(logn,*) 'Revision nuber: ', icable_rev
      WRITE(logn,*) ''
      WRITE(logn,*)'This is the latest revision of you workin copy as sourced ' 
      WRITE(logn,*)'by the SVN INFO command at build time. Please note that the' 
      WRITE(logn,*)'accuracy of this number is dependent on how recently you ' 
      WRITE(logn,*)'used SVN UPDATE.'
   
      ! get svn status (see WRITE comments)
      ! (jhan: make this output prettier & not limitted to 200 chars) 
      WRITE(logn,*)'SVN STATUS indicates that you have (at least) the following'
      WRITE(logn,*)'local changes: '
      IF(ioerror==0) then 
         READ(440,'(A)',IOSTAT=ioerror) icable_status
         WRITE(logn,*) TRIM(icable_status)
         WRITE(logn,*) ''
      else   
         WRITE(logn,*) '.cable_rev file does not exist,' 
         WRITE(logn,*) 'suggesting you did not build libcable here' 
         WRITE(logn,*) ''
      endif 

   CLOSE(440)

END SUBROUTINE report_version_no



END MODULE cable_common_module

