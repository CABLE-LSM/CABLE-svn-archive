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
  !jhan:make this logical
  CHARACTER(LEN=3) :: diag_soil_resp=''

  CHARACTER(LEN=20) :: fwsoil_switch=''

  ! Ticket #56
  !jhan:options?
  CHARACTER(LEN=20) :: gs_switch=''

  !INH - new switch for revised coupling on implicit step of ACCESS-CM2 Ticket #132
  LOGICAL :: l_revised_coupling = .FALSE.

  !INH -apply revised sensitvity/correction terms to soilsnow energy balance
  LOGICAL :: l_rev_corr = .FALSE.     !switch to revert to unchanged code

  !ticket#179
  LOGICAL :: soil_thermal_fix = .FALSE.

  !jhan:options?
  CHARACTER(LEN=3) :: ssnow_potev=''
 
   !jhan: this is redundant now we all use filename%veg?
  CHARACTER(LEN=200) ::                                                       &
       veg_pars_file  !

  CHARACTER(LEN=20) ::                                                        &
       phenology_switch = 'MODIS'   ! alternative is 'climate'
  !--- LN ------------------------------------------[

  CHARACTER(LEN=10) :: RunIden       = 'STANDARD'  !
  CHARACTER(LEN=6)  :: MetType       = ' ' !
  CHARACTER(LEN=20) :: soil_struc    = "default" ! 'default' or 'sli'
  CHARACTER(LEN=3)  :: POP_out       = 'rst' ! POP output type ('epi' or 'rst')
  CHARACTER(LEN=50) :: POP_rst       = ' ' !
  CHARACTER(LEN=8)  :: casa_out_freq = 'annually' ! 'daily', 'monthly', 'annually'
  CHARACTER(LEN=10)  :: vcmax = 'standard' ! "standard" or "Walker2014"
  CHARACTER(LEN=10)  :: POPLUC_RunType = 'static' ! 'static', 'init', 'restart'

  LOGICAL ::                                                                  &
       call_pop               = .FALSE., & !
       POP_fromZero           = .FALSE.,                                      &
       CALL_Climate           = .FALSE.,                                      &
       Climate_fromZero       = .FALSE.,                                      &
       CASA_fromZero          = .FALSE.,                                      &
       popluc                 = .FALSE.

  INTEGER  ::                                                                 &
       casa_spin_startyear = 1950,                                            &
       casa_spin_endyear   = 1960,                                            &
       yearstart           = 0,                                               &
       yearend             = 0,                                               &
       casa_nrep           = 1
  !--- LN ------------------------------------------]

  CHARACTER(LEN=5) ::                                                         &
       run_diag_level  !

      CHARACTER(LEN=3) ::                                                         &
       !H!DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical)
       leaf_respiration    ! either ON or OFF (jhan:Make Logical)

  ! Custom soil respiration - see Ticket #42
  CHARACTER(LEN=10) ::                                                        &
       smrf_name,   & ! Soil Moist Respiration Function
       strf_name      ! Soil Temp Respiration Function

  LOGICAL ::                                                                  &
       initialize_mapping    = .FALSE., & !
       consistency_check     = .FALSE., & !
       casa_dump_read        = .FALSE., & !
       casa_dump_write       = .FALSE., & !
       !CBL3cable_runtime_coupled = .TRUE. , & !
       CABLE_RUNTIME_COUPLED = .FALSE., & !
       LogWorker             = .TRUE. , & ! Write Output of each worker
                             ! L.Stevens - Test Switches
       l_new_roughness_soil  = .FALSE., & !
       l_new_runoff_speed    = .FALSE., & !
       l_new_reduce_soilevp  = .FALSE., & !

                             ! Switch for customized soil respiration - see Ticket #42
       srf = .FALSE.,                                                         &

                             !! vh_js !!
       litter = .FALSE.

  !MD
  LOGICAL :: gw_model = .FALSE.
  LOGICAL :: alt_forcing = .FALSE.

  !using GSWP3 forcing?
  LOGICAL :: gswp3 = .FALSE.
  LOGICAL :: or_evap = .FALSE.
  LOGICAL :: test_new_gw = .FALSE.
  LOGICAL :: sync_nc_file = .FALSE.
  INTEGER :: max_spins = -1
  LOGICAL :: fix_access_roots = .FALSE.  !use pft dependent roots in ACCESS
  !ACCESS roots
  LOGICAL :: access13roots = .FALSE.     !switch to use ACCESS1.3 %froot

  LOGICAL :: l_limit_labile = .FALSE.    ! #237: limit Labile in spinup
  LOGICAL :: NtilesThruMetFile = .FALSE. ! #199: Specify Ntiles thru met file 

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
   TYPE soilin_type

      REAL, DIMENSION(:),ALLOCATABLE ::                                        &
         silt,    & !
         clay,    & !
         sand,    & !
         swilt,   & !
         sfc,     & !
         ssat,    & !
         bch,     & !
         hyds,    & !
         sucs,    & !
         rhosoil, & !
         css,     & !
         c3         !
   
   END TYPE soilin_type
 

   TYPE vegin_type

      REAL, DIMENSION(:),ALLOCATABLE ::                                        &
         canst1,     & !
         dleaf,      & !
         length,     & !
         width,      & !
         vcmax,      & !
         ejmax,      & !
         hc,         & !
         xfang,      & !
         rp20,       & !
         rpcoef,     & !
         rs20,       & !
         wai,        & ! 
         rootbeta,   & ! 
         shelrb,     & !
         vegcf,      & !  
         frac4,      & !
         xalbnir,    & !
         extkn,      & ! 
         tminvj,     & !
         tmaxvj,     & !
         vbeta         !
      
      REAL, DIMENSION(:,:),ALLOCATABLE ::                                      &
         froot,      & !
         cplant,     & !
         csoil,      & !
         ratecp,     & !
         ratecs,     & !
         refl,     & !
         taul        !
      
   END TYPE vegin_type

   CHARACTER(LEN=70), DIMENSION(:), POINTER ::                                 &
      veg_desc,   & ! decriptions of veg type
      soil_desc     ! decriptns of soil type 

   TYPE(soilin_type), SAVE  :: soilin
   TYPE(vegin_type),  SAVE  :: vegin

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


SUBROUTINE get_type_parameters(logn,vegparmnew, classification)

   ! Gets parameter values for each vegetation type and soil type.
   
   USE cable_def_types_mod, ONLY : mvtype, ms, ncs, ncp, mstype, nrb 

   INTEGER,INTENT(IN) :: logn     ! log file unit number
   
   CHARACTER(LEN=4), INTENT(INOUT), OPTIONAL :: classification
   
   LOGICAL,INTENT(IN)      :: vegparmnew ! new format input file 
   
   CHARACTER(LEN=80) :: comments 
   CHARACTER(LEN=10) :: vegtypetmp                   
   CHARACTER(LEN=25) :: vegnametmp                   
   
   REAL    :: notused                           
   INTEGER :: ioerror ! input error integer
   INTEGER :: a, jveg ! do loop counter


 
   !================= Read in vegetation type specifications: ============
   OPEN(40,FILE=filename%veg,STATUS='old',ACTION='READ',IOSTAT=ioerror)
      
      IF(ioerror/=0) then 
         STOP 'CABLE_log: Cannot open veg type definitions.'
      ENDIF
     
      IF (vegparmnew) THEN
         
         ! assume using IGBP/CSIRO vegetation types
         READ(40,*) comments
         READ(40,*) mvtype
         IF( present(classification) )                                         &
            WRITE(classification,'(a4)') comments(1:4)
      
      ELSE
         
         ! assume using CASA vegetation types
         !classification = 'CASA'
         READ(40,*)
         READ(40,*)
         READ(40,*) mvtype ! read # vegetation types
         READ(40,*)
         READ(40,*)
         comments = 'CASA'
      
      END IF
         
      WRITE(logn, '(A31,I3,1X,A10)') '  Number of vegetation types = ',        &
                  mvtype,TRIM(comments)
   
    
      ! Allocate memory for type-specific vegetation parameters:
      ALLOCATE (                                                               &
         vegin%canst1( mvtype ), vegin%dleaf( mvtype ),                        &
         vegin%length( mvtype ), vegin%width( mvtype ),                        &
         vegin%vcmax( mvtype ),  vegin%ejmax( mvtype ),                        &
         vegin%hc( mvtype ), vegin%xfang( mvtype ),                            &
         vegin%rp20( mvtype ), vegin%rpcoef( mvtype ),                         &
         vegin%rs20( mvtype ), vegin%wai( mvtype ),                            &
         vegin%rootbeta( mvtype ), vegin%shelrb( mvtype ),                     &
         vegin%vegcf( mvtype ), vegin%frac4( mvtype ),                         &
         vegin%xalbnir( mvtype ), vegin%extkn( mvtype ),                       &
         vegin%tminvj( mvtype ), vegin%tmaxvj( mvtype ),                       &
         vegin%vbeta( mvtype ), vegin%froot( ms, mvtype ),                     &
         vegin%cplant( ncp, mvtype ), vegin%csoil( ncs, mvtype ),              &
         vegin%ratecp( ncp, mvtype ), vegin%ratecs( ncs, mvtype ),             &
         vegin%refl( nrb, mvtype ), vegin%taul( nrb, mvtype ),             &
         veg_desc( mvtype ) )
      
      
      IF( vegparmnew ) THEN    ! added to read new format (BP dec 2007)
            
         ! Read in parameter values for each vegetation type:
         DO a = 1,mvtype 
            
            READ(40,*) jveg, vegtypetmp, vegnametmp
                 
            IF( jveg .GT. mvtype )                                             &
               STOP 'jveg out of range in parameter file'
               
            veg_desc(jveg) = vegnametmp 
               
            READ(40,*) vegin%hc(jveg), vegin%xfang(jveg), vegin%width(jveg),   &
                        &   vegin%length(jveg), vegin%frac4(jveg)
            ! only refl(1:2) and taul(1:2) used
            READ(40,*) vegin%refl(1:3,jveg) ! rhowood not used ! BP may2011
            READ(40,*) vegin%taul(1:3,jveg) ! tauwood not used ! BP may2011
            READ(40,*) notused, notused, notused, vegin%xalbnir(jveg)
            READ(40,*) notused, vegin%wai(jveg), vegin%canst1(jveg),           &
               vegin%shelrb(jveg), vegin%vegcf(jveg), vegin%extkn(jveg)
            READ(40,*) vegin%vcmax(jveg), vegin%rp20(jveg),                    &
                       vegin%rpcoef(jveg),                                     &
                       vegin%rs20(jveg)
            READ(40,*) vegin%tminvj(jveg), vegin%tmaxvj(jveg),                 &
                       vegin%vbeta(jveg), vegin%rootbeta(jveg)
            READ(40,*) vegin%cplant(1:3,jveg), vegin%csoil(1:2,jveg)
            ! rates not currently set to vary with veg type
            READ(40,*) vegin%ratecp(1:3,jveg), vegin%ratecs(1:2,jveg)

         END DO

      ELSE

         DO a = 1,mvtype 
            READ(40,'(8X,A70)') veg_desc(a) ! Read description of each veg type
         END DO

         READ(40,*); READ(40,*) 
         READ(40,*) vegin%canst1
         READ(40,*) vegin%width
         READ(40,*) vegin%length
         READ(40,*) vegin%vcmax
         READ(40,*) vegin%hc
         READ(40,*) vegin%xfang
         READ(40,*) vegin%rp20
         READ(40,*) vegin%rpcoef
         READ(40,*) vegin%rs20
         READ(40,*) vegin%shelrb
         READ(40,*) vegin%frac4
         READ(40,*) vegin%wai
         READ(40,*) vegin%vegcf
         READ(40,*) vegin%extkn
         READ(40,*) vegin%tminvj
         READ(40,*) vegin%tmaxvj
         READ(40,*) vegin%vbeta
         READ(40,*) vegin%xalbnir
         READ(40,*) vegin%rootbeta
         READ(40,*) vegin%cplant(1,:)
         READ(40,*) vegin%cplant(2,:)
         READ(40,*) vegin%cplant(3,:)
         READ(40,*) vegin%csoil(1,:)
         READ(40,*) vegin%csoil(2,:)
         READ(40,*) 
         READ(40,*) vegin%ratecp(:,1)
            
         ! Set ratecp to be the same for all veg types:
         vegin%ratecp(1,:)=vegin%ratecp(1,1)
         vegin%ratecp(2,:)=vegin%ratecp(2,1)
         vegin%ratecp(3,:)=vegin%ratecp(3,1)
         READ(40,*) 
         READ(40,*) vegin%ratecs(:,1)
         vegin%ratecs(1,:)=vegin%ratecs(1,1)
         vegin%ratecs(2,:)=vegin%ratecs(2,1)
         
         ! old table does not have taul and refl ! BP may2011
         vegin%taul(1,:) = 0.07
         vegin%taul(2,:) = 0.425
         vegin%taul(3,:) = 0.0
         vegin%refl(1,:) = 0.07
         vegin%refl(2,:) = 0.425
         vegin%refl(3,:) = 0.0

      ENDIF

      WRITE(6,*)'CABLE_log:Closing veg params file: ',trim(filename%veg)
      
   CLOSE(40)
      
   ! new calculation dleaf since April 2012 (cable v1.8 did not use width)
   vegin%dleaf = SQRT(vegin%width * vegin%length)
    
        
    
  !================= Read in soil type specifications: ============
   OPEN(40,FILE=filename%soil,STATUS='old',ACTION='READ',IOSTAT=ioerror)

      IF(ioerror/=0) then 
           STOP 'CABLE_log: Cannot open soil type definitions.'
      ENDIF

      READ(40,*); READ(40,*)
      READ(40,*) mstype ! Number of soil types
      READ(40,*); READ(40,*)
  
      ALLOCATE ( soil_desc(mstype) )
      ALLOCATE ( soilin%silt(mstype), soilin%clay(mstype), soilin%sand(mstype) )
      ALLOCATE ( soilin%swilt(mstype), soilin%sfc(mstype), soilin%ssat(mstype) )
      ALLOCATE ( soilin%bch(mstype), soilin%hyds(mstype), soilin%sucs(mstype) )
      ALLOCATE ( soilin%rhosoil(mstype), soilin%css(mstype) )
     
      DO a = 1,mstype 
         READ(40,'(8X,A70)') soil_desc(a) ! Read description of each soil type
      END DO

      READ(40,*); READ(40,*) 
      READ(40,*) soilin%silt
      READ(40,*) soilin%clay
      READ(40,*) soilin%sand
      READ(40,*) soilin%swilt
      READ(40,*) soilin%sfc
      READ(40,*) soilin%ssat
      READ(40,*) soilin%bch
      READ(40,*) soilin%hyds
      READ(40,*) soilin%sucs
      READ(40,*) soilin%rhosoil
      READ(40,*) soilin%css

   CLOSE(40)

END SUBROUTINE get_type_parameters


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

