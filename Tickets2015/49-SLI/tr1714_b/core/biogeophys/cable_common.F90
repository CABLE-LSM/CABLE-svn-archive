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
   INTEGER, SAVE :: CurYear  ! current year of multiannual run
   
   ! user switches turned on/off by the user thru namelists

   ! trunk modifications protected by these switches 
   TYPE hide_switches
      LOGICAL ::                                                               & 
         ! L.Stevens - Test Switches 
         L_NEW_ROUGHNESS_SOIL  = .FALSE., & ! from Ticket? 
         L_NEW_RUNOFF_SPEED    = .FALSE., & ! from Ticket?
         L_NEW_REDUCE_SOILEVP  = .FALSE., & ! from Ticket?
         Ticket46 = .FALSE.,              & !
         !jhan: default should be FALSE, bu set up nml etc
         Ticket49Bug1 = .true.,           & ! 
         Ticket49Bug2 = .true.,           & ! 
         Ticket49Bug3 = .true.,           & ! 
         Ticket49Bug4 = .true.,           & ! 
         Ticket49Bug5 = .true.,           & ! 
         Ticket49Bug6 = .true.              ! 

      END TYPE hide_switches 

   ! instantiate internal switches 
   TYPE (hide_switches), SAVE :: hide
   
   
   ! set from environment variable $HOME
   CHARACTER(LEN=200) ::                                                       & 
      myhome

   ! switch to calc sil albedo using soil colour - Ticket #27
   LOGICAL :: calcsoilalbedo = .FALSE. 
   !---Lestevens Sept2012
   !---CASACNP switches and cycle index
   LOGICAL, SAVE :: l_casacnp,l_laiFeedbk,l_vcmaxFeedbk
   
   !---CABLE runtime switches def in this type
   TYPE kbl_internal_switches
      LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
            um_radiation = .FALSE.
      LOGICAL :: offline = .FALSE., mk3l = .FALSE.
   END TYPE kbl_internal_switches 

   ! instantiate internal switches 
   TYPE(kbl_internal_switches), SAVE :: cable_runtime

   ! user switches turned on/off by the user thru namelists
   ! CABLE-2.0 user switches all in single namelist file cable.nml
   ! clean these up for new namelist(s) format	
   TYPE kbl_user_switches
      !jhan: this is redundant now we all use filename%veg?
      CHARACTER(LEN=200) ::                                                    &
         VEG_PARS_FILE  ! 
      
      CHARACTER(LEN=20) ::                                                     &
         FWSOIL_SWITCH     !
      
     CHARACTER(LEN=10):: RunIden  !
     CHARACTER(LEN=4) :: MetType  !
     CHARACTER(LEN=20) :: CANOPY_STRUC !
     CHARACTER(LEN=20) :: SOIL_STRUC !
     CHARACTER(LEN=3)  :: POP_out = 'rst' ! POP output type ('epi' or 'rst')
     CHARACTER(LEN=50) :: POP_rst = ' ' !

   LOGICAL ::                                                               &
          CALL_POP               = .FALSE., & !
          POP_fromZero           = .FALSE.

     INTEGER  :: &
          CASA_SPIN_STARTYEAR = 1950, &
          CASA_SPIN_ENDYEAR   = 1960, &
          YEARSTART           = 1950, &
          YEAREND             = 1960, &
          CASA_OUT_FREQ       = 365, &
          CASA_NREP           = 1

      CHARACTER(LEN=5) ::                                                      &
         RUN_DIAG_LEVEL  !
      
      CHARACTER(LEN=3) ::                                                      &
         SSNOW_POTEV,      & !
         DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical) 
         LEAF_RESPIRATION    ! either ON or OFF (jhan:Make Logical) 

      ! Custom soil respiration - see Ticket #42
      CHARACTER(LEN=10) ::                                                     &
         SMRF_NAME,   & ! Soil Moist Respiration Function
         STRF_NAME      ! Soil Temp Respiration Function

      LOGICAL ::                                                               &
         INITIALIZE_MAPPING = .FALSE., & ! 
         CONSISTENCY_CHECK = .FALSE.,  & !
         CASA_DUMP_READ = .FALSE.,     & !
         CASA_DUMP_WRITE = .FALSE.,    & !
         CABLE_RUNTIME_COUPLED = .TRUE., & !
         ! L.Stevens - Test Switches
         L_NEW_ROUGHNESS_SOIL  = .FALSE., & !
         L_NEW_RUNOFF_SPEED    = .FALSE., & !
         L_NEW_REDUCE_SOILEVP  = .FALSE., & !

	     ! Switch for customized soil respiration - see Ticket #42
         SRF = .FALSE.
         
   END TYPE kbl_user_switches

   ! instantiate internal switches 
   TYPE(kbl_user_switches), SAVE :: cable_user

   ! external files read/written by CABLE
   TYPE filenames_type

   CHARACTER(LEN=200) ::                                                        &
      met,        & ! name of file for CABLE input
      path,       & ! path for output and restart files for CABLE and CASA
      out,        & ! name of file for CABLE output
      log,        & ! name of file for execution log
      restart_in, & ! name of restart file to read
      restart_out,& ! name of restart file to read
      LAI,        & ! name of file for default LAI
      type,       & ! file for default veg/soil type
      veg,        & ! file for vegetation parameters
      soil,       & ! name of file for soil parameters
      soilcolor,  & ! file for soil color(soilcolor_global_1x1.nc)
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

  SUBROUTINE HANDLE_ERR( status )
    use netcdf
    INTEGER, INTENT(IN) :: status
    IF(status /= NF90_noerr) THEN
       PRINT*,"netCDF error:"
       PRINT*, TRIM(NF90_strerror(status))
       STOP "Stopped"
    END IF
  END SUBROUTINE HANDLE_ERR

  FUNCTION IS_LEAPYEAR( YYYY )
    IMPLICIT NONE
    INTEGER :: YYYY
    LOGICAL :: IS_LEAPYEAR

    IS_LEAPYEAR = .FALSE.
    IF ( ( ( MOD( YYYY,  4 ) .EQ. 0 .AND. MOD( YYYY, 100 ) .NE. 0 ) .OR. &
         MOD( YYYY,400 ) .EQ. 0 ) ) IS_LEAPYEAR = .TRUE.

  END FUNCTION IS_LEAPYEAR

  SUBROUTINE YMDHMS2DOYSOD( YYYY,MM,DD,HOUR,MINUTE,SECOND,DOY,SOD )

    ! Compute Day-of-year and second-of-day from given date and time or
    ! reverse (if REV=.TRUE.)

    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: YYYY,MM,DD,HOUR,MINUTE,SECOND
    INTEGER,INTENT(OUT) :: DOY,SOD

    !  LOGICAL :: IS_LEAPYEAR
    INTEGER, DIMENSION(12) :: MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    IF ( IS_LEAPYEAR( YYYY ) ) MONTH(2) = 29

    IF ( DD .GT. MONTH(MM) .OR. DD .LT. 1 .OR. &
         MM .GT. 12 .OR. MM .LT. 1 ) THEN
       WRITE(*,*)"Wrong date entered in YMDHMS2DOYSOD "
       WRITE(*,*)"DATE : ",YYYY,MM,DD
       STOP
    ENDIF
    DOY = DD
    IF ( MM .GT. 1 ) DOY = DOY + SUM( MONTH( 1:MM-1 ) )
    SOD = HOUR * 3600 + MINUTE * 60 + SECOND

  END SUBROUTINE YMDHMS2DOYSOD

  SUBROUTINE DOYSOD2YMDHMS( YYYY,DOY,SOD,MM,DD,HOUR,MINUTE,SECOND )

    ! Compute Day-of-year and second-of-day from given date and time or
    ! reverse (if REV=.TRUE.)

    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: YYYY,DOY,SOD
    INTEGER,INTENT(OUT) :: MM,DD,HOUR,MINUTE,SECOND

    !  LOGICAL :: IS_LEAPYEAR
    INTEGER :: MON, i
    INTEGER, DIMENSION(12) :: MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    IF ( IS_LEAPYEAR( YYYY ) ) MONTH(2) = 29

    IF ( SOD .GE. 86400 .OR. SOD .LT. 0 .OR. &
         DOY .GT. SUM(MONTH) .OR. DOY .LT. 1 ) THEN
       WRITE(*,*)"Wrong date entered in DOYSOD2YMDHMS "
       WRITE(*,*)"DOYSOD : ",DOY,SOD
       STOP
    ENDIF

    MON = 0
    DO i = 1, 12
       IF ( MON + MONTH(i) .LT. DOY ) THEN
          MON = MON + MONTH(i)
       ELSE
          MM  = i
          DD  = DOY - MON
          EXIT
       ENDIF
    END DO
    HOUR   = INT( REAL(SOD)/3600. )
    MINUTE = INT( ( REAL(SOD) - REAL(HOUR)*3600.) / 60. )
    SECOND = SOD - HOUR*3600 - MINUTE*60

  END SUBROUTINE DOYSOD2YMDHMS

! get svn revision number and status
SUBROUTINE report_version_no( logn )

#ifdef NAG
   USE F90_UNIX_ENV, only: getenv
#endif
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

