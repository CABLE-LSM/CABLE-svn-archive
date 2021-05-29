!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
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
! ==============================================================================

MODULE cable_common_module

USE cable_runtime_opts_mod ,ONLY : cable_user
USE cable_runtime_opts_mod ,ONLY : satuparam
USE cable_runtime_opts_mod ,ONLY : wiltparam

  IMPLICIT NONE

  !---allows reference to "gl"obal timestep in run (from atm_step)
  !---total number of timesteps, and processing node
  INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl
!jhan:rm this?
  LOGICAL :: L_fudge = .FALSE.

  INTEGER, SAVE :: CurYear  ! current year of multiannual run

  ! set from environment variable $HOME
  CHARACTER(LEN=200) ::                                                       &
       myhome

  ! switch to calc sil albedo using soil colour - Ticket #27
  LOGICAL :: calcsoilalbedo = .FALSE.
  !---Lestevens Sept2012
  !---CASACNP switches and cycle index
  LOGICAL, SAVE :: l_casacnp,l_laiFeedbk,l_vcmaxFeedbk
!jhan: from ESM1.5
   LOGICAL :: l_luc = .FALSE.
   LOGICAL :: l_thinforest = .FALSE.

  !---CABLE runtime switches def in this type
  TYPE kbl_internal_switches
     LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
          um_radiation = .FALSE., um_hydrology = .FALSE.
     LOGICAL :: offline = .FALSE., mk3l = .FALSE.
  END TYPE kbl_internal_switches

  ! instantiate internal switches
  TYPE(kbl_internal_switches), SAVE :: cable_runtime

  ! external files read/written by CABLE
  TYPE filenames_type

     CHARACTER(LEN=500) ::                                                        &
          met,        & ! name of file for CABLE input
          path='./',       & ! path for output and restart files for CABLE and CASA
          out,        & ! name of file for CABLE output
          log,        & ! name of file for execution log
          restart_in = ' ', & ! name of restart file to read
          restart_out,& ! name of restart file to read
          LAI,        & ! name of file for default LAI
          TYPE,       & ! file for default veg/soil type
          veg,        & ! file for vegetation parameters
          soil,       & ! name of file for soil parameters
          soilcolor,  & ! file for soil color(soilcolor_global_1x1.nc)
          inits,      & ! name of file for initialisations
          soilIGBP,   & ! name of file for IGBP soil map
          gw_elev       !name of file for gw/elevation data

  END TYPE filenames_type

  TYPE(filenames_type), SAVE :: filename

  ! hydraulic_redistribution switch _soilsnow module
  LOGICAL :: redistrb = .FALSE.  

  TYPE organic_soil_params
     !Below are the soil properties for fully organic soil

     REAL ::    &
          hyds_vec_organic  = 1.0e-4,&
          sucs_vec_organic = 10.3,   &
          clappb_organic = 2.91,     &
          ssat_vec_organic = 0.9,    &
          watr_organic   = 0.1,     &
          sfc_vec_hk      = 1.157407e-06, &
          swilt_vec_hk      = 2.31481481e-8

  END TYPE organic_soil_params

  TYPE gw_parameters_type

     REAL ::                   &
          MaxHorzDrainRate=2e-4,  & !anisintropy * q_max [qsub]
          EfoldHorzDrainRate=2.0, & !e fold rate of q_horz
          MaxSatFraction=2500.0,     & !parameter controll max sat fraction
          hkrz=0.5,               & !hyds_vec variation with z
          zdepth=1.5,             & !level where hyds_vec(z) = hyds_vec(no z)
          frozen_frac=0.05,       & !ice fraction to determine first non-frozen layer for qsub
          SoilEvapAlpha = 1.0,    & !modify field capacity dependence of soil evap limit
          IceAlpha=3.0,           &
          IceBeta=1.0

     REAL :: ice_impedence=5.0

     TYPE(organic_soil_params) :: org

     INTEGER :: level_for_satfrac = 6
     LOGICAL :: ssgw_ice_switch = .FALSE.

     LOGICAL :: subsurface_sat_drainage = .TRUE.

  END TYPE gw_parameters_type

  TYPE(gw_parameters_type), SAVE :: gw_params

  REAL, SAVE ::        &!should be able to change parameters!!!
       max_glacier_snowd=1100.0,&
       snow_ccnsw = 2.0, &
                                !jh!an:clobber - effectively force single layer snow
                                !snmin = 100.0,      & ! for 1-layer;
       snmin = 1.,          & ! for 3-layer;
       max_ssdn = 750.0,    & !
       max_sconds = 2.51,   & !
       frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)

contains 

!jhan:Start from ESM1.5
SUBROUTINE get_type_parameters(logn,vegparmnew, classification)
USE cable_params_mod,         ONLY: vegin
USE cable_params_mod,         ONLY: soilin

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
  !jhan:fudge
  CHARACTER(LEN=70) :: veg_desc(17)
  CHARACTER(LEN=70) :: soil_desc(17)

 
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
!      ALLOCATE (                                                               &
!         vegin%canst1( mvtype ), vegin%dleaf( mvtype ),                        &
!         vegin%length( mvtype ), vegin%width( mvtype ),                        &
!         vegin%vcmax( mvtype ),  vegin%ejmax( mvtype ),                        &
!         vegin%hc( mvtype ), vegin%xfang( mvtype ),                            &
!         vegin%rp20( mvtype ), vegin%rpcoef( mvtype ),                         &
!         vegin%rs20( mvtype ), vegin%wai( mvtype ),                            &
!         vegin%rootbeta( mvtype ), vegin%shelrb( mvtype ),                     &
!         vegin%vegcf( mvtype ), vegin%frac4( mvtype ),                         &
!         vegin%xalbnir( mvtype ), vegin%extkn( mvtype ),                       &
!         vegin%tminvj( mvtype ), vegin%tmaxvj( mvtype ),                       &
!         vegin%vbeta( mvtype ), vegin%froot( ms, mvtype ),                     &
!         vegin%cplant( ncp, mvtype ), vegin%csoil( ncs, mvtype ),              &
!         vegin%ratecp( ncp, mvtype ), vegin%ratecs( ncs, mvtype ),             &
!         vegin%refl( nrb, mvtype ), vegin%taul( nrb, mvtype ),             &
!         veg_desc( mvtype ) )
      
      
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
  
      !ALLOCATE ( soil_desc(mstype) )
      !ALLOCATE ( soilin%silt(mstype), soilin%clay(mstype), soilin%sand(mstype) )
      !ALLOCATE ( soilin%swilt(mstype), soilin%sfc(mstype), soilin%ssat(mstype) )
      !ALLOCATE ( soilin%bch(mstype), soilin%hyds(mstype), soilin%sucs(mstype) )
      !ALLOCATE ( soilin%rhosoil(mstype), soilin%css(mstype) )
     
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

!! get svn revision number and status
!SUBROUTINE report_version_no( logn )
!   INTEGER, INTENT(IN) :: logn
!   ! set from environment variable $HOME
!   CHARACTER(LEN=200) ::                                                       & 
!      myhome,       & ! $HOME (POSIX) environment/shell variable
!      fcablerev,    & ! recorded svn revision number at build time
!      icable_status   ! recorded svn STATUS at build time (ONLY 200 chars of it)
!
!   
!   INTEGER :: icable_rev, ioerror
!    
!   CALL getenv("HOME", myhome) 
!   fcablerev = TRIM(myhome)//TRIM("/.cable_rev")
!   
!   OPEN(440,FILE=TRIM(fcablerev),STATUS='old',ACTION='READ',IOSTAT=ioerror)
!
!      IF(ioerror==0) then 
!         ! get svn revision number (see WRITE comments)
!         READ(440,*) icable_rev
!      ELSE 
!         icable_rev=0 !default initialization
!         PRINT *, "We'll keep running but the generated revision number "     
!         PRINT *, " in the log & file will be meaningless."     
!      ENDIF
!      
!      
!      WRITE(logn,*) ''
!      WRITE(logn,*) 'Revision nuber: ', icable_rev
!      WRITE(logn,*) ''
!      WRITE(logn,*)'This is the latest revision of you workin copy as sourced ' 
!      WRITE(logn,*)'by the SVN INFO command at build time. Please note that the' 
!      WRITE(logn,*)'accuracy of this number is dependent on how recently you ' 
!      WRITE(logn,*)'used SVN UPDATE.'
!   
!      ! get svn status (see WRITE comments)
!      ! (jhan: make this output prettier & not limitted to 200 chars) 
!      WRITE(logn,*)'SVN STATUS indicates that you have (at least) the following'
!      WRITE(logn,*)'local changes: '
!      IF(ioerror==0) then 
!         READ(440,'(A)',IOSTAT=ioerror) icable_status
!         WRITE(logn,*) TRIM(icable_status)
!         WRITE(logn,*) ''
!      else   
!         WRITE(logn,*) '.cable_rev file does not exist,' 
!         WRITE(logn,*) 'suggesting you did not build libcable here' 
!         WRITE(logn,*) ''
!      endif 
!
!   CLOSE(440)
!
!END SUBROUTINE report_version_no

!jhan:End   from ESM1.5

  ELEMENTAL FUNCTION IS_LEAPYEAR( YYYY )
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: YYYY
    LOGICAL :: IS_LEAPYEAR

    IS_LEAPYEAR = .FALSE.
    IF ( ( ( MOD( YYYY,  4 ) .EQ. 0 .AND. MOD( YYYY, 100 ) .NE. 0 ) .OR. &
         MOD( YYYY,400 ) .EQ. 0 ) ) IS_LEAPYEAR = .TRUE.

  END FUNCTION IS_LEAPYEAR

  FUNCTION LEAP_DAY( YYYY )
    IMPLICIT NONE
    INTEGER :: YYYY, LEAP_DAY

    IF ( IS_LEAPYEAR ( YYYY ) ) THEN
       LEAP_DAY = 1
    ELSE
       LEAP_DAY = 0
    END IF
  END FUNCTION LEAP_DAY



END MODULE cable_common_module
