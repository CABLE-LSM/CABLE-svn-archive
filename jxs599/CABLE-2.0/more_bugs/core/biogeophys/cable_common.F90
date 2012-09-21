!========================================================================!
!=== PURPOSE: to enable accessibility of fundamental vars in global   ===!
!=== model thru out program (mainly CABLE).                           ===!
!=== USE: use module in subroutines (SRs) at top level of global model===!
!--- 2 define cable_timestep_ data with vars peculiar to global model ===!
!=== optionally call alias_ SR to give name smore familiar to cable   ===!
!=== people. then again use this module in  any subsequent SR in which===!
!=== you want to access this data.                                    ===!  
!========================================================================!

MODULE cable_common_module
   IMPLICIT NONE 

   !---allows reference to "gl"obal timestep in run (from atm_step)
   !---total number of timesteps, and processing node 
   INTEGER, SAVE :: ktau_gl, kend_gl, knode_gl, kwidth_gl
   
   !---CABLE runtime switches def in this type
   TYPE kbl_internal_switches
      LOGICAL :: um = .FALSE., um_explicit = .FALSE., um_implicit = .FALSE.,   &
            um_radiation = .FALSE., um_hydrology = .FALSE.
      LOGICAL :: offline = .FALSE., mk3l = .FALSE.
   END TYPE kbl_internal_switches 

   TYPE(kbl_internal_switches), SAVE :: cable_runtime

   !---CABLE runtime switches def in this type
   TYPE kbl_user_switches
      
      CHARACTER(LEN=200) ::                                                    &
         VEG_PARS_FILE,       & ! 
         LEAF_RESPIRATION,    & !
         FWSOIL_SWITCH          !
      
      CHARACTER(LEN=20) :: DIAG_SOIL_RESP !
      CHARACTER(LEN=5) :: RUN_DIAG_LEVEL  !
      CHARACTER(LEN=3) :: SSNOW_POTEV     !
      LOGICAL ::                                                               &
         INITIALIZE_MAPPING = .FALSE., & ! 
         CONSISTENCY_CHECK = .FALSE.,  & !
         CASA_DUMP_READ = .FALSE.,     & !
         CASA_DUMP_WRITE = .FALSE.,    & !
         CABLE_RUNTIME_COUPLED  = .FALSE.!


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
      redistrb      ! Turn on/off the hydraulic redistribution
   
   ! hydraulic_redistribution parameters _soilsnow module
   REAL :: wiltParam, satuParam


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
      
   ! new calculation dleaf sin ce April 2012 
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


END MODULE cable_common_module

