module cable_driver_module
   
   implicit none
   

contains

subroutine cable_offline_driver( met, air, canopy, rad, rough, &
                                 ssnow, soil, veg, bal,        &
                                 sum_flux, bgc, &
                                 casabiome, casapool, casaflux, &
                                 casamet, casabal, phen, dup, buf ) 
                         
   use cable_driverData_mod
   use cable_def_types_mod, only : alloc_cbm_var
   
   use cable_DUPlicate_types_mod, only : alloc_DupVars, set_Dupvars, Reset_Dupvars
   use cable_DUPlicate_types_mod, only : TDupVars
   
   !use cable_buffer_types_mod, only : alloc_BufVars, set_Bufvars, test_Bufvars
   use cable_Buffer_types_mod, only : TBufVars

   ! CABLE variables
   TYPE (met_type)       :: met     ! met input variables
   TYPE (air_type)       :: air     ! air property variables
   TYPE (canopy_type)    :: canopy  ! vegetation variables
   TYPE (radiation_type) :: rad     ! radiation variables
   TYPE (roughness_type) :: rough   ! roughness varibles
   TYPE (balances_type)  :: bal     ! energy and water balance variables
   TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
   
   ! CABLE parameters
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
   TYPE (driver_type)    :: C         ! constants used locally  
   
   TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
   
   ! CASA-CNP variables 
   TYPE (casa_biome)     :: casabiome
   TYPE (casa_pool)      :: casapool
   TYPE (casa_flux)      :: casaflux
   TYPE (casa_met)       :: casamet
   TYPE (casa_balance)   :: casabal
   TYPE (phen_variable)  :: phen 
   
   type (TdupVars) :: dup 
   type (TbufVars) :: buf 
   
   logical, save :: first_call = .true.
   integer, save :: iDiag0 
   
   ! Vars for standard for quasi-bitwise reproducability b/n runs
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   CHARACTER(len=30), PARAMETER ::                                             &
      Ftrunk_sumbal  = ".trunk_sumbal",                                        &
      Fnew_sumbal    = "new_sumbal"

   DOUBLE PRECISION ::                                                                     &
      trunk_sumbal = 0.0, & !
      new_sumbal = 0.0
      integer :: ioerror

   ! END header

   ! Open, read and close the namelist file.
   OPEN( 10, FILE = CABLE_NAMELIST )
      READ( 10, NML=CABLE )   !where NML=CABLE defined above
   CLOSE(10)
   CALL CABLE_error_log( "Namelist finished read (driver)" )

   ! Open, read and close the consistency check file.
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   IF(cable_user%consistency_check) THEN 
      OPEN( 11, FILE = Ftrunk_sumbal,STATUS='old',ACTION='READ',IOSTAT=ioerror )
         IF(ioerror==0) then
            READ( 11, * ) trunk_sumbal  ! written by previous trunk version
         ENDIF
      CLOSE(11)
   ENDIF
   
   ! Open log file:
   OPEN(logn,FILE=filename%log)
 
   CALL report_version_no( logn )
    
   IF( IARGC() > 0 ) THEN
      CALL GETARG(1, filename%met)
      CALL GETARG(2, casafile%cnpipool)
      CALL CABLE_error_log( " CLI args to CABLE: ") 
   ENDIF

   cable_runtime%offline = .TRUE.
   
   ! associate pointers used locally with global definitions
   CALL point2constants( C )
   
   ! StopGaps 
   IF( l_casacnp  .AND. ( icycle == 0 .OR. icycle > 3 ) )                   &
      STOP 'icycle must be 1 to 3 when using casaCNP'
   IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) .AND. ( .NOT. l_casacnp ) )       &
      STOP 'casaCNP required to get prognostic LAI or Vcmax'
   IF( l_vcmaxFeedbk .AND. icycle < 2 )                                     &
      STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
   IF( icycle > 0 .AND. ( .NOT. soilparmnew ) )                             &
      STOP 'casaCNP must use new soil parameters'

   ! Check for gswp run
   IF (ncciy /= 0) THEN
      PRINT *, 'Looking for global offline run info.'
      IF (ncciy < 1986 .OR. ncciy > 1995) THEN
         PRINT *, 'Year ', ncciy, ' outside range of dataset!'
         STOP 'Please check input in namelist file.'
      ELSE
         CALL prepareFiles(ncciy)
      ENDIF
   ENDIF
   

   ! Open met data and get site information from netcdf file.
   ! This retrieves time step size, number of timesteps, starting date,
   ! latitudes, longitudes, number of sites. 
   CALL open_met_file( dels, kend, spinup, C%TFRZ )
 
   ! Checks where parameters and initialisations should be loaded from.
   ! If they can be found in either the met file or restart file, they will 
   ! load from there, with the met file taking precedence. Otherwise, they'll
   ! be chosen from a coarse global grid of veg and soil types, based on 
   ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
   CALL load_parameters( met, air, ssnow, veg, bgc,                            &
                         soil, canopy, rough, rad, sum_flux,                   &
                         bal, logn, vegparmnew, casabiome, casapool,           &
                         casaflux, casamet, casabal, phen, C%EMSOIL,        &
                         C%TFRZ )
!JHAN: allocate duplicate vars as mp, kend                          
   if( first_call ) then 
      call alloc_DupVars(dup, kend)
      ! can add another copy of vars to Buf vars
      !call alloc_BufVars(buf, kend)
      first_call = .false.
   endif
            
   ! Open output file:
   CALL open_output_file( dels, soil, veg, bgc, rough )
   CALL CABLE_error_log( "output file opened" )
 
   ssnow%otss_0 = ssnow%tgg(:,1)
   ssnow%otss = ssnow%tgg(:,1)
   canopy%fes_cor = 0.
   canopy%fhs_cor = 0.
   met%ofsd = 0.1
   
   ! outer loop - spinup loop no. ktau_tot :
   ktau_tot = 0 
   
!JHAN: take out main DO loop 
!@   DO

      ! globally (WRT code) accessible kend through USE cable_common_module
      ktau_gl = 0
      kend_gl = kend
      knode_gl = 0
!JHAN: original code, timestep starts here   
      ! time step loop over ktau
      DO ktau=kstart, kend 
         
         ! increment total timstep counter
         ktau_tot = ktau_tot + 1
         
         ! globally (WRT code) accessible kend through USE cable_common_module
         ktau_gl = ktau_gl + 1
         
         ! somethings (e.g. CASA-CNP) only need to be done once per day  
         ktauday=int(24.0*3600.0/dels)
         idoy = mod(ktau/ktauday,365)
         IF(idoy==0) idoy=365
         
         ! needed for CASA-CNP
         nyear =INT((kend-kstart+1)/(365*ktauday))
   
         canopy%oldcansto=canopy%cansto
   
         ! Get met data and LAI, set time variables.
         ! Rainfall input may be augmented for spinup purposes:
         met%ofsd = met%fsd(:,1) + met%fsd(:,2) 
         CALL get_met_data( spinup, spinConv, met, soil,                    &
                            rad, veg, kend, dels, C%TFRZ, ktau ) 

         ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
         IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                                                casapool, casamet )
   
         IF (l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

#ifdef PROJECT 
!JHAN: NEW CODE BLOCK STARTS HERE: this is where _cbm is called originally
!JHAN: make a copy of met data per ktau as they are read in and 
!JHAN: store in arrays with ktau element
         call set_DupVars( ktau, dup, met, bgc, soil, veg, ssnow, canopy, &
                        rad, rough, air )

         !print *,"jhan:STOP"
         !STOP
!JHAN: below re-sets these vars to NaN. "!!" therefore suggests they cannot be reset safely
!JHAN: and need to be recorded at ktau level
         CALL alloc_cbm_var(air,    mp)
         !CALL alloc_cbm_var(bgc,   mp)
         CALL alloc_cbm_var(canopy,mp)
         CALL alloc_cbm_var(met,   mp)
         !CALL alloc_cbm_var(bal,   mp)
         CALL alloc_cbm_var(rad,   mp)
         CALL alloc_cbm_var(rough, mp)
         !CALL alloc_cbm_var(soil,  mp)
         !CALL alloc_cbm_var(ssnow, mp)
         !CALL alloc_cbm_var(veg,   mp)
         !CALL alloc_cbm_var(sum_flux, mp)

!JHAN: THIS IS COMMENTED OLD CODE in attempt to do only once
!@         ! CALL land surface scheme for this timestep, all grid points:
!@         CALL cbm( dels, air, bgc, canopy, met,                             &
!@                   bal, rad, rough, soil, ssnow,                            &
!@                   sum_flux, veg )
!@         
!@         ! resets met data per ktau from the copies made previously 
!@         ! store in arrays with ktau element
!@         call reset_DupVars( ktau, dup, met, bgc, soil,veg, ssnow, canopy, & 
!@                        rad, rough, air )
!@   
!@         ! CALL land surface scheme for this timestep, all grid points:
!@         CALL cbm( dels, air, bgc, canopy, met,                             &
!@                   bal, rad, rough, soil, ssnow,                            &
!@                   sum_flux, veg )
!@         
!@         CALL CABLE_error_log( "cbm finished" )
!@         
!@         ssnow%smelt = ssnow%smelt*dels
!@         ssnow%rnof1 = ssnow%rnof1*dels
!@         ssnow%rnof2 = ssnow%rnof2*dels
!@         ssnow%runoff = ssnow%runoff*dels
!@   
!@   
!@         !jhan this is insufficient testing. condition for 
!@         !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
!@         IF(icycle >0) THEN
!@            call bgcdriver( ktau, kstart, kend, dels, met,                     &
!@                            ssnow, canopy, veg, soil, casabiome,               &
!@                            casapool, casaflux, casamet, casabal,              &
!@                            phen, spinConv, spinup, ktauday, idoy,             &
!@                            .FALSE., .FALSE. )
!@         ENDIF 
!@   
!@         ! sumcflux is pulled out of subroutine cbm
!@         ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
!@         CALL sumcflux( ktau, kstart, kend, dels, bgc,                         &
!@                        canopy, soil, ssnow, sum_flux, veg,                    &
!@                        met, casaflux, l_vcmaxFeedbk )
!@   
!@         ! Write time step's output to file if either: we're not spinning up 
!@         ! or we're spinning up and the spinup has converged:
!@         IF((.NOT.spinup).OR.(spinup.AND.spinConv))                            &
!@            CALL write_output( dels, ktau, met, canopy, ssnow,                 &
!@                               rad, bal, air, soil, veg, C%SBOLTZ,             &
!@                               C%EMLEAF, C%EMSOIL )
!@   
!@         ! dump bitwise reproducible testing data
!@         IF( cable_user%RUN_DIAG_LEVEL == 'zero') THEN
!@            IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
!@               call cable_diag( iDiag0, "FLUXES", mp, kend, ktau,              &
!@                                knode_gl, "FLUXES",                            &
!@                          canopy%fe + canopy%fh )
!@         ENDIF
!@         
!JHAN: END 	the timestep here - shows we cantread data and then call cbm in new loop - WHY?
      END DO ! END Do loop over timestep ktau


!JHAN: times set before DO loop over timestep
      ktau_tot = 0 
      ktau_gl = 0
      kend_gl = kend
      knode_gl = 0


!JHAN: restart timestep loopp
!@      ! time step loop over ktau
      DO ktau=kstart, kend 
         
         ! increment total timstep counter
         ktau_tot = ktau_tot + 1
         
         ! globally (WRT code) accessible kend through USE cable_common_module
         ktau_gl = ktau_gl + 1
         
         ! somethings (e.g. CASA-CNP) only need to be done once per day  
         ktauday=int(24.0*3600.0/dels)
         idoy = mod(ktau/ktauday,365)
         IF(idoy==0) idoy=365
         
         ! needed for CASA-CNP
         nyear =INT((kend-kstart+1)/(365*ktauday))
   
!@         canopy%oldcansto=canopy%cansto

!JHAN: resets met data per ktau from the copies made previously 
!JHAN: store in arrays with ktau element
         call reset_DupVars( ktau, dup, met, bgc, soil,veg, ssnow, canopy, & 
                        rad, rough, air ) 
!JHAN: this is ssecond time this is called
          ! Get met data and LAI, set time variables.
!         ! Rainfall input may be augmented for spinup purposes:
!!@          met%ofsd = met%fsd(:,1) + met%fsd(:,2) 
!!@         CALL get_met_data( spinup, spinConv, met, soil,                    &
!!@                            rad, veg, kend, dels, C%TFRZ, ktau ) 
!!@
         ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
         IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                                                casapool, casamet )
   
         IF (l_laiFeedbk) veg%vlai(:) = casamet%glai(:)
!JHAN: NEW CODE BLOCK ENDS HERE
#endif 

         ! CALL land surface scheme for this timestep, all grid points:
         CALL cbm( dels, air, bgc, canopy, met,                             &
                   bal, rad, rough, soil, ssnow,                            &
                   sum_flux, veg )
         
         !call Set_BufVars( ktau, buf, met, bgc, soil,veg, ssnow, canopy, &
         !               rad, rough, air )

         CALL CABLE_error_log( "cbm finished" )
         

         !call test_BufVars( ktau, buf, met, bgc, soil,veg, ssnow, canopy, &
         !               rad, rough, air )


         ssnow%smelt = ssnow%smelt*dels
         ssnow%rnof1 = ssnow%rnof1*dels
         ssnow%rnof2 = ssnow%rnof2*dels
         ssnow%runoff = ssnow%runoff*dels
   
   
         !jhan this is insufficient testing. condition for 
         !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
         IF(icycle >0) THEN
            call bgcdriver( ktau, kstart, kend, dels, met,                     &
                            ssnow, canopy, veg, soil, casabiome,               &
                            casapool, casaflux, casamet, casabal,              &
                            phen, spinConv, spinup, ktauday, idoy,             &
                            .FALSE., .FALSE. )
         ENDIF 
   
         ! sumcflux is pulled out of subroutine cbm
         ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
         CALL sumcflux( ktau, kstart, kend, dels, bgc,                         &
                        canopy, soil, ssnow, sum_flux, veg,                    &
                        met, casaflux, l_vcmaxFeedbk )
   
         ! Write time step's output to file if either: we're not spinning up 
         ! or we're spinning up and the spinup has converged:
         IF((.NOT.spinup).OR.(spinup.AND.spinConv))                            &
            CALL write_output( dels, ktau, met, canopy, ssnow,                 &
                               rad, bal, air, soil, veg, C%SBOLTZ,             &
                               C%EMLEAF, C%EMSOIL )
   
         ! dump bitwise reproducible testing data
         IF( cable_user%RUN_DIAG_LEVEL == 'zero') THEN
            IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
               call cable_diag( iDiag0, "FLUXES", mp, kend, ktau,              &
                                knode_gl, "FLUXES",                            &
                          canopy%fe + canopy%fh )
         ENDIF
         
      END DO ! END Do loop over timestep ktau

     
   ! Close met data input file:
   CALL close_met_file
 
   ! Close output file and deallocate main variables:
   CALL close_output_file( bal, air, bgc, canopy, met,                         &
                           rad, rough, soil, ssnow,                            &
                           sum_flux, veg )

  
   ! Check this run against standard for quasi-bitwise reproducability
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   IF(cable_user%consistency_check) THEN 
      
      new_sumbal = SUM(bal%wbal_tot) + SUM(bal%ebal_tot)                       &
                       + SUM(bal%ebal_tot_cncheck)
  
      IF( new_sumbal == trunk_sumbal) THEN

         print *, ""
         print *, &
         "Internal check shows this version reproduces the trunk sumbal"
      
      ELSE

         print *, ""
         print *, &
         "Internal check shows in this version new_sumbal != trunk sumbal"
         print *, &
         "Writing new_sumbal to the file:", TRIM(Fnew_sumbal)
               
         OPEN( 12, FILE = Fnew_sumbal )
            WRITE( 12, * ) new_sumbal  ! written by previous trunk version
         CLOSE(12)
      
      ENDIF   
      
   ENDIF

   ! Close log file
   CLOSE(logn)

END subroutine cable_offline_driver


SUBROUTINE prepareFiles(ncciy)
  USE cable_IO_vars_module, ONLY: logn,gswpfile
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncciy

  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy

  CALL renameFiles(logn,gswpfile%rainf,16,ncciy,'rainf')
  CALL renameFiles(logn,gswpfile%snowf,16,ncciy,'snowf')
  CALL renameFiles(logn,gswpfile%LWdown,16,ncciy,'LWdown')
  CALL renameFiles(logn,gswpfile%SWdown,16,ncciy,'SWdown')
  CALL renameFiles(logn,gswpfile%PSurf,16,ncciy,'PSurf')
  CALL renameFiles(logn,gswpfile%Qair,14,ncciy,'Qair')
  CALL renameFiles(logn,gswpfile%Tair,14,ncciy,'Tair')
  CALL renameFiles(logn,gswpfile%wind,15,ncciy,'wind')

END SUBROUTINE prepareFiles


SUBROUTINE renameFiles(logn,inFile,nn,ncciy,inName)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: logn
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: ncciy
  CHARACTER(LEN=99), INTENT(INOUT) :: inFile
  CHARACTER(LEN=*),  INTENT(IN)    :: inName
  INTEGER :: idummy

  READ(inFile(nn:nn+3),'(i4)') idummy
  IF (idummy < 1983 .OR. idummy > 1995) THEN
    PRINT *, 'Check position of the year number in input gswp file', inFile
    STOP
  ELSE
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
  ENDIF

END SUBROUTINE renameFiles


end module cable_driver_module




