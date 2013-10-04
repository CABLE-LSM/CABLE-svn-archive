!===COPYRIGHT==================================================================
! The source codes are part of the australian 
! Community Atmosphere Biosphere Land Exchange (CABLE) model. 
! Please register online at xxx and sign the agreement before use 
! contact: whox@xxxx.yyy about registration user agreement
!==============================================================================


!==============================================================================
! Name: cable_driver
! Purpose: offline driver for CABLE model
! CALLed from: executed PROGRAM 
! MODULEs used:   define_dimensions
!                 define_types
!                 io_variables
!                 cable_common_module
!                 input_module
!                 output_module
!                 cbm_module
!                 casadimension
!                 casavariable
! 
! CALLs:       open_met_file
!              load_parameters
!              open_output_file
!              spincasacnp
!              get_met_data
!              casa_feedback
!              cbm
!              bgcdriver
!              sumcflux
!              write_output
!              casa_poolout
!              casa_fluxout
!              create_restart
!              close_met_file
!              close_output_file
!              prepareFiles
!
! Major contribution: land surface modeling team, CSIRO, Aspendale
!
! input  file: ￼[SiteName].nc
!              ￼poolcnpIn[SiteName].csv -- for CASA-CNP only
!              ￼gridinfo_CSIRO_1x1.nc
!              ￼def_veg_params.txt
!              ￼def_soil_params.txt -- nearly redundant, can be switched on
!              ￼restart_in.nc -- not strictly required
!
! output file: ￼log_cable.txt
!              ￼out_cable.nc
!              ￼restart_out.nc
!              ￼poolcnpOut.csv -- from CASA-CNP

!==============================================================================


!==============================================================================
! changes since version release on 
! changes made by who on date
!
!==============================================================================


PROGRAM offline_driver
   
   USE define_dimensions, ONLY:r_1,i_d,ms,mp,mvtype,mstype
   USE define_types
   USE io_variables, ONLY: logn,filename,gswpfile,ncciy,leaps, &
        verbose, fixedCO2,output,check,patchout,patch_type,soilparmnew, &
        redistrb, wiltParam, satuParam
   USE cable_common_module, only : ktau_gl, cable_user, cable_runtime
   USE input_module, ONLY: open_met_file,load_parameters, &
        get_met_data,close_met_file
   USE output_module, ONLY: create_restart,open_output_file, &
        write_output,close_output_file
   USE cbm_module
   ! modules related to CASA-CNP
   USE casadimension, only : icycle 
   USE casavariable, only : casafile 
   
   IMPLICIT NONE

   INTEGER(i_d)          :: kend ! no. of time steps in run
   TYPE (air_type)       :: air  ! air property variables
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
   TYPE (canopy_type)    :: canopy ! vegetation variables
   TYPE (met_type)       :: met  ! met input variables
   TYPE (balances_type)  :: bal  ! energy and water balance variables
   TYPE (radiation_type) :: rad  ! radiation variables
   TYPE (roughness_type) :: rough ! roughness varibles
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (soil_snow_type) :: ssoil ! soil and snow variables
   TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
   TYPE (veg_parameter_type) :: veg  ! vegetation parameters	 
   TYPE (casa_biome)     :: casabiome
   TYPE (casa_pool)      :: casapool
   TYPE (casa_flux)      :: casaflux
   TYPE (casa_met)       :: casamet
   TYPE (casa_balance)   :: casabal
   TYPE (phen_variable)  :: phen 
  
   LOGICAL :: &
      vegparmnew, &     ! using new format input file (BP dec 2007)
      spinup, &         ! should the model spinup to soil state equilibrium?
      spinConv, &       ! has spinup converged?
      spincasainput, &    ! TRUE: input required to spin CASA-CNP wil be saved;
                        ! FALSE: read input to spin CASA-CNP 1000 years
      spincasa, &       ! TRUE: CASA-CNP Will spin mloop times,
                        ! FALSE: no spin up
      l_casacnp, &      ! using CASA-CNP with CABLE
      l_laiFeedbk, &    ! using prognostic LAI
      l_vcmaxFeedbk     ! using prognostic Vcmax
   
   INTEGER(i_d), parameter :: & 
      kstart = 1   ! start of simulation #
   
   INTEGER(i_d) :: & 
      kstart, &   ! start of simulation #
      tstep, &    ! time step counter for spinup
      mloop, &    ! # spinup loops for CASA-CNP
      ktauday, &  !  day counter for CASA-CNP
      nyear, &    ! year counter for CASA-CNP
      ktau        !loop increment (corresponds to timestep)

   REAL(r_1) :: &  
      dels, &     ! time step size in seconds
      delsoilM, & ! allowed variation in soil moisture for spin up
      delsoilT    ! allowed variation in soil temperature for spin up
   
   REAL(r_1), ALLOCATABLE, DIMENSION(:,:)  :: & 
      soilMtemp, &   ! temporary storage for spin up
      soilTtemp      ! temporary storage for spin up

   NAMELIST/CABLE/filename, &
                  vegparmnew,&
                  soilparmnew,&
                  spinup,delsoilM,delsoilT,&
                  output,&
                  patchout,&
                  check,&
                  verbose,leaps,logn,fixedCO2, &
                  spincasainput,     &
                  spincasa,          &
                  l_casacnp, l_laiFeedbk, l_vcmaxFeedbk, &
                  icycle,            &
                  casafile,          &
                  ncciy,             &
                  gswpfile,          &
                  redistrb, wiltParam, satuParam, &
                  cable_user
   !___END header


   !___ initialize etc

   ! there are still switches in core code which distinguish 
   ! between mode in which CABLE is being run
   cable_runtime%offline = .TRUE.
   
   ! Open, read and close the namelist file.
   OPEN(10,FILE='cable.nml')
      READ(10,NML=CABLE)
   CLOSE(10)
 
   IF (l_casacnp  .AND. (icycle == 0 .OR. icycle > 3)) &
                STOP 'icycle must be 1 to 3 when using CASA-CNP'
   IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) &
                STOP 'CASA-CNP required to get prognostic LAI or Vcmax'
   IF (l_vcmaxFeedbk .AND. icycle < 2) &
                STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
   IF (icycle > 0 .AND. (.NOT. soilparmnew)) &
                STOP 'CASA-CNP must use new soil parameters'
 
   ! Open log file:
   OPEN(logn,FILE=filename%log)
 
   ! IF this is a gswp run, check we have input
   IF (ncciy /= 0) THEN
     PRINT *, 'Looking for global offline run info.'
     IF (ncciy < 1986 .OR. ncciy > 1995) THEN
       PRINT *, 'Year ', ncciy, ' outside range of dataset!'
       PRINT *, 'Please check input in namelist file.'
       STOP
     ELSE
       CALL prepareFiles(ncciy)
     ENDIF
   ENDIF
 
   ! Open met data and get site information from netcdf file.
   ! This retrieves time step size, number of timesteps, starting date,
   ! latitudes, longitudes, number of sites. 
   CALL open_met_file(dels,kend,spinup)
 
   ! Checks where parameters and initialisations should be loaded from.
   ! If they can be found in either the met file or restart file, they will 
   ! load from there, with the met file taking precedence. Otherwise, they'll
   ! be chosen from a coarse global grid of veg and soil types, based on 
   ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
   CALL load_parameters(met,air,ssoil,veg,bgc,soil,canopy, &
        rough,rad,sum_flux,bal,logn,vegparmnew, &
        casabiome,casapool,casaflux,casamet,casabal,phen)
 
   ! Open output file:
   CALL open_output_file(dels,soil,veg,bgc,rough)
   
   tstep = 0          ! initialise spinup time step
   spinConv = .FALSE. ! initialise spinup convergence variable
 
   !___ END initialize etc


 
   !___ run model: including spinup if .TRUE.i

   ! Note: outermost DO loop facilitates EXIT in ELSE condition
   DO
      
      ktau_gl=0
      
      ! time step loop:
      DO ktau=kstart, kend ! time step loop
         
         ! increment total timstep counter
         tstep = tstep + 1
         ktau_gl= ktau_gl + 1

         canopy%oldcansto=canopy%cansto
 
         ! Get met data and LAI, set time variables.
         ! Rainfall input may be augmented for spinup purposes:
         CALL get_met_data(spinup,spinConv,met,soil,rad,veg,kend,dels) 
 
         ! Feedback prognostic vcmax and daily LAI from CASA-CNP to CABLE
         IF (l_vcmaxFeedbk) CALL casa_feedback( ktau_gl, veg, casabiome, &
                                       casapool, casamet )

         IF (l_laiFeedbk)   veg%vlai(:) = casamet%glai(:)
 
         ! CALL land surface scheme for this timestep, all grid points:
         CALL cbm(dels, air, bgc, canopy, met, bal, rad, rough, soil, &
                     ssoil, sum_flux, veg)
 
         if(icycle >0) then
            call bgcdriver(ktau_gl,kstart,kend,dels,met,ssoil,canopy,veg,soil, &
                   casabiome,casapool,casaflux,casamet,casabal,phen)
         endif 
 
         ! sumcflux is pulled out of subroutine cbm
         ! so that CASA-CNP can be called before adding the fluxes (Feb 2008, YP)
         CALL sumcflux(ktau_gl, kstart, kend, dels, bgc, canopy,  &
                        soil, ssoil, sum_flux, veg, met, casaflux, &
                        l_vcmaxFeedbk )
 
         ! Write time step's output to file if either: we're not spinning up 
         ! or we're spinning up and the spinup has converged:
         IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
              & (dels,met,canopy,ssoil,rad,bal,air,soil,veg)
 
      END DO ! END timestep loop
     
      
      ! see if spinup (if conducting one) has converged:
      IF(spinup.AND..NOT.spinConv) THEN
         
         ! Write to screen and log file:
         WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
                                    ' of data set complete...'
         
         WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
                                       ' of data set complete...'
       
         ! IF not 1st run through whole dataset:
         IF(INT(tstep/kend)>1) THEN 
            
            ! evaluate spinup
            IF(ANY(ABS(ssoil%wb-soilMtemp)>delsoilM).OR. &
                           ANY(ABS(ssoil%tgg-soilTtemp)>delsoilT)) THEN
           
               ! No complete convergence yet
               PRINT *, 'ssoil%wb : ', ssoil%wb
               PRINT *, 'soilMtemp: ', soilMtemp
               PRINT *, 'ssoil%tgg: ', ssoil%tgg
               PRINT *, 'soilTtemp: ', soilTtemp
         
            ELSE ! spinup has converged
         
               spinConv = .TRUE.
               
               ! Write to screen and log file:
               WRITE(*,'(A33)') ' Spinup has converged - final run'
               WRITE(logn,'(A52)') &
                ' Spinup has converged - final run - writing all data'
               WRITE(logn,'(A37,F7.5,A28)') &
                        ' Criteria: Change in soil moisture < ', &
                        delsoilM, ' in any layer over whole run'
               WRITE(logn,'(A40,F7.5,A28)' ) & 
                         '           Change in soil temperature < ', &
                           delsoilT, ' in any layer over whole run'
         
            END IF !end evaluate spin up
            
         ELSE ! allocate variables for storage
         
            ALLOCATE(soilMtemp(mp,ms), &
                & soilTtemp(mp,ms))
       
         END IF !END NOT 1st run through dataset

       ! store soil moisture and temperature
       soilTtemp = ssoil%tgg
       soilMtemp = REAL(ssoil%wb,r_1)

      ELSE
         
         ! if not spinning up, or spin up has converged, exit:
         EXIT

      END IF ! END if spinning up test convergence 
   
   END DO ! END run model loop outermost DO

   !___ END run model


   !___conditionally call CASA-CNP 
   IF (icycle > 0) THEN
     CALL casa_poolout(ktau_gl,veg,soil,casabiome,casapool,casaflux,casamet, &
                         casabal,phen)
     ktauday=int(24.0*3600.0/dels)
     nyear =INT((kend-kstart+1)/(365*ktauday))
     CALL casa_fluxout(nyear,veg,soil,casabal,casamet)
   END IF
   !___END call CASA-CNP 


   !___CABLE IO 
   ! Write restart file if requested:
   IF(output%restart) CALL create_restart(logn,dels,&
        soil,veg,ssoil,canopy,rough,rad,bgc,bal)
 
   ! Close met data input file:
   CALL close_met_file
   ! Close output file and deallocate main variables:
   CALL close_output_file(bal, air, &
        bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)
 
   ! Close log file
   CLOSE(logn)
   !___END CABLE IO 
 
   STOP

END PROGRAM offline_driver






!___SUBROUTINES used IF gswp run 

SUBROUTINE prepareFiles(ncciy)

!==============================================================================
! Name: prepareFiles
! Purpose: Temporary measure in development to ensure proper input files read 
! CALLed from: cable_driver
! MODULEs used:  define_dimensions 
! 
! CALLs: renameFiles
!
! Major contribution: land surface modeling team, CSIRO, Aspendale
!==============================================================================

   USE define_dimensions, ONLY: i_d
   USE io_variables, ONLY: logn,gswpfile
   IMPLICIT NONE
   INTEGER(i_d), INTENT(IN) :: ncciy

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

!==============================================================================
! Name: renameFiles
! Purpose: Temporary measure in development to ensure proper input files read 
! CALLed from: prepareFiles
! MODULEs used:   define_dimensions
!                 io_variables
!
! CALLs:
!
! Major contribution: land surface modeling team, CSIRO, Aspendale
!==============================================================================


   USE define_dimensions, ONLY: i_d
   IMPLICIT NONE
   INTEGER(i_d), INTENT(IN) :: logn
   INTEGER(i_d), INTENT(IN) :: nn
   INTEGER(i_d), INTENT(IN) :: ncciy
   CHARACTER(LEN=99), INTENT(INOUT) :: inFile
   CHARACTER(LEN=*),  INTENT(IN)    :: inName
   INTEGER(i_d) :: idummy
 
   READ(inFile(nn:nn+3),'(i4)') idummy
   IF (idummy < 1983 .OR. idummy > 1995) THEN
      PRINT *, 'Check position of the year number in input gswp file', inFile
      STOP
   ELSE
      WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
      WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
   ENDIF
    
END SUBROUTINE renameFiles

!___END SUBROUTINES used IF gswp run 
