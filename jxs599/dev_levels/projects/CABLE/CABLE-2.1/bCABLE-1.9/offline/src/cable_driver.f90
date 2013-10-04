! cable_driver.f90
!
! Netcdf offline driver for CABLE land surface scheme, May 2009.
! Gab Abramowitz, University of New South Wales, gabsun@gmail.com
!
! Thanks to Peter Isaac (Monash) for introducing the namelist file (Oct 2007)
!
! Modified and added casaCNP, Bernard Pak, Sep 2010.
!
! Please send bug reports to Bernard.Pak@csiro.au

PROGRAM offline_driver
  use cable_data_module
  USE cbm_module
  USE define_dimensions, ONLY:r_1,i_d,ms,mp,mvtype,mstype
  USE define_types
  USE io_variables, ONLY: logn,filename,gswpfile,ncciy,leaps, &
       verbose, fixedCO2,output,check,patchout,patch_type,soilparmnew, &
       redistrb, wiltParam, satuParam
  USE input_module, ONLY: open_met_file,load_parameters, &
       get_met_data,close_met_file
  USE output_module, ONLY: create_restart,open_output_file, &
       write_output,close_output_file
  ! new modules related to casaCNP
  USE casa_cnp   ! icycle declared in casadimension which is used in casa_cnp
                 ! casafile declared in casavariable, also used in casa_cnp
  use cable_common_module, only : ktau_gl, kend_gl, cable_user, cable_runtime
   
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

  !jhan: new data structure
  !type (model_cable), save :: cable
  !type (cable_auto), save :: cable_auto
  !type (model_canopy), save :: canopy_data
  
  REAL(r_1)         :: dels ! time step size in seconds
  INTEGER(i_d)      :: kstart ! start of simulation #
  LOGICAL    :: vegparmnew   ! using new format input file (BP dec 2007)
  LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
  LOGICAL    :: spinConv ! has spinup converged?
  LOGICAL    :: spincasainput
                   ! TRUE: input required to spin casacnp wil be saved;
                   ! FALSE: input will be read in to spin casacnp 1000 years
  LOGICAL    :: spincasa ! TRUE: casacnp will spin mloop times,
                         ! FALSE: no spin up
  LOGICAL    :: l_casacnp     ! using casaCNP with CABLE
  LOGICAL    :: l_laiFeedbk   ! using prognostic LAI
  LOGICAL    :: l_vcmaxFeedbk ! using prognostic Vcmax
  REAL(r_1)  :: delsoilM ! allowed variation in soil moisture for spin up
  REAL(r_1)  :: delsoilT ! allowed variation in soil temperature for spin up
  REAL(r_1),POINTER  :: soilMtemp(:,:) ! temporary storage for spin up
  REAL(r_1),POINTER  :: soilTtemp(:,:) ! temporary storage for spin up
  INTEGER(i_d) :: tstep  ! time step counter for spinup
  INTEGER(i_d) :: mloop  ! # spinup loops for casaCNP
  INTEGER(i_d) :: ktauday  !  day counter for casaCNP
  INTEGER(i_d) :: nyear    ! year counter for casaCNP
!  INTEGER(i_d) :: newktau   ! = ktau + krestart
!  INTEGER(i_d) :: krestart  ! index for restart timestep (for gswp run)
   integer(i_d) :: ktau !loop increment (corresponds to timestep)
  NAMELIST/CABLENML/filename, &
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

   cable_runtime%offline = .true. 

  !===================================================================!
  ! Open, read and close the namelist file.
  OPEN(10,FILE='cable.nml')
  READ(10,NML=CABLENML)
  CLOSE(10)
  !=====================================================================!

  IF (l_casacnp  .AND. (icycle == 0 .OR. icycle > 3)) &
               STOP 'icycle must be 1 to 3 when using casaCNP'
  IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) &
               STOP 'casaCNP required to get prognostic LAI or Vcmax'
  IF (l_vcmaxFeedbk .AND. icycle < 2) &
               STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
  IF (icycle > 0 .AND. (.NOT. soilparmnew)) &
               STOP 'casaCNP must use new soil parameters'

  ! Open log file:
  OPEN(logn,FILE=filename%log)

  ! Check for gswp run
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
  

  !jhan: alloc and point to old struture
  call allocate_cable_types( met, rad, soil, ssoil, canopy, rough, bal, veg)
  
  
  kstart = 1

!  if(icycle>0) then
!    if (spincasa) then
!      mloop=500
!      call spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
!                       casaflux,casamet,casabal,phen)
!    endif
!  endif

  tstep = 0          ! initialise spinup time step
  spinConv = .FALSE. ! initialise spinup convergence variable
  ! spinup loop:

  
DO
   ktau_gl=0
   kend_gl=kend
   ! time step loop:

   DO ktau=kstart, kend ! time step loop
      ! increment total timstep counter
      tstep = tstep + 1
      ktau_gl= ktau_gl + 1
      
      canopy%oldcansto=canopy%cansto

      !jhan: alloc and point to old struture
      call allocate_cable_auto(met, rad, soil, ssoil, canopy, rough, bal, veg)

      ! Get met data and LAI, set time variables.
      ! Rainfall input may be augmented for spinup purposes:
      CALL get_met_data(spinup,spinConv,met,soil,rad,veg,kend,dels) 

      ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
      IF (l_vcmaxFeedbk) CALL casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
      IF (l_laiFeedbk)   veg%vlai(:) = casamet%glai(:)
!      IF (icycle > 0) THEN
!        ! vcmax feedback
!        IF (icycle > 1) CALL casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
!        ! LAI feedback
!        veg%vlai(:) = casamet%glai(:)
!      ENDIF

      ! CALL land surface scheme for this timestep, all grid points:
      !CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
      CALL cbm(dels, air, bgc, canopy, met, &
             & bal, rad, rough, soil, ssoil, sum_flux, veg)
!      CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
!             & bal, rad, rough, soil, ssoil, sum_flux, veg, mvtype, mstype)

      if(icycle >0) then
        call bgcdriver(ktau_gl,kstart,kend,dels,met,ssoil,canopy,veg,soil, &
                  casabiome,casapool,casaflux,casamet,casabal,phen)

      ! sumcflux is pulled out of subroutine cbm
      ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
      CALL sumcflux(ktau_gl, kstart, kend, dels, bgc, canopy,  &
                  & soil, ssoil, sum_flux, veg, met, casaflux)
      endif 
      ! Write time step's output to file if either: we're not spinning up 
      ! or we're spinning up and the spinup has converged:
      IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
             & (dels,met,canopy,ssoil,rad,bal,air,soil,veg)


      !jhan: alloc and point to old struture
      call deallocate_cable_auto()

   
   END DO
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
        END IF
      ELSE ! allocate variables for storage
        ALLOCATE(soilMtemp(mp,ms), &
               & soilTtemp(mp,ms))
      END IF
      ! store soil moisture and temperature
      soilTtemp = ssoil%tgg
      soilMtemp = REAL(ssoil%wb,r_1)
    ELSE
      ! if not spinning up, or spin up has converged, exit:
      EXIT
    END IF
    
    
  END DO

  IF (icycle > 0) THEN
    CALL casa_poolout(ktau_gl,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
    ktauday=int(24.0*3600.0/dels)
    nyear =INT((kend-kstart+1)/(365*ktauday))
    CALL casa_fluxout(nyear,veg,soil,casabal,casamet)
  END IF

  ! Write restart file if requested:
  IF(output%restart) CALL create_restart(logn,dels,&
       soil,veg,ssoil,canopy,rough,rad,bgc,bal)
!  IF(output%restart) CALL create_restart(logn,ktau_gl,dels,&
!       soil,veg,ssoil,canopy,rough,rad,bgc,bal,mvtype,mstype)
! print *, sum_flux%sumpn, sum_flux%sumrp, sum_flux%sumrd, bal%wbal_tot, bal%ebal_tot

  ! Close met data input file:
  CALL close_met_file
  ! Close output file and deallocate main variables:
  CALL close_output_file(bal, air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)

   stop
END PROGRAM offline_driver


SUBROUTINE prepareFiles(ncciy)
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
