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
!
PROGRAM offline_driver
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
use cable_common_module, only : ktau_gl, cable_user, jhprint
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
  REAL(r_1)         :: dels ! time step size in seconds
  INTEGER(i_d)      :: kstart ! start of simulation #
  INTEGER(i_d)      :: ktau    ! index of time step = 1 ..  kend
  LOGICAL    :: vegparmnew   ! using new format input file (BP dec 2007)
  LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
  LOGICAL    :: spinConv ! has spinup converged?
  LOGICAL    :: spincasainput
                   ! TRUE: input required to spin casacnp wil be saved;
                   ! FALSE: input will be read in to spin casacnp 1000 years
  LOGICAL    :: spincasa ! TRUE: casacnp will spin mloop times,
                         ! FALSE: no spin up
  REAL(r_1)  :: delsoilM ! allowed variation in soil moisture for spin up
  REAL(r_1)  :: delsoilT ! allowed variation in soil temperature for spin up
  REAL(r_1),allocatable, save :: soilMtemp(:,:) ! temporary storage for spin up
  REAL(r_1),allocatable, save :: soilTtemp(:,:) ! temporary storage for spin up
  INTEGER(i_d) :: tstep  ! time step counter for spinup
  INTEGER(i_d) :: mloop  ! # spinup loops for casaCNP
   integer :: jhi =1
!  INTEGER(i_d) :: newktau   ! = ktau + krestart
!  INTEGER(i_d) :: krestart  ! index for restart timestep (for gswp run)
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
                 icycle,            &
                 casafile,          &
                 ncciy,             &
                 gswpfile,          &
                 redistrb, wiltParam, satuParam, &
                 cable_user
  !===================================================================!
  ! Open, read and close the namelist file.
  OPEN(10,FILE='cable.nml')
  READ(10,NML=CABLE)
  CLOSE(10)
  !=====================================================================!

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
       rough,rad,sum_flux,bal,logn,vegparmnew)

  ! Open output file:
  CALL open_output_file(dels,soil,veg,bgc,rough)

!   call new_init_temp(met,air,ssoil,veg,bgc,soil,canopy, &
!       rough,rad,sum_flux,bal,logn,vegparmnew)
!
!   subroutine new_init_temp(met,air,ssoil,veg,bgc,soil,canopy, &
!       rough,rad,sum_flux,bal,logn,vegparmnew)
!      implicit none
!      TYPE (met_type), INTENT(INOUT) :: met
!      TYPE (air_type), INTENT(INOUT) :: air
!      TYPE (soil_snow_type), INTENT(OUT) :: ssoil
!      TYPE (veg_parameter_type), INTENT(OUT)  :: veg
!      TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
!      TYPE (soil_parameter_type), INTENT(OUT) :: soil
!      TYPE (canopy_type), INTENT(OUT)    :: canopy
!      TYPE (roughness_type), INTENT(OUT) :: rough
!      TYPE (radiation_type),INTENT(OUT)  :: rad
!      TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
!      TYPE (balances_type), INTENT(OUT)  :: bal
!      INTEGER,INTENT(IN) :: logn     ! log file unit number
!      LOGICAL,INTENT(IN) :: vegparmnew  ! are we using the new format?
!
!      rough%zref_uv = rough%z
!      canopy%zetar = 
!
!      return
!   end subroutine new_init_temp

!  print *, 'mp mstype mvtype = ',mp,mstype,mvtype
  if(icycle>0) then
    mloop=500
    kstart=1
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    call alloc_phenvariable(phen,mp)

    call casa_readpoint(veg,soil,casaflux,casamet,rad)
    call casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    call casa_readphen(veg,casamet,phen)
    call casa_init(casapool,casabal,veg)
!  print *, 'mp mstype mvtype = ',mp,mstype,mvtype
!    if (spincasa) then
!      call spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
!                       casaflux,casamet,casabal,phen)
!    endif
  endif

  kstart = 1
  tstep = 0          ! initialise spinup time step
  spinConv = .FALSE. ! initialise spinup convergence variable
  ! spinup loop:
call jhprint('kstart', kstart)
call jhprint('kend', kend)
DO
   ktau_gl=0
   ! time step loop:
   call jhprint('outer loop', jhi)
   DO ktau = kstart, kend ! time step loop
      ! increment total timstep counter
      tstep = tstep + 1
      ktau_gl= ktau_gl + 1
! print *, 'tstep = ', tstep
      canopy%oldcansto=canopy%cansto
      if (jhi>1)  call jhprint('ktau_gl', ktau_gl ) 
      ! Get met data and LAI, set time variables.
      ! Rainfall input may be augmented for spinup purposes:
      CALL get_met_data(spinup,spinConv,met,soil,rad,veg,kend,dels) 
        
      ! CALL land surface scheme for this timestep, all grid points:
      !CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
      CALL cbm(dels, air, bgc, canopy, met, &
             & bal, rad, rough, soil, ssoil, sum_flux, veg)
!      CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
!             & bal, rad, rough, soil, ssoil, sum_flux, veg, mvtype, mstype)

      if(icycle >0) then
        call bgcdriver(ktau,kstart,kend,dels,met,ssoil,canopy,veg,soil, &
                  casabiome,casapool,casaflux,casamet,casabal,phen)
      endif 

      ! sumcflux is pulled out of subroutine cbm
      ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
      CALL sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
                  & soil, ssoil, sum_flux, veg, met, casaflux)

      ! Write time step's output to file if either: we're not spinning up 
      ! or we're spinning up and the spinup has converged:
      IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
             & (dels,met,canopy,ssoil,rad,bal,air,soil,veg)

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
         
         if (.NOT. allocated(soilMtemp) ) &
            ALLOCATE(soilMtemp(mp,ms), &
                     soilTtemp(mp,ms))
      END IF
   
      ! store soil moisture and temperature

      call jhprint ('shape soilTtemp ', shape(soilTtemp) )  
      call jhprint ('shape soilMtemp ', shape(soilMtemp) )  
      call jhprint ('shape ssoil%tgg ', shape(ssoil%tgg) )  
      call jhprint ('shape ssoil%wb ', shape(ssoil%wb) )  
      
      soilTtemp = ssoil%tgg
      soilMtemp = REAL(ssoil%wb,r_1)
      call jhprint ('ssoil%wb done ', real(ssoil%wb) )  

   ELSE
      print *, "are about to exit spinup cond"
      ! if not spinning up, or spin up has converged, exit:
      EXIT
   END IF
   jhi = jhi+1
END DO

  ! Write restart file if requested:
  IF(output%restart) CALL create_restart(logn,dels,&
       soil,veg,ssoil,canopy,rough,rad,bgc,bal)
!  IF(output%restart) CALL create_restart(logn,ktau,dels,&
!       soil,veg,ssoil,canopy,rough,rad,bgc,bal,mvtype,mstype)
! print *, sum_flux%sumpn, sum_flux%sumrp, sum_flux%sumrd, bal%wbal_tot, bal%ebal_tot

  ! Write final pool sizes for casaCNP
  IF (icycle>0) call casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux, &
                                  casamet,casabal,phen)
  
  ! Close met data input file:
  CALL close_met_file
  ! Close output file and deallocate main variables:
  CALL close_output_file(bal, air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)

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
