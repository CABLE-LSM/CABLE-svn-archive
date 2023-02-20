! this file includes codes for initilization, input and output related to mic_cnp.F90
! Strategy for implementing MIMICS into CABLE
!
!(1) allocate all the variables at the start
!(2) read the parameter lookup table
!    make the following parameters with (mp) dimensions and add to the "allocate"
!    xav,xak,xdesorp,xbeta,xdiffsoc and rootbetax
!(3) assign all the parameter values to 1:mp
!(4) each timestep: get the litter input from the previous day, and current time 
!    of soil temperature, moisture
!(5) output soil respiration and update pool sizes
!(6) add "mic" as a switch (=1,2 and 3 for C, CN and CNP cycles)
!(7) call "vmic" from "biogeochem" in "casa_inout.F90"
! strategy based on casa-cnp in CABLE
! (8) check the unit of "deltvmic" for consistency
!
!===========================================================================================
MODULE vmic_inout_mod
  USE cable_def_types_mod    !, ONLY : mp,ms,r_2,mvtype,mstype
  USE casavariable,             ONLY : casa_biome
  USE vmic_constant_mod
  USE vmic_variable_mod
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE vmic_allocate(micparam,micinput,micoutput,miccpool,micnpool)
    implicit none
    TYPE(mic_parameter)  :: micparam
    TYPE(mic_input)      :: micinput
    TYPE(mic_cpool)      :: miccpool
    TYPE(mic_npool)      :: micnpool
    TYPE(mic_output)     :: micoutput

      call vmic_allocate_parameter(micparam)
      call vmic_allocate_input(micinput)
      call vmic_allocate_output(micoutput)
      call vmic_allocate_cpool(miccpool)
      call vmic_allocate_npool(micnpool)

  end SUBROUTINE vmic_allocate 


  SUBROUTINE vmic_parameter(veg,soil,casabiome,micparam,micinput,micfile)
  ! read the parameter lookup table and assigm them to (1:mp)
    implicit none
    TYPE(mic_parameter)                    :: micparam
    TYPE(mic_input)                        :: micinput
    TYPE(micfile_type)                     :: micfile
    TYPE (veg_parameter_type),  INTENT(IN) :: veg
    TYPE (soil_parameter_type), INTENT(IN) :: soil
    TYPE (casa_biome),          INTENT(IN) :: casabiome
    real(r_2), dimension(17)               :: xav_pft,      &
                                              xak_pft,      &
                                              xavexp_pft,   &
                                              xbeta_pft,    &
                                              xdiffsoc_pft, &
                                              xdesorpt_pft, &
                                              rootbeta_pft, & 
                                              xcnleaf,      &
                                              xcnwood,      &
                                              xcnroot,      &
                                              fligleaf,     &
                                              fligwood,     &
                                              fligroot
								
    ! parameter estimates based on ICP-China sites for 4 PFTs based optimizinfg 5 parameters
    ! xav      xak       xavexp    xbeta    xdiffsoc   
    ! 14.524    15.091     1.008     1.255    18.329
    !  5.475    13.271     1.011     1.446     4.636
    !  6.125     9.869     1.324     1.582     1.051
    !  4.757    19.013     1.271     1.955     0.263
    !
    !  7.72     14.31      1.15      1.56      6.07

!    data      xav_pft/ 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72/
!    data      xak_pft/14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31/
!    data   xavexp_pft/ 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15/
!    data    xbeta_pft/ 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56/
!    data xdiffsoc_pft/ 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07/
    data rootbeta_pft/ 0.96, 0.96, 0.96, 0.98, 0.96, 0.97, 0.94, 0.96, 0.97, 0.94, 0.97, 0.97, 0.97, 0.97, 0.96, 0.96, 0.96/   

    ! local variables
    integer ivt,np,ns,isoil
   
      xdesorpt_pft(:) = 1.0
      open(201,file=micfile%micbiome)
      print *, 'reading micparameter file', micfile%micbiome  
      read(201,*)
      read(201,*) (xav_pft(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xak_pft(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xavexp_pft(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xdiffsoc_pft(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xcnleaf(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xcnwood(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (xcnroot(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (fligleaf(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (fligwood(ivt),ivt=1,mvtype)
      read(201,*)
      read(201,*) (fligroot(ivt),ivt=1,mvtype)      
      close(201)
     
        do np=1,mp
           ivt          = veg%iveg(np)
	       isoil        = soil%isoilm(np)
           micparam%xav(np)      = xav_pft(ivt)
           micparam%xak(np)      = xak_pft(ivt)
           micparam%xavexp(np)   = xavexp_pft(ivt)	  
           micparam%xbeta(np)    = xbeta_pft(ivt)
           micparam%xdiffsoc(np) = xdiffsoc_pft(ivt)
           micparam%rootbetax(np)= rootbeta_pft(ivt)
           micparam%xdesorp(np)  = xdesorpt_pft(ivt)
           micparam%xcnleaf(np)  = xcnleaf(ivt)
           micparam%xcnroot(np)  = xcnroot(ivt)	  
           micparam%xcnwood(np)  = xcnwood(ivt)
           micparam%fligleaf(np) = fligleaf(ivt)
           micparam%fligroot(np) = fligroot(ivt)
           micparam%fligwood(np) = fligwood(ivt)
           ! soil parameters	  
!           micparam%bulkd(np,1:ms) = soil%rhosoil(np)
!           micparam%clay(np,1:ms)  = soil%clay(np)	
           ! vegetation parameters
       
!           micparam%pft(np)     =  veg%iveg(np)	   
!           micparam%xcnleaf(np) = 1.0/casabiome%ratioNCplantmax(ivt,1)
!           micparam%xcnroot(np) = 1.0/casabiome%ratioNCplantmax(ivt,3)	  
!           micparam%xcnwood(np) = 1.0/casabiome%ratioNCplantmax(ivt,2)
!           micparam%fligleaf(np)=  casabiome%fracligninplant(ivt,1)
!           micparam%fligroot(np)=  casabiome%fracligninplant(ivt,3)
!           micparam%fligwood(np)=  casabiome%fracligninplant(ivt,2)
        enddo  
!        do ns=1,ms
!           micparam%zse(ns) = soil%zse(ns)
!        enddo
	
    end SUBROUTINE vmic_parameter

  SUBROUTINE vmic_param_constant(veg,soil,micparam)
    USE vmic_constant_mod
    USE vmic_variable_mod
	USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type
    implicit none
    TYPE (veg_parameter_type),  INTENT(IN)    :: veg	
    TYPE (soil_parameter_type), INTENT(IN)    :: soil	
    TYPE (mic_parameter),       INTENT(INout) :: micparam
    !local variables   
    integer np,ns
    real(r_2) Kor, Kok, Q1, Q2, fm, fs

    !xdiffsoc=xopt(5); rootbetax=xopt(7)
  
      Kor = 4.0; Kok = 4.0; Q1= Kor; Q2  = Kok; fm = 0.05; fs= 0.05
      micparam%Q1(:,:) = Q1; micparam%Q2(:,:) =Q2; micparam%fm(:,:)=fm; micparam%fs(:,:)=fs

      ! calculate mp by ms all parameter values
      do np=1, mp
!         do ns=1,ms
!            micparam%sdepth(np,ns)   = soil%zse(ns)
!            micparam%fracroot(np,ns) = veg%froot(np,ns)
!          enddo !"ns"
          micparam%diffsocx(np) = micparam%xdiffsoc(np) * diffsoc  !"diffsoc" from mic_constant
      enddo    ! "np=1,mp"
  
      if(diag==1) then
!         print *, micparam%fracroot(outp,:) 
!         print *, micparam%sdepth(outp,:)
         print *, micparam%diffsocx(outp)
      endif
   END SUBROUTINE vmic_param_constant  

   SUBROUTINE vmic_init(micparam,micinput,miccpool,micnpool)
    USE vmic_constant_mod
    USE vmic_variable_mod
    implicit none

    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_cpool),     INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    ! local variables
    ! for numerical solution
    real(r_2),    parameter            :: tol = 1.0E-04
    real(r_2),    parameter            :: tolx = 0.0001
    real(r_2),    parameter            :: tolf = 0.000001
    integer,     parameter             :: ntrial = 100
    integer np,ns,ip
    real(r_2),    dimension(mcpool)    :: cpooldef,xpool0,y

!	  print *, 'calling vmic_init'

      cpooldef(1) = 16.5*0.1;     cpooldef(2) = 16.5*0.1
      cpooldef(3) = 16.5*0.025;   cpooldef(4) = 16.5*0.025
      cpooldef(5) = 16.5*0.1125;  cpooldef(6) = 16.5*0.375;  cpooldef(7) = 16.5*0.2625

      do ip=1,mcpool
         miccpool%cpool(:,:,ip) = cpooldef(ip)
      enddo
  
!      call vmic_param_time(micparam,micinput,micnpool)
!  
!      do np=1,mp
!            do ns=1,ms
!
!               ! initial pool sizes
!               do ip=1,mcpool
!                  xpool0(ip) = miccpool%cpool(np,ns,ip)
!               enddo
!
!!!!               call mnewt(ntrial,np,ns,kinetics,micparam,micinput,xpool0,tolx,tolf)
!               call vmic_c(np,ns,micparam,micinput,xpool0,y)
!               if(maxval(xpool0(1:mcpool))>1.0e4.or.minval(xpool0(1:mcpool))<0.0) then
!                  xpool0 = cpooldef
!               endif
!               do ip=1,mcpool
!                  miccpool%cpool(:,:,ip) = xpool0(ip)
!               enddo
!            enddo
!       enddo
  

  END SUBROUTINE vmic_init

   SUBROUTINE vmic_input(dleaf,dwood,droot,nsoilmin,veg,casamet,micparam,micinput,micnpool)
     USE casavariable,                    ONLY : casa_met
     USE casaparm,                        ONLY : tkzeroc
     implicit none
     TYPE (veg_parameter_type),INTENT(IN)     :: veg  ! vegetation parameters
     TYPE (casa_met),          INTENT(IN)     :: casamet
     TYPE(mic_parameter),      INTENT(IN)     :: micparam
     TYPE(mic_input),          INTENT(INOUT)  :: micinput
     TYPE(mic_npool),          INTENT(INOUT)  :: micnpool   
     real(r_2),  dimension(mp)                :: dleaf,dwood,droot,fcnpp,nsoilmin
     real(r_2),  dimension(17)                :: fcnpp_pft
     data fcnpp_pft/500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0/
 
     ! local variables
     integer np,ivt
  
       do np=1,mp
      !
      !    ivt =micparam%pft(np)
          ivt = veg%iveg(np)
          micinput%fcnpp(np) = fcnpp_pft(ivt)         ! gc/m2/year for calculating microbial turnover rate (constant!!)
          micinput%dleaf(np) = dleaf(np)              ! gc/m2/deltvmic
          micinput%dwood(np) = dwood(np)              ! gc/m2/deltvmic
          micinput%droot(np) = droot(np)              ! gc/m2/deltvmic
          micinput%tavg(np,:)= casamet%tsoil(np,:)- tkzeroc
          micinput%wavg(np,:)= casamet%moist(np,:)
          micnpool%mineralN(np,1:ms) = nsoilmin(np)   ! temporary solution
   enddo
  end SUBROUTINE vmic_input

  SUBROUTINE vmic_driver(ktau,dels,idoY,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                         casamet,casabal,phen,micparam,micinput,miccpool,micnpool,micoutput)
    USE casadimension
    USE casa_cnp_module
    USE casavariable
    use vmic_carbon_cycle_mod
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: ktau
    REAL,    INTENT(IN)    :: dels
    INTEGER, INTENT(IN)    :: idoy
    INTEGER, INTENT(IN)    :: LALLOC
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (phen_variable),         INTENT(INOUT) :: phen
    TYPE (mic_parameter),         INTENT(INOUT) :: micparam
    TYPE (mic_input),             INTENT(INOUT) :: micinput
    TYPE (mic_cpool),             INTENT(INOUT) :: miccpool
    TYPE (mic_npool),             INTENT(INOUT) :: micnpool
    TYPE (mic_output),            INTENT(INOUT) :: micoutput

    ! local variables added by ypwang following Chris Lu 5/nov/2012

    REAL, DIMENSION(mp)  :: cleaf2met,cleaf2str,                                  &
                            croot2met,croot2str,cwood2cwd,                        &
                            nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,    &
                            pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd
    REAL(r_2), dimension(mp)           :: dleaf,dwood,droot
    ! local variables
    real(r_2),    dimension(mcpool)    :: xpool0,xpool1
    real(r_2),    dimension(ms)        :: ypooli,ypoole,fluxsoc
    REAL(r_2),    DIMENSION(mp)        :: xnplimit,xNPuptake
    REAL(r_2),    DIMENSION(mp)        :: xklitter,xksoil,xkNlimiting
    REAL(r_2),    DIMENSION(mp)        :: xkleafcold,xkleafdry,xkleaf
    REAL(r_2),    DIMENSION(mp)        :: nsoilmin
    REAL(r_2)                          :: delty,timex,diffsocxx
    real(r_2),    dimension(ms)        :: zsex
    integer ndeltvmic,ntime
    INTEGER  np,ip,ns,is
    integer ivegx,isox


    CALL phenology(idoy,veg,phen)
    CALL avgsoil(veg,soil,casamet)
    CALL casa_rplant1(veg,casabiome,casapool,casaflux,casamet)

    IF (.NOT.cable_user%CALL_POP) THEN
       CALL casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
    ENDIF

    CALL casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
         casamet,phen)
    CALL casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
         casaflux,casamet,phen)

    ! changed by ypwang following Chris Lu on 5/nov/2012
    CALL casa_delplant(veg,casabiome,casapool,casaflux,casamet,            &
                       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
 
    call vmic_param_time(veg,soil,micparam,micinput,micnpool)
! litter fall flux in g C/m2/delt, where "delt" is one day in casa-cnp, and timesetp in vmic is hourly
    dleaf(:) = (cleaf2met(:)+cleaf2str(:))/24.0
    dwood(:) =  cwood2cwd(:)/24.0
    droot(:) = (croot2met(:)+croot2str(:))/24.0
    nsoilmin(:) = 2.0    ! mg N/kg soil

    call vmic_input(dleaf,dwood,droot,nsoilmin,veg,casamet,micparam,micinput,micnpool)
    
    do np=1,mp
       ! do MIMICS for each soil layer
        do ns=1,ms
           do ip=1,mcpool
              xpool0(ip) = miccpool%cpool(np,ns,ip)
           enddo

           ndeltvmic=24  ! 24-hourly
           delty = real(ndeltvmic) * deltvmic  ! time step in rk4 in "delt"
           timex = real(idoy)
   
           do ntime=1,1
             call rk4modelx(timex,delty,np,ns,micparam,micinput,xpool0,xpool1)
             ! the following used to avoid poolsize<0.0
             do ip=1,mcpool
                xpool0(ip) = max(1.0e-8,xpool1(ip))
             enddo
           enddo

           xpool1=xpool0 
           do ip=1,mcpool
              miccpool%cpool(np,ns,ip) = xpool1(ip)
           enddo

           ! soil respiration in mg c/cm3/deltvmic		   
           micoutput%rsoil(np,ns) = micinput%dleaf(np)+micinput%dwood(np)+micinput%droot(np) &
                                  - sum(xpool1(1:mcpool)-sum(xpool0(1:mcpool)))
        enddo ! "ns"
  
        if(diag==1) then  
           print *, 'np1', outp,micparam%diffsocx(outp)
           do ns=1,ms
              print *, ns, miccpool%cpool(outp,ns,:) 
           enddo  
        endif
  
        do ip=1,mcpool
           do ns=1,ms
              ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
           enddo  !"ns"

           fluxsoc(:) = 0.0  ! This flux is added in "modelx"
           diffsocxx= micparam%diffsocx(np)
  
           ! only do every 24*deltvmic  
!           call bioturb(ndeltvmic,ms,micparam%zse(1:ms),deltvmic,diffsocxx,fluxsoc,ypooli,ypoole) 
           zsex = soil%zse 
           call bioturb(ndeltvmic,ms,zsex,deltvmic,diffsocxx,fluxsoc,ypooli,ypoole)  
           do ns=1,ms
              miccpool%cpool(np,ns,ip) = ypoole(ns)
           enddo
        enddo !"ip"

        enddo !"np" 
  
  END SUBROUTINE vmic_driver
  
  SUBROUTINE write_mic_output_nc ( miccpool, micnpool, micoutput, micfile, &
       ctime, FINAL )

    USE CASAVARIABLE
    USE CABLE_COMMON_MODULE
    USE casa_ncdf_module, ONLY: HANDLE_ERR
    USE vmic_constant_mod, ONLY: mcpool
    USE vmic_variable_mod, ONLY: mic_cpool, mic_npool, mic_output, micfile_type


    USE cable_def_types_mod, ONLY: veg_parameter_type

    USE netcdf

    IMPLICIT NONE

    TYPE(mic_cpool), INTENT(IN)      :: miccpool
    TYPE(mic_npool), INTENT(IN)      :: micnpool
    TYPE(mic_output), INTENT(IN)     :: micoutput
    TYPE(micfile_type), INTENT(IN)   :: micfile

    INTEGER   :: STATUS, ctime
    INTEGER   :: mp_ID, soil_ID, miccarb_ID, t_ID, i
    LOGICAL   :: FINAL
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    INTEGER, SAVE :: time_ID, rsoil_ID, nmic_ID, cmic_ID
    INTEGER, SAVE :: FILE_ID, CNT = 0
    LOGICAL   :: EXRST
    CHARACTER(len=50) :: RecordDimName


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       IF (LEN( TRIM(micfile%micoutput) ) .GT. 0) THEN
          fname=TRIM(micfile%micoutput)
       ELSE
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_mic_out.nc'
       ENDIF
       INQUIRE( FILE=TRIM( fname ), EXIST=EXRST )
       EXRST = .FALSE.
       IF ( EXRST ) THEN
          STATUS = NF90_open(fname, mode=nf90_write, ncid=FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          CALL1 = .FALSE.

          STATUS = nf90_inq_dimid(FILE_ID, 'time', t_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = nf90_inq_varid(FILE_ID, 'time', time_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = nf90_inq_varid(FILE_ID,'mic_rsoil' , rsoil_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = nf90_inq_varid(FILE_ID,'nmic_ID' , nmic_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = nf90_inq_varid(FILE_ID,'cmic_ID', cmic_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ELSE
          ! Create NetCDF file:
          STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Put the file in define mode:
          STATUS = NF90_redef(FILE_ID)

          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle"   , icycle  )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

          ! Define dimensions:
          ! Land (number of points)
          STATUS = NF90_def_dim(FILE_ID, 'mp'   , mp   , mp_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'soil'  , ms  , soil_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'mic_carbon_pools' , mcpool , miccarb_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define variables
          STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),time_ID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = NF90_def_var(FILE_ID, 'mic_rsoil' ,NF90_FLOAT,(/mp_ID,soil_ID,t_ID/),rsoil_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = NF90_def_var(FILE_ID, 'mic_npool' ,NF90_FLOAT,(/mp_ID,soil_ID,t_ID/),nmic_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          STATUS = NF90_def_var(FILE_ID, 'mic_cpool' ,NF90_FLOAT, &
               (/mp_ID,soil_ID,miccarb_ID,t_ID/),cmic_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! End define mode:
          STATUS = NF90_enddef(FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          CALL1 = .FALSE.
       ENDIF !( EXRST )
    ENDIF

    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, time_ID, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, rsoil_ID, REAL(micoutput%rsoil,4),   &
            start=(/ 1,1,CNT /), count=(/ mp,ms,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, nmic_ID, REAL(micnpool%mineralN,4),   &
            start=(/ 1,1,CNT /), count=(/ mp,ms,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, cmic_ID, REAL(miccpool%cpool,4),   &
            start=(/ 1,1,1,CNT /), count=(/ mp,ms,mcpool,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " Casa Output written to ",fname
    ENDIF

  END SUBROUTINE write_mic_output_nc
END MODULE vmic_inout_mod
