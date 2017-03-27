SUBROUTINE get_casa_restart(casamet,casapool,casabal,phen)
!SUBROUTINE get_casa_restart(casabiome,casamet,casapool,casabal,veg,phen)
! read pool sizes from restart file, then initialize a few more related vars.
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_def_types_mod
  USE cable_io_vars_module, ONLY: landpt, patch, max_vegpatches
  USE cable_common_module,  ONLY: filename
  USE cable_read_module
  USE cable_abort_module,  ONLY: nc_abort
  USE netcdf

  IMPLICIT NONE
!  TYPE (casa_biome),   INTENT(IN)    :: casabiome
  TYPE (casa_met),     INTENT(INOUT) :: casamet
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
!  TYPE (veg_parameter_type), INTENT(IN) :: veg
  TYPE (phen_variable),   INTENT(INOUT) :: phen

  ! local variables
  INTEGER :: ncid, ncok, poolP, plantID, npt
  ! added by ypw to reset some points with wrong pool sizes (NaN values)
  INTEGER  nvt
  INTEGER, dimension(mp)             :: iveg_avg
  INTEGER, dimension(mvtype)         :: npt_avg
  REAL,    dimension(mvtype)         :: avglai,avgclab,avgnsoilmin,avgplab,avgpsorb,avgpocc
  REAL,    dimension(mvtype,mplant)  :: avgcplant,  avgnplant,  avgpplant
  REAL,    dimension(mvtype,mlitter) :: avgclitter, avgnlitter, avgplitter
  REAL,    dimension(mvtype,msoil)   :: avgcsoil,   avgnsoil,   avgpsoil
  REAL, parameter :: crootmin = 0.01       ! gc m-2
  REAL, parameter :: crootmax = 100000.0  ! gc m-2

  LOGICAL ::                                                                  &
       from_restart = .TRUE., & ! insist variables/params load
       dummy        = .TRUE.    ! To replace completeSet in readpar; unused

  ncok = NF90_OPEN(filename%restart_in,0,ncid) ! opened once before, less checks
  ncok = NF90_INQ_DIMID(ncid,'pools_plant',plantID)
  ncok = NF90_INQUIRE_DIMENSION(ncid,plantID,len=poolP)
  IF(ncok /= NF90_NOERR) CALL nc_abort                                       &
       (ncok,'Error finding number of plant pools in restart file '            &
       //TRIM(filename%restart_in)//' (SUBROUTINE get_casa_restart)')
  IF(poolP /= mplant) CALL abort('Number of plant pools in '//               &
        'restart file '//TRIM(filename%restart_in)//                         &
        ' differs from number in CASA_dimension')

  CALL readpar(ncid,'iveg',dummy,iveg_avg,filename%restart_in, &
                max_vegpatches,'def',from_restart,mp)
  CALL readpar(ncid,'LAI',dummy,casamet%glai,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'phase',dummy,phen%phase,filename%restart_in, &
                max_vegpatches,'def',from_restart,mp)
  CALL readpar(ncid,'Clabile',dummy,casapool%Clabile,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'CASA_Cplant',dummy,casapool%Cplant,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Clitter',dummy,casapool%Clitter,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'CASA_Csoil',dummy,casapool%Csoil,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Nplant',dummy,casapool%Nplant,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Nlitter',dummy,casapool%Nlitter,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Nsoil',dummy,casapool%Nsoil,filename%restart_in, & 
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Nsoilmin',dummy,casapool%Nsoilmin,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Pplant',dummy,casapool%Pplant,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Plitter',dummy,casapool%Plitter,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Psoil',dummy,casapool%Psoil,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Psoillab',dummy,casapool%Psoillab,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Psoilsorb',dummy,casapool%Psoilsorb,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)
  CALL readpar(ncid,'Psoilocc',dummy,casapool%Psoilocc,filename%restart_in, &
                max_vegpatches,'cnp',from_restart,mp)

  ncok = NF90_CLOSE(ncid)
  IF(ncok/=NF90_NOERR) CALL nc_abort(ncok,'Error closing restart file '     &
        //TRIM(filename%restart_in)// '(SUBROUTINE get_casa_restart)')

  if(initcasa==0) then
     do npt=1,mp 
        casamet%lon(npt) = patch(npt)%longitude 
        casamet%lat(npt) = patch(npt)%latitude 
     enddo   
  endif
  
 ! npt=26493
 ! print *, 'pool sizes for 26493 from restart: ', casapool%cplant(npt,:), casapool%nplant(npt,:)
  ! calculate PFT mean pool sizes
  !print *, 'calcualting PGT means and assign means to some NaN land points'
  npt_avg =0
  avgcplant =0.0; avgclitter=0.0; avgcsoil=0.0; avglai = 0.0; avgclab=0.0
  avgnplant =0.0; avgnlitter=0.0; avgnsoil=0.0; avgnsoilmin=0.0
  avgpplant =0.0; avgplitter=0.0; avgpsoil=0.0; avgplab=0.0; avgpsorb=0.0; avgpocc=0.0 
  avglai    =0.0
  do npt=1,mp
     nvt = iveg_avg(npt)
     if(nvt>0.and.nvt<=mvtype-7) then
     if(.not.(isnan(casapool%cplant(npt,1)).or.isnan(casapool%cplant(npt,3)) &
        .or.isnan(casapool%nplant(npt,1)).or.isnan(casapool%nplant(npt,3)) &
        .or.isnan(casapool%pplant(npt,1)).or.isnan(casapool%pplant(npt,3)))) then
     if(min(casapool%cplant(npt,1),casapool%cplant(npt,3)) >crootmin ) then
         npt_avg(nvt) = npt_avg(nvt) +1
         avgcplant(nvt,:)  = avgcplant(nvt,:)   + casapool%cplant(npt,:)
         avgnplant(nvt,:)  = avgnplant(nvt,:)   + casapool%nplant(npt,:)
         avgpplant(nvt,:)  = avgpplant(nvt,:)   + casapool%pplant(npt,:)
         avgclitter(nvt,:) = avgclitter(nvt,:)  + casapool%clitter(npt,:)
         avgnlitter(nvt,:) = avgnlitter(nvt,:)  + casapool%nlitter(npt,:)
         avgplitter(nvt,:) = avgplitter(nvt,:)  + casapool%plitter(npt,:)
         avgcsoil(nvt,:)   = avgcsoil(nvt,:)    + casapool%csoil(npt,:)
         avgnsoil(nvt,:)   = avgnsoil(nvt,:)    + casapool%nsoil(npt,:)
         avgpsoil(nvt,:)   = avgpsoil(nvt,:)    + casapool%psoil(npt,:)
         avgclab(nvt)      = avgclab(nvt)       + casapool%clabile(npt)
         avgnsoilmin(nvt)  = avgnsoilmin(nvt)   + casapool%nsoilmin(npt)
         avgplab(nvt)      = avgplab(nvt)       + casapool%psoillab(npt)
         avgpsorb(nvt)     = avgpsorb(nvt)      + casapool%psoilsorb(npt)
         avgpocc(nvt)      = avgpocc(nvt)       + casapool%psoilocc(npt)
         avglai(nvt)       = avglai(nvt)        + casamet%glai(npt)
     endif
     endif
     endif
  enddo

!  where (npt_avg>1) 
     avgcplant  = avgcplant/max(1.0,real(spread(npt_avg,2,mplant)))
     avgnplant  = avgnplant/max(1.0,real(spread(npt_avg,2,mplant)))
     avgpplant  = avgpplant/max(1.0,real(spread(npt_avg,2,mplant)))
     avgclitter = avgclitter/max(1.0,real(spread(npt_avg,2,mlitter)))
     avgnlitter = avgnlitter/max(1.0,real(spread(npt_avg,2,mlitter)))
     avgplitter = avgplitter/max(1.0,real(spread(npt_avg,2,mlitter)))
     avgcsoil   = avgcsoil/max(1.0,real(spread(npt_avg,2,msoil)))
     avgnsoil   = avgnsoil/max(1.0,real(spread(npt_avg,2,msoil)))
     avgpsoil   = avgpsoil/max(1.0,real(spread(npt_avg,2,msoil)))
     avgclab    = avgclab/max(1.0,real(npt_avg))
     avgnsoilmin  = avgnsoilmin/max(1.0,real(npt_avg))
     avgplab      = avgplab/max(1.0,real(npt_avg))
     avgpsorb     = avgpsorb/max(1.0,real(npt_avg))
     avgpocc      = avgpocc/max(1.0,real(npt_avg))
     avglai       = avglai/max(1.0,real(npt_avg))
       
!  endwhere    

  do npt=1,mp
     nvt = iveg_avg(npt)
     if(nvt>0.and.nvt<=mvtype-7) then
     if(isnan(casapool%cplant(npt,1)).or.isnan(casapool%cplant(npt,3)).or. &
        isnan(casapool%nplant(npt,1)).or.isnan(casapool%nplant(npt,3)).or. &
        isnan(casapool%pplant(npt,1)).or.isnan(casapool%pplant(npt,3)).or. &
        min(casapool%cplant(npt,1),casapool%cplant(npt,3)) <=crootmin) then

    !    print *, 'PFT mean pool sizes used for ', npt, 'vegtype ',nvt, &
    !   ' before: ', casapool%cplant(npt,:), 'after: ', avgcplant(nvt,:)
        
        casapool%cplant(npt,:)  = avgcplant(nvt,:)
        casapool%nplant(npt,:)  = avgnplant(nvt,:)
        casapool%pplant(npt,:)  = avgpplant(nvt,:)

        casapool%clitter(npt,:) = avgclitter(nvt,:)
        casapool%nlitter(npt,:) = avgnlitter(nvt,:)
        casapool%plitter(npt,:) = avgplitter(nvt,:)

        casapool%csoil(npt,:)   = avgcsoil(nvt,:)
        casapool%nsoil(npt,:)   = avgnsoil(nvt,:)
        casapool%psoil(npt,:)   = avgpsoil(nvt,:)

        casapool%clabile(npt)   = avgclab(nvt)
        casapool%nsoilmin(npt)  = avgnsoilmin(nvt)
        casapool%psoillab(npt)  = avgplab(nvt)
        casapool%psoilsorb(npt) = avgpsorb(nvt)
        casapool%psoilocc(npt)  = avgpocc(nvt)

        casamet%glai(npt)      =  avglai(nvt)

     endif
     endif
  enddo


  ! check pool sizes
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0
  casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
  casabal%FCnppyear=0.0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%Nsoil(:,:) / casapool%ratioNPsoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF

  IF (icycle >1) THEN
    casapool%nplant     = MAX(1.e-6,casapool%nplant)
    casapool%nlitter    = MAX(1.e-6,casapool%nlitter)
    casapool%nsoil      = MAX(1.e-6,casapool%nsoil)
    casapool%nsoilmin   = MAX(1.e-6,casapool%nsoilmin)
    casabal%nplantlast  = casapool%nplant
    casabal%nlitterlast = casapool%nlitter
    casabal%nsoillast   = casapool%nsoil
    casabal%nsoilminlast= casapool%nsoilmin
    casabal%sumnbal     = 0.0
    casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
    casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
  ENDIF

  IF (icycle >2) THEN
    casapool%pplant       = MAX(1.0e-7,casapool%pplant)
    casapool%plitter      = MAX(1.0e-7,casapool%plitter)
    casapool%psoil        = MAX(1.0e-7,casapool%psoil)
    casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0,  YP
    casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
    casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
    casabal%pplantlast    = casapool%pplant
    casabal%plitterlast   = casapool%plitter
    casabal%psoillast     = casapool%psoil
    casabal%psoillablast  = casapool%psoillab
    casabal%psoilsorblast = casapool%psoilsorb
    casabal%psoilocclast  = casapool%psoilocc
    casabal%sumpbal       = 0.0
    casabal%FPweayear=0.0;casabal%FPdustyear=0.0;casabal%FPsnetyear=0.0
    casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
  ENDIF

!  npt=26493
!  print *, 'pool sizes for 26493 after: ', casamet%glai(npt),casapool%cplant(npt,:),  &
!            casapool%nplant(npt,:), casapool%pplant(npt,:),                           &
!            casapool%clitter(npt,:), casapool%nlitter(npt,:),casapool%plitter(npt,:), &
!            casapool%csoil(npt,:),   casapool%nsoil(npt,:), casapool%psoil(npt,:),    &
!            casapool%nsoilmin(npt), casapool%psoillab(npt),casapool%psoilsorb(npt), casapool%psoilocc(npt)
END SUBROUTINE get_casa_restart


