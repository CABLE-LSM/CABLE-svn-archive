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
! Purpose: bgcdriver - interface between casacnp and cable
!          sumcflux  - accumulating carbon fluxes (not required for UM)
!
! Called from: cable_driver for offline version
!              Not currently called/available for ACCESS version
!
! Contact: Yingping.Wang@csiro.au
!
! History: Model development by Yingping Wang, coupling to Mk3L by Bernard Pak
!          ssoil changed to ssnow
!
! Sep 2015 Vanessa Haverd: Modifications to optionally call POP annually; replace default wood turnover
! with woody biomass mortality; track sapwood biomass for use in casacnp autotrophic
! respiration; track sapwood cross-sectional area for use in carbon allocation when LALLOC=3
! ==============================================================================

!#     define UM_BUILD YES 
SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssnow,canopy,veg,soil, &
                     casabiome,casapool,casaflux,casamet,casabal,phen, &
                     pop, spinConv, spinup, ktauday, idoy,loy, dump_read,   &
                     dump_write, gpp_ann )

   USE cable_def_types_mod
   USE cable_common_module, only: cable_runtime
   USE casadimension
   USE casaparm
   USE casavariable
   USE phenvariable
   USE cable_common_module,  ONLY: CurYear, CABLE_USER
   USE TypeDef,              ONLY: i4b, dp
   USE POPMODULE,            ONLY: POPStep
   USE POP_TYPES,            ONLY: POP_TYPE

   IMPLICIT NONE
 
   INTEGER,      INTENT(IN) :: ktau ! integration step number
   INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
   INTEGER,      INTENT(IN) :: kend ! total # timesteps in run
   
   INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
   INTEGER,      INTENT(IN)                  :: ktauday
   logical,      INTENT(IN) :: spinConv, spinup
   logical,      INTENT(IN) :: dump_read, dump_write 
        
   REAL,         INTENT(IN) :: dels ! time setp size (s)
   REAL, DIMENSION(mp)    ,    INTENT(IN) :: gpp_ann
   TYPE (met_type), INTENT(INOUT)       :: met  ! met input variables
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow ! soil and snow variables
   TYPE (canopy_type), INTENT(INOUT) :: canopy ! vegetation variables
   TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
   TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
   TYPE (casa_biome),          INTENT(INOUT) :: casabiome
   TYPE (casa_pool),           INTENT(INOUT) :: casapool
   TYPE (casa_flux),           INTENT(INOUT) :: casaflux
   TYPE (casa_met),            INTENT(INOUT) :: casamet
   TYPE (casa_balance),        INTENT(INOUT) :: casabal
   TYPE (phen_variable),       INTENT(INOUT) :: phen
   TYPE(POP_TYPE),             INTENT(INOUT) :: POP
   !! vh_js !! LALLOC declared and value assigned here instead of in casa_variable
   INTEGER :: LALLOC 

   ! local variables added ypwang 5/nov/2012
   real,      dimension(mp)  :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
   real,      dimension(mp)  :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
   real,      dimension(mp)  :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
   real(r_2), dimension(mp)  :: xnplimit,  xkNlimiting, xklitter, xksoil ,xkleaf,xkleafcold,xkleafdry

   INTEGER                                   :: it, nit
   REAL(dp)                               :: StemNPP(mp,2)
   REAL(dp), allocatable, save ::  LAImax(:)    , Cleafmean(:),  Crootmean(:)
   REAL(dp), allocatable :: NPPtoGPP(:)
   CHARACTER                                 :: cyear*4
   CHARACTER                                 :: ncfile*99
   INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles
   
         
   if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
   if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
   if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
   if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
   if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))
      

   !! vh_js !!
    IF (cable_user%CALL_POP) THEN
       LALLOC = 3 ! POP setting (requires casaflux%sapwood_area)
       Iw = POP%Iwood
    ELSE
       LALLOC = 0  ! default setting
    ENDIF
      
   
   IF ( .NOT. dump_read ) THEN  ! construct casa met and flux inputs from current CABLE run
      IF(ktau == kstart) THEN
         casamet%tairk  = 0.0
         casamet%tsoil  = 0.0
         casamet%moist  = 0.0
         casaflux%cgpp  = 0.0
         ! add initializations (BP jul2010)
         casaflux%Crsoil   = 0.0
         casaflux%crgplant = 0.0
         casaflux%crmplant = 0.0
         casaflux%clabloss = 0.0
         ! casaflux%crmplant(:,leaf) = 0.0
         ! end changes (BP jul2010)

      ENDIF

      IF(MOD(ktau,ktauday)==1) THEN
         casamet%tairk = met%tk
         casamet%tsoil = ssnow%tgg
         casamet%moist = ssnow%wb
         casaflux%cgpp = (-canopy%fpn+canopy%frday)*dels
         casaflux%crmplant(:,leaf) = canopy%frday*dels
      ELSE
         Casamet%tairk  =casamet%tairk + met%tk
         casamet%tsoil = casamet%tsoil + ssnow%tgg
         casamet%moist = casamet%moist + ssnow%wb
         casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dels
         casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + canopy%frday*dels
      ENDIF

     
      IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)
   
         IF ( icycle .GT. 0 ) THEN

            ! if calling POP, set woody biomass turnover to zero: it will be set using POP mortality
            IF (cable_user%CALL_POP)  casabiome%plantrate(:,2) = 0.0

            CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd, gpp_ann)

            IF (cable_user%CALL_POP) THEN ! CALL_POP

               ! accumulate annual variables for use in POP
               IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                  casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
                  LAImax = casamet%glai
                  Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
                  Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
               ELSE
                  casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                  LAImax = max(casamet%glai, LAImax)
                  Cleafmean = Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
                  Crootmean = Crootmean +casapool%cplant(:,3)/real(LOY)/1000.
               ENDIF
        
                    
               IF(MOD((ktau-kstart+1)/ktauday,LOY)==0) THEN ! end of year

                  StemNPP(:,1) = casaflux%stemnpp !/float(ktauday*LOY)
                  StemNPP(:,2) = 0.0
                  WHERE (casabal%FCgppyear > 1.e-5)
                     NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
                  ELSEWHERE
                     NPPtoGPP = 0.5
                  ENDWHERE
                  CALL POPStep(pop, max(StemNPP(Iw,:)/1000.,0.01), int(veg%disturbance_interval(Iw,:), i4b),&
                   real(veg%disturbance_intensity(Iw,:),dp)      ,&
                   LAImax(Iw), Cleafmean(Iw), Crootmean(Iw), NPPtoGPP(Iw))

                  WHERE (pop%pop_grid(:)%cmass_sum_old.gt.1.e-12)       
                     casapool%CLitter(Iw,3) = casapool%CLitter(Iw,3) + &
                          (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality+POP%pop_grid(:)%cat_mortality &
                          + POP%pop_grid(:)%fire_mortality + POP%pop_grid(:)%res_mortality  ) * &
                          casapool%Cplant(Iw,2)/POP%pop_grid(:)%cmass_sum_old 

                     casapool%Cplant(Iw,2) = casapool%Cplant(Iw,2) -  &
                          (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality+POP%pop_grid(:)%cat_mortality &
                          + POP%pop_grid(:)%fire_mortality+ POP%pop_grid(:)%res_mortality  ) * &
                          casapool%Cplant(:,2)/POP%pop_grid(:)%cmass_sum_old 
                     
                     casaflux%frac_sapwood(Iw) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
                     casaflux%sapwood_area(Iw) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)       
                     
                  ENDWHERE
               
               ENDIF  ! end of year
            ENDIF ! CALL_POP
               
         ENDIF  ! icycle .gt. 0
   
         IF( (.NOT.spinup).OR.(spinup.AND.spinConv)) THEN 
            IF ( dump_write )  THEN
               !CLN CHECK FOR LEAP YEAR
               WRITE(CYEAR,FMT="(I4)") CurYear + INT((ktau-kstart)/(LOY*ktauday))
               ncfile = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
               CALL write_casa_dump( ncfile, casamet , casaflux, idoy, &
                                 kend/ktauday )

            ENDIF
         ENDIF

      ENDIF  ! end of day

   ELSE ! dump_read: ! use casa met and flux inputs from dumpfile

      IF( MOD((ktau-kstart+1),ktauday) == 0 ) THEN  ! end of day

         IF (cable_user%CALL_POP)  casabiome%plantrate(:,2) = 0.0
       
         CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
              casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
              cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
              nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
              pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd, gpp_ann)

         IF (cable_user%CALL_POP) THEN

            ! accumulate annual variables for use in POP
            IF(MOD(ktau/ktauday,LOY)==1) THEN
               casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 ! (assumes 70% of wood NPP is allocated above ground)
               LAImax = casamet%glai
               Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
               Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
            ELSE
               casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
               LAImax = max(casamet%glai, LAImax)
               Cleafmean = Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
               Crootmean = Crootmean +casapool%cplant(:,3)/real(LOY)/1000.
            ENDIF
        
       
            IF(MOD((ktau-kstart+1)/ktauday,LOY)==0) THEN  ! end of year
               StemNPP(:,1) = casaflux%stemnpp ! (assumes 70% of wood NPP is allocated above ground & static alloc)
               StemNPP(:,2) = 0.0
               WHERE (casabal%FCgppyear > 1.e-5)
                  NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
               ELSEWHERE
                  NPPtoGPP = 0.5
               ENDWHERE

               CALL POPStep(pop, max(StemNPP(Iw,:)/1000.,0.01), int(veg%disturbance_interval(Iw,:), i4b),&
                    real(veg%disturbance_intensity(Iw,:),dp)      ,&
                    LAImax(Iw), Cleafmean(Iw), Crootmean(Iw), NPPtoGPP(Iw))
               
               WHERE (pop%pop_grid(:)%cmass_sum_old.gt.1.e-12)       
                  casapool%CLitter(Iw,3) = casapool%CLitter(Iw,3) + &
                       (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality+POP%pop_grid(:)%cat_mortality &
                       + POP%pop_grid(:)%fire_mortality + POP%pop_grid(:)%res_mortality  ) * &
                       casapool%Cplant(Iw,2)/POP%pop_grid(:)%cmass_sum_old 
                  
                  casapool%Cplant(Iw,2) = casapool%Cplant(Iw,2) -  &
                       (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality+POP%pop_grid(:)%cat_mortality &
                          + POP%pop_grid(:)%fire_mortality+ POP%pop_grid(:)%res_mortality  ) * &
                          casapool%Cplant(:,2)/POP%pop_grid(:)%cmass_sum_old 
                  
                  casaflux%frac_sapwood(Iw) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
                  casaflux%sapwood_area(Iw) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)       
                  
               ENDWHERE
               
               ENDIF ! end of year

         ENDIF ! CALL_POP

      ENDIF ! end of day

   ENDIF ! dump_read

END SUBROUTINE bgcdriver

SUBROUTINE read_casa_dump(  ncfile, casamet, casaflux, ncall, kend, allATonce )
      USE netcdf
      USE cable_def_types_mod,   ONLY : r_2,ms,mp
      USE casadimension,         ONLY : mplant,mdyear
      USE casavariable,          ONLY : casa_met, casa_flux
#     ifndef UM_BUILD            
      USE cable_diag_module,     ONLY : get_var_ncr2, &
                                        get_var_ncr3, stderr_nc
#     endif
      IMPLICIT NONE

      TYPE (casa_flux), INTENT(INOUT) :: casaflux
      TYPE (casa_met), INTENT(inout)  :: casamet
      INTEGER, INTENT(in)             :: kend, ncall
      CHARACTER(len=*), INTENT(in)    :: ncfile
      LOGICAL, INTENT(in)             :: allATonce

      !netcdf IDs/ names
      INTEGER, PARAMETER :: num_vars=7 
      INTEGER, PARAMETER :: num_dims=3 
      INTEGER, SAVE                        :: ncrid  ! netcdf file ID
      INTEGER , DIMENSION(num_vars)        :: varrID ! (1) tvair, (2) pmb

      !vars 
      CHARACTER(len=*), DIMENSION(num_vars), PARAMETER :: &
            var_name =  (/  "lat          ", &
                            "lon          ", &
                            "casamet_tairk", & 
                            "tsoil        ", &
                            "moist        ", &
                            "cgpp         ", &
                            "crmplant     " /)

      REAL     , DIMENSION(mp)        :: lat, lon
      REAL(r_2), DIMENSION(mp)        :: tairk,  cgpp
      REAL(r_2), DIMENSION(mp,ms)     :: tsoil, moist
      REAL(r_2), DIMENSION(mp,mplant) :: crmplant
      INTEGER :: ncok,  idoy

#     ifndef UM_BUILD            
      IF ( allATonce .OR. ncall .EQ. 1 ) THEN
         ncok = NF90_OPEN(TRIM(ncfile), nf90_nowrite, ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'re-opening ', ncfile)
      ENDIF
#     endif
      IF ( allATonce ) THEN
         DO idoy=1,mdyear
            
            CALL get_var_ncr2(ncrid, var_name(3), tairk   , idoy )
            CALL get_var_ncr3(ncrid, var_name(4), tsoil   , idoy ,ms)
            CALL get_var_ncr3(ncrid, var_name(5), moist   , idoy ,ms)
            CALL get_var_ncr2(ncrid, var_name(6), cgpp    , idoy )
            CALL get_var_ncr3(ncrid, var_name(7), crmplant, idoy ,mplant)
            
            casamet%Tairkspin(:,idoy) = tairk
            casamet%cgppspin (:,idoy) = cgpp
            casamet%crmplantspin_1(:,idoy) = crmplant(:,1)
            casamet%crmplantspin_2(:,idoy) = crmplant(:,2)
            casamet%crmplantspin_3(:,idoy) = crmplant(:,3)
            casamet%Tsoilspin_1(:,idoy)    = tsoil(:,1)
            casamet%Tsoilspin_2(:,idoy)    = tsoil(:,2)
            casamet%Tsoilspin_3(:,idoy)    = tsoil(:,3)
            casamet%Tsoilspin_4(:,idoy)    = tsoil(:,4)
            casamet%Tsoilspin_5(:,idoy)    = tsoil(:,5)
            casamet%Tsoilspin_6(:,idoy)    = tsoil(:,6)
            casamet%moistspin_1(:,idoy)    = moist(:,1)
            casamet%moistspin_2(:,idoy)    = moist(:,2)
            casamet%moistspin_3(:,idoy)    = moist(:,3)
            casamet%moistspin_4(:,idoy)    = moist(:,4)
            casamet%moistspin_5(:,idoy)    = moist(:,5)
            casamet%moistspin_6(:,idoy)    = moist(:,6)
            
         END DO
      ELSE

         CALL get_var_ncr2(ncrid, var_name(3), tairk   ,ncall )
         CALL get_var_ncr3(ncrid, var_name(4), tsoil   ,ncall , ms)
         CALL get_var_ncr3(ncrid, var_name(5), moist   ,ncall , ms)
         CALL get_var_ncr2(ncrid, var_name(6), cgpp    ,ncall )
         CALL get_var_ncr3(ncrid, var_name(7), crmplant,ncall , mplant)
         
         casamet%tairk     = tairk
         casamet%tsoil     = tsoil
         casamet%moist     = moist
         casaflux%cgpp     = cgpp
         casaflux%crmplant = crmplant

      ENDIF

#     ifndef UM_BUILD            
      IF ( allATonce .OR. ncall .EQ. kend ) THEN
         ncok = NF90_CLOSE(ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'closing ', ncfile)
      ENDIF
#     endif
   END SUBROUTINE read_casa_dump

!! DOES THIS NEED TO BE DELETED FOR NOW - REPLACED WITH BP CODE (LATER?)

SUBROUTINE write_casa_dump( ncfile, casamet, casaflux, n_call, kend )
  USE netcdf
  USE cable_def_types_mod,   ONLY : r_2,ms,mp
  USE cable_common_module,   ONLY : kend_gl
#     ifndef UM_BUILD            
  USE cable_diag_module,     ONLY : def_dims, def_vars, def_var_atts, &
       put_var_ncr1, put_var_ncr2,       &
       put_var_ncr3, stderr_nc
#     endif
  USE casavariable,          ONLY : CASA_MET, CASA_FLUX
  USE casadimension,         ONLY : mplant

  IMPLICIT NONE
     
  INTEGER, INTENT(in) :: &
       n_call, &         ! this timestep #
       kend              ! final timestep of run


  !number of instances. dummied here and so=1
  !integer :: inst =1

  !netcdf IDs/ names
  CHARACTER(len=*)   :: ncfile
  INTEGER, PARAMETER :: num_vars=7
  INTEGER, PARAMETER :: num_dims=3
  INTEGER, SAVE :: ncid       ! netcdf file ID

  !vars
  CHARACTER(len=*), DIMENSION(num_vars), PARAMETER :: &
       var_name =  (/  "lat          ", &
       "lon          ", &
       "casamet_tairk", &
       "tsoil        ", &
       "moist        ", &
       "cgpp         ", &
       "crmplant     " /)

  INTEGER, DIMENSION(num_vars) :: varID ! (1) tvair, (2) pmb

  !dims
  CHARACTER(len=*), DIMENSION(num_dims), PARAMETER :: &
       dim_name =  (/ "pnt ", &
       "soil", &
       "time" /)

  INTEGER, PARAMETER :: soil_dim = 6

  INTEGER, DIMENSION(soil_dim), PARAMETER  :: soil = (/ 1,2,3,4,5,6 /)

  INTEGER, DIMENSION(num_dims)  :: &
       dimID   ! (1) x, (2) y, (3) time

  INTEGER, DIMENSION(num_dims)  :: &
                                !x,y generally lat/lon BUT for single site = 1,1
       dim_len = (/-1,soil_dim,-1/)  ! (1) x, (2) y, (3) soil, (4) time [re-set]

  TYPE(CASA_MET)  :: CASAMET
  TYPE(CASA_FLUX) :: CASAFLUX

  !local only
  INTEGER :: ncok      !ncdf return status

  ! END header
#ifndef UM_BUILD 
  dim_len(1)        = mp
  dim_len(num_dims) = kend

  IF (n_call == 1) THEN
     ! create netCDF dataset: enter define mode
     ncok = nf90_create(path = TRIM(ncfile), cmode = nf90_clobber, ncid = ncid)
     IF (ncok /= nf90_noerr) CALL stderr_nc(ncok,'ncdf creating ', ncfile)

     ! define dimensions: from name and length
     CALL def_dims(num_dims, ncid, dimID, dim_len, dim_name )

     ! define variables: from name, type, dims
     CALL def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

     ! define variable attributes
     !CLN LATER!             CALL def_var_atts( ncfile, ncid, varID )

     ncok = nf90_enddef(ncid)

     CALL put_var_ncr1(ncid, var_name(1), REAL(casamet%lat)  )
     CALL put_var_ncr1(ncid, var_name(2), REAL(casamet%lon)  )

  ENDIF

  CALL put_var_ncr2(ncid, var_name(3), casamet%tairk    ,n_call )
  CALL put_var_ncr3(ncid, var_name(4), casamet%tsoil    ,n_call, ms )
  CALL put_var_ncr3(ncid, var_name(5), casamet%moist    ,n_call, ms )
  CALL put_var_ncr2(ncid, var_name(6), casaflux%cgpp    ,n_call )
  CALL put_var_ncr3(ncid, var_name(7), casaflux%crmplant,n_call, mplant )

  IF (n_call == kend ) &
       ncok = nf90_close(ncid)            ! close: save new netCDF dataset
#endif
END SUBROUTINE write_casa_dump

 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx  ! local variables
  real, dimension(17)                   ::  xnslope
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  ! first initialize 
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf) 
  npleafx = 14.2 

  DO np=1,mp
    ivt=veg%iveg(np)
    IF (casamet%iveg2(np)/=icewater & 
        .AND. casamet%glai(np)>casabiome%glaimin(ivt)  &
        .AND. casapool%cplant(np,leaf)>0.0) THEN
      ncleafx(np) = MIN(casabiome%ratioNCplantmax(ivt,leaf), &
                        MAX(casabiome%ratioNCplantmin(ivt,leaf), &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf)))
      IF (icycle>2 .AND. casapool%pplant(np,leaf)>0.0) THEN
        npleafx(np) = MIN(30.0,MAX(8.0,real(casapool%nplant(np,leaf) &
                /casapool%pplant(np,leaf))))
      ENDIF
    ENDIF

    IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
      IF (ivt/=2) THEN
        veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                        + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
      ELSE
        IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
          veg%vcmax(np) = ( casabiome%nintercept(ivt)  &
                          + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ELSE
          veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                          + casabiome%nslope(ivt)*ncleafx(np)/casabiome%sla(ivt) )*1.0e-6
        ENDIF
      ENDIF
      veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)
    ENDIF

!    veg%vcmax(np) = ( nintercept(ivt)  &
!                  + nslope(ivt)*(0.4+8.5/npleafx(np)) &
!                  * ncleafx(np)/casabiome%sla(ivt))*(1.0e-6)
!    veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)

!write(*,991) np, ivt,veg%vlai(np),veg%vcmax(np)*1.0e6
!write(*,891) np,ivt,casapool%cplant(np,leaf),casapool%nplant(np,leaf),casapool%pplant(np,leaf)
!891 format(2(i6),3(f9.3,2x))
  ENDDO

  veg%ejmax = 2.0 * veg%vcmax
!991 format(i6,2x,i4,2x,2(f9.3,2x))
 END SUBROUTINE casa_feedback


SUBROUTINE sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
     soil, ssnow, sum_flux, veg, met, casaflux, l_vcmaxFeedbk)

  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ktau ! integration step number
  INTEGER, INTENT(IN)    :: kstart ! starting value of ktau
  INTEGER, INTENT(IN)    :: kend ! total # timesteps in run
!  INTEGER, INTENT(IN)    :: mvtype  ! Number of veg types
!  INTEGER, INTENT(IN)    :: mstype ! Number of soil types
  REAL,    INTENT(IN)    :: dels ! time setp size (s)
  TYPE (bgc_pool_type),       INTENT(INOUT) :: bgc
  TYPE (canopy_type),         INTENT(INOUT) :: canopy
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil
  TYPE (soil_snow_type),      INTENT(INOUT) :: ssnow
  TYPE (sum_flux_type),       INTENT(INOUT) :: sum_flux
  TYPE (met_type),            INTENT(IN)    :: met    
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  LOGICAL, INTENT(IN)   :: l_vcmaxFeedbk ! using prognostic Vcmax

!   if(icycle<=0) then
!     these are executed in cbm
!      CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)
!      CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)
!   else
    if(icycle>0) then
       canopy%frp(:) = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,froot) &
                        +casaflux%crgplant(:))/86400.0
       canopy%frs(:) = casaflux%Crsoil(:)/86400.0
       canopy%frpw(:)= casaflux%crmplant(:,wood)/86400.0
       canopy%frpr(:)= casaflux%crmplant(:,froot)/86400.0
    endif
    if(ktau == kstart) then
       sum_flux%sumpn  = canopy%fpn*dels
       sum_flux%sumrd  = canopy%frday*dels
       sum_flux%dsumpn = canopy%fpn*dels
       sum_flux%dsumrd = canopy%frday*dels
       sum_flux%sumrpw = canopy%frpw*dels
       sum_flux%sumrpr = canopy%frpr*dels
       sum_flux%sumrp  = canopy%frp*dels
       sum_flux%dsumrp = canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = canopy%frs*dels
    else
       sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dels
       sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dels
       sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dels
       sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dels
       sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dels
       sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dels
       sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dels
       sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = sum_flux%sumrs+canopy%frs*dels
    endif
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp
    IF (icycle <= 1) THEN
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    ELSE
      IF (l_vcmaxFeedbk) THEN
        canopy%fnee = canopy%fpn + canopy%frs + canopy%frp &
                    + casaflux%clabloss(:)/86400.0
      ELSE
        canopy%fnee = (casaflux%Crsoil-casaflux%cnpp+casaflux%clabloss)/86400.0
      ENDIF
    ENDIF

END SUBROUTINE sumcflux

  SUBROUTINE totcnppools(kloop,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
                               bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)
  ! this subroutine is temporary, and its needs to be modified for multiple tiles within a cell
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,                     INTENT(IN)  :: kloop
  TYPE (veg_parameter_type),   INTENT(IN)  :: veg  ! vegetation parameters
  TYPE(casa_pool),             INTENT(IN)  :: casapool
  TYPE(casa_met),              INTENT(IN)  :: casamet
  real,      dimension(5,mvtype,mplant)    :: bmcplant,  bmnplant,  bmpplant
  real,      dimension(5,mvtype,mlitter)   :: bmclitter, bmnlitter, bmplitter
  real,      dimension(5,mvtype,msoil)     :: bmcsoil,   bmnsoil,   bmpsoil
  real,      dimension(5,mvtype)           :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  real,      dimension(mvtype)             :: bmarea
  ! local variables
  INTEGER  npt,nvt


      bmcplant(kloop,:,:)  = 0.0;  bmnplant(kloop,:,:)  = 0.0; bmpplant(kloop,:,:)  = 0.0
      bmclitter(kloop,:,:) = 0.0;  bmnlitter(kloop,:,:) = 0.0; bmplitter(kloop,:,:) = 0.0
      bmcsoil(kloop,:,:)   = 0.0;  bmnsoil(kloop,:,:)   = 0.0; bmpsoil(kloop,:,:)   = 0.0
      bmnsoilmin(kloop,:)  = 0.0;  bmpsoillab(kloop,:)  = 0.0; bmpsoilsorb(kloop,:) = 0.0;  bmpsoilocc(kloop,:) = 0.0

      bmarea(:) = 0.0

      do npt=1,mp
         nvt=veg%iveg(npt)
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)   + casapool%cplant(npt,:) * casamet%areacell(npt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)   + casapool%nplant(npt,:) * casamet%areacell(npt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)   + casapool%pplant(npt,:) * casamet%areacell(npt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:) + casapool%clitter(npt,:) * casamet%areacell(npt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:) + casapool%nlitter(npt,:) * casamet%areacell(npt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:) + casapool%plitter(npt,:) * casamet%areacell(npt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)     + casapool%csoil(npt,:) * casamet%areacell(npt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)     + casapool%nsoil(npt,:) * casamet%areacell(npt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)     + casapool%psoil(npt,:) * casamet%areacell(npt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)   + casapool%nsoilmin(npt) * casamet%areacell(npt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)   + casapool%psoillab(npt) * casamet%areacell(npt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)  + casapool%psoilsorb(npt) * casamet%areacell(npt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)   + casapool%psoilocc(npt) * casamet%areacell(npt)
         bmarea(nvt)  = bmarea(nvt) + casamet%areacell(npt)
      enddo

      do nvt=1,mvtype
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)/bmarea(nvt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)/bmarea(nvt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)/bmarea(nvt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:)/bmarea(nvt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:)/bmarea(nvt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:)/bmarea(nvt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)/bmarea(nvt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)/bmarea(nvt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)/bmarea(nvt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)/bmarea(nvt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)/bmarea(nvt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)/bmarea(nvt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)/bmarea(nvt)
      enddo

  END SUBROUTINE totcnppools

  SUBROUTINE analyticpool(kend,veg,soil,casabiome,casapool,                                 &
                          casaflux,casamet,casabal,phen,                                    &
                          avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd, &
                          avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd, &
                          avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd, &
                          avgcgpp, avgcnpp, avgnuptake, avgpuptake,                         &
                          avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil,                 &
                          avgratioNCsoilmic,avgratioNCsoilslow,avgratioNCsoilpass,         &
                          avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc)
  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kend
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen


  ! local variables
  real, dimension(mp)      :: avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd
  real, dimension(mp)      :: avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd
  real, dimension(mp)      :: avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd
  real, dimension(mp)      :: avgcgpp, avgcnpp, avgnuptake, avgpuptake
  real(r_2), dimension(mp) :: avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil
  real, dimension(mp)      :: avgratioNCsoilmic,  avgratioNCsoilslow,  avgratioNCsoilpass
  real,      dimension(mp) :: avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol     1       61.3
  ! Andisol     2       103.9
  ! Aridisol    3       92.8
  ! Entisol     4       136.9
  ! Gellisol    5       98.2
  ! Histosol    6       107.6
  ! Inceptisol  7       84.1
  ! Mollisol    8       110.1
  ! Oxisol      9       35.4
  ! Spodosol    10      41.0
  ! Ultisol     11      51.5
  ! Vertisol    12      190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

  INTEGER :: nyear,iyear
  real year

  ! compute the mean litter input in g(C, N and P)/day from plant pools
  casaflux%fromLtoS = 0.0
  casaflux%fromStoS = 0.0

  casabal%sumcbal(:)   = 0.0
  casabal%sumnbal(:)   = 0.0
  casabal%sumpbal(:)   = 0.0

  do npt=1,mp
  if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) THEN
    casaflux%fromLtoS(npt,mic,metb)   = 0.45
                                          ! metb -> mic
    casaflux%fromLtoS(npt,mic,str)   = 0.45*(1.0-casabiome%fracLigninplant(veg%iveg(npt),leaf))
                                          ! str -> mic
    casaflux%fromLtoS(npt,slow,str)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),leaf)
                                          ! str -> slow
    casaflux%fromLtoS(npt,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood))
                                          ! CWD -> fmic
    casaflux%fromLtoS(npt,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood)
                                          ! CWD -> slow
!! set the following two backflow to set (see Bolker 199x)
!    casaflux%fromStoS(npt,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(npt))
!    casaflux%fromStoS(npt,mic,pass)  = 0.45

    casaflux%fromStoS(npt,slow,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.997 - 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.003 + 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(npt) )

    casaflux%klitter(npt,metb) = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),metb)
    casaflux%klitter(npt,str)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),str)&
                               * exp(-3.0*casabiome%fracLigninplant(veg%iveg(npt),leaf))
    casaflux%klitter(npt,cwd)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),cwd)

    casaflux%ksoil(npt,mic)    = avgxksoil(npt)*casabiome%soilrate(veg%iveg(npt),mic)  &
                               * (1.0 - 0.75 *(soil%silt(npt)+soil%clay(npt)))
    casaflux%ksoil(npt,slow)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),slow)
    casaflux%ksoil(npt,pass)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),pass)

    if(veg%iveg(npt)==cropland) THEN     ! for cultivated land type
       casaflux%ksoil(npt,mic)  = casaflux%ksoil(npt,mic) * 1.25
       casaflux%ksoil(npt,slow) = casaflux%ksoil(npt,slow)* 1.5
       casaflux%ksoil(npt,pass) = casaflux%ksoil(npt,pass)* 1.5
    endif
  endif
  enddo


  do npt=1,mp       
   if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) then
    ! compute steady-state litter and soil C pool sizes
     casapool%clitter(npt,metb) = (avgcleaf2met(npt)+avgcroot2met(npt))/casaflux%klitter(npt,metb)
     casapool%clitter(npt,str) = (avgcleaf2str(npt)+avgcroot2str(npt))/casaflux%klitter(npt,str)
     casapool%clitter(npt,cwd) = (avgcwood2cwd(npt))/casaflux%klitter(npt,cwd)
     casapool%csoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                 +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                 +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                               /casaflux%ksoil(npt,mic)
      casapool%csoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                 + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                 + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                 + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                /casaflux%ksoil(npt,slow)
      casapool%csoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                  +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                /casaflux%ksoil(npt,pass)
      if(icycle <=1) then
         casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) * casapool%clitter(npt,:)
         casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   * casapool%Csoil(npt,:)
         casapool%nsoilmin(npt) = 2.0
         casabal%sumnbal(npt)   = 0.0
       ELSE
       ! compute steady-state litter and soil N pool sizes
         casapool%nlitter(npt,metb) = (avgnleaf2met(npt)+avgnroot2met(npt))/casaflux%klitter(npt,metb)
         casapool%nlitter(npt,str) = (avgnleaf2str(npt)+avgnroot2str(npt))/casaflux%klitter(npt,str)
         casapool%nlitter(npt,cwd) = (avgnwood2cwd(npt))/casaflux%klitter(npt,cwd)

         casapool%nsoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * avgratioNCsoilmic(npt)/casaflux%ksoil(npt,mic)
         casapool%nsoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * avgratioNCsoilslow(npt)/casaflux%ksoil(npt,slow)
         casapool%nsoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   * avgratioNCsoilpass(npt)/casaflux%ksoil(npt,pass)
          casapool%Nsoilmin(npt)    = avgnsoilmin(npt)

        ENDIF

        IF (icycle<=2) THEN
            totpsoil(npt)          = psorder(casamet%isorder(npt)) *xpsoil50(casamet%isorder(npt))
           casapool%plitter(npt,:)= casapool%Nlitter(npt,:)/casapool%ratioNPlitter(npt,:)
            casapool%psoil(npt,:)  = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
            ! why is this commented here but used in UM
            ! casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:)  * casapool%clitter(npt,:)
            ! casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)    * casapool%Csoil(npt,:)
            casapool%psoillab(npt) = totpsoil(npt) *fracpLab(casamet%isorder(npt))
            casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                    /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
            casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(casamet%isorder(npt))
        ELSE
        ! compute the steady-state litter and soil P pools
          casapool%plitter(npt,metb) = (avgpleaf2met(npt)+avgproot2met(npt))/casaflux%klitter(npt,metb)
          casapool%plitter(npt,str) = (avgpleaf2str(npt)+avgproot2str(npt))/casaflux%klitter(npt,str)
          casapool%plitter(npt,cwd) = (avgpwood2cwd(npt))/casaflux%klitter(npt,cwd)

          casapool%psoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * (casapool%ratioNCsoil(npt,mic)/casapool%ratioNPsoil(npt,mic))/casaflux%ksoil(npt,mic)
          casapool%psoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * (casapool%ratioNCsoil(npt,slow)/casapool%ratioNPsoil(npt,slow))/casaflux%ksoil(npt,slow)
          casapool%psoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   *  (casapool%ratioNCsoil(npt,pass)/casapool%ratioNPsoil(npt,pass))/casaflux%ksoil(npt,pass)         
          ! assign the mineral pools
          casapool%psoillab(npt)      = avgpsoillab(npt)
          casapool%psoilsorb(npt)     = avgPsoilsorb(npt)
          casapool%psoilocc(npt)      = avgPsoilocc(npt)
        ENDIF
  ENDIF
  ENDDO

  END SUBROUTINE analyticpool


