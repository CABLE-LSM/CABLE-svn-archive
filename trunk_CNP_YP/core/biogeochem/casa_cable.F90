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
! ==============================================================================
!CABLE_LSM:This has to be commented for offline
!#define UM_BUILD YES

module casa_cable

contains

SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssnow,canopy,veg,soil, &
                     casabiome,casapool,casaflux,casamet,casabal,phen, &
                     spinConv, spinup, ktauday, idoy, dump_read, dump_write )

   USE cable_def_types_mod
   USE casadimension
   USE casaparm
   USE casavariable
   USE phenvariable
   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: ktau ! integration step number
   INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
   INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

   INTEGER,      INTENT(IN)                  :: idoy ! day of year (1-365)
   INTEGER,      INTENT(IN)                  :: ktauday
   logical,      INTENT(IN) :: spinConv, spinup
   logical,      INTENT(IN) :: dump_read, dump_write

   REAL,         INTENT(IN) :: dels ! time setp size (s)
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

   ! local variables added ypwang 5/nov/2012
   real,      dimension(mp)  :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
   real,      dimension(mp)  :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
   real,      dimension(mp)  :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
   real(r_2), dimension(mp)  :: xnplimit,  xkNlimiting, xklitter, xksoil ,xkleaf,xkleafcold,xkleafdry


   !    phen%phase = 2

   IF ( .NOT. dump_read ) then
      if(ktau == kstart) then
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
      IF(mod(ktau,ktauday)==1) THEN
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
         casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + &
                                       canopy%frday*dels
      ENDIF

      IF(mod((ktau-kstart+1),ktauday)==0) THEN

         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)

         ! added ypwang 5/nov/2012
         if(ktau/ktauday .le. 365)then
           casamet%Tairkspin     (:,idoy) = casamet%tairk(:)
           casamet%cgppspin      (:,idoy) = casaflux%cgpp(:)
           casamet%crmplantspin_1(:,idoy) = casaflux%crmplant(:,1)
           casamet%crmplantspin_2(:,idoy) = casaflux%crmplant(:,2)
           casamet%crmplantspin_3(:,idoy) = casaflux%crmplant(:,3)
           casamet%Tsoilspin_1   (:,idoy) = casamet%tsoil(:,1)
           casamet%Tsoilspin_2   (:,idoy) = casamet%tsoil(:,2)
           casamet%Tsoilspin_3   (:,idoy) = casamet%tsoil(:,3)
           casamet%Tsoilspin_4   (:,idoy) = casamet%tsoil(:,4)
           casamet%Tsoilspin_5   (:,idoy) = casamet%tsoil(:,5)
           casamet%Tsoilspin_6   (:,idoy) = casamet%tsoil(:,6)
           casamet%moistspin_1   (:,idoy) = casamet%moist(:,1)
           casamet%moistspin_2   (:,idoy) = casamet%moist(:,2)
           casamet%moistspin_3   (:,idoy) = casamet%moist(:,3)
           casamet%moistspin_4   (:,idoy) = casamet%moist(:,4)
           casamet%moistspin_5   (:,idoy) = casamet%moist(:,5)
           casamet%moistspin_6   (:,idoy) = casamet%moist(:,6)
         end if

         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                         casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,  &
                         xksoil,xkleaf,xkleafcold,xkleafdry,&
                         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,   &
                         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,   &
                         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

       ! modified ypwang 5/nov/2012
       !  CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
       !             casamet,casabal,phen)
       !  IF((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN
       !     IF ( dump_write ) &
       !        comment out by ypwang 5/nov/2012
       !        call ncdf_dump( casamet%tairk, casamet%tsoil, casamet%moist, &
       !                        casaflux%cgpp, casaflux%crmplant, idoy, &
       !                        kend/ktauday )
       !  ENDIF

      ENDIF

   ELSE



      IF( mod((ktau-kstart+1),ktauday) == 0 )  then
      ! modified yp wang 5/nov/2012

       CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                       casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,  &
                       xksoil,xkleaf,xkleafcold,xkleafdry,&
                       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,   &
                       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,   &
                       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

      !CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
      !                casamet,casabal,phen)

      ENDIF

   ENDIF

END SUBROUTINE bgcdriver
! ==============================================================================

SUBROUTINE POPdriver(casaflux,casabal,veg, POP)

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


  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_flux),           INTENT(IN) :: casaflux
  TYPE (casa_balance),        INTENT(IN) :: casabal
  TYPE(POP_TYPE),             INTENT(INOUT) :: POP

  INTEGER                                   :: it, nit
  REAL(dp)                               :: StemNPP(mp,2)
  REAL(dp), allocatable :: NPPtoGPP(:)
  REAL(dp), allocatable ::  LAImax(:)  , Cleafmean(:),  Crootmean(:)
  CHARACTER                                 :: cyear*4
  CHARACTER                                 :: ncfile*99
  !! vh_js !!
  INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles

  ! INTEGER, INTENT(IN) :: wlogn
  INTEGER , parameter :: wlogn=6

  if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
  if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
  if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
  if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))

  IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP
     Iw = POP%Iwood

     StemNPP(:,1) = casaflux%stemnpp
     StemNPP(:,2) = 0.0
     WHERE (casabal%FCgppyear > 1.e-5 .and. casabal%FCnppyear > 1.e-5  )
        NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
     ELSEWHERE
        NPPtoGPP = 0.5
     ENDWHERE
     LAImax = casabal%LAImax
     Cleafmean = casabal%cleafmean
     Crootmean = casabal%Crootmean

     CALL POPStep(pop, max(StemNPP(Iw,:)/1000.,0.0001), int(veg%disturbance_interval(Iw,:), i4b),&
          real(veg%disturbance_intensity(Iw,:),dp)      ,&
          max(LAImax(Iw),0.001), Cleafmean(Iw), Crootmean(Iw), NPPtoGPP(Iw))


  ENDIF ! CALL_POP

END SUBROUTINE POPdriver
! ==============================================================================
subroutine read_casa_dump(ncfile, casamet, casaflux, ktau, kend)
      use netcdf
      use cable_def_types_mod,   only : r_2,ms
      use casadimension,         only : mplant,mdyear
      USE casa_cnp_module
      use casa_dump_module,     only : get_var_nc, stderr_nc

      TYPE (casa_flux), intent(inout) :: casaflux
      TYPE (casa_met), intent(inout)  :: casamet
      integer, intent(in) :: kend, ktau

      !netcdf IDs/ names
      character(len=*)  ncfile
      integer, parameter :: num_vars=5
      integer, parameter :: num_dims=4
      integer:: ncid       ! netcdf file ID

      !vars
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/  "casamet_tairk", &
                            "tsoil        ", &
                            "moist        ", &
                            "cgpp         ", &
                            "crmplant     " /)

      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb

      real(r_2), dimension(mp) :: &
         tairk,  &
         cgpp

      real(r_2), dimension(mp,ms) :: &
         tsoil, &
         moist

      real(r_2), dimension(mp,mplant) :: &
         crmplant

!      write(89,*)'opening file'
      ncok = NF90_OPEN(ncfile, nf90_nowrite, ncid)
         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile)

      do idoy=1,mdyear

!         write(89,*)'get tairk'
         call get_var_nc(ncid, var_name(1), tairk   , idoy, kend )
!         write(89,*)'get tsoil'
         call get_var_nc(ncid, var_name(2), tsoil   , idoy, kend ,ms)
!         write(89,*)'get moist'
         call get_var_nc(ncid, var_name(3), moist   , idoy, kend ,ms)
!         write(89,*)'get cgpp'
         call get_var_nc(ncid, var_name(4), cgpp    , idoy, kend )
!         write(89,*)'get crmplant'
         call get_var_nc(ncid, var_name(5), crmplant, idoy, kend ,mplant)


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

      end do

      ncok = NF90_CLOSE(ncid)
         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile)


   end subroutine read_casa_dump


SUBROUTINE write_casa_dump( ncfile, casamet, casaflux, phen, climate, n_call, kend )
  USE netcdf
  USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
  USE cable_common_module,   ONLY : kend_gl
#     ifndef UM_BUILD
  USE cable_diag_module,     ONLY : def_dims, def_vars, def_var_atts, &
       put_var_ncr1, put_var_ncr2,       &
       put_var_ncr3, stderr_nc
#     endif
  USE casavariable,          ONLY : CASA_MET, CASA_FLUX
  USE casadimension,         ONLY : mplant
  USE phenvariable
  USE cable_common_module,  ONLY:  CABLE_USER

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
       n_call, &         ! this timestep #
       kend              ! final timestep of run

  TYPE (casa_flux),             INTENT(IN) :: casaflux
  TYPE (casa_met),              INTENT(IN) :: casamet
  TYPE (phen_variable),         INTENT(IN) :: phen
  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

  !number of instances. dummied here and so=1
  !integer :: inst =1

  !netcdf IDs/ names
  CHARACTER(len=*)   :: ncfile
  INTEGER            :: num_vars
  INTEGER, PARAMETER :: num_dims=3
  INTEGER, SAVE :: ncid       ! netcdf file ID

  !vars
  CHARACTER, DIMENSION(:), POINTER :: var_name*15


  INTEGER, DIMENSION(:), POINTER :: varID ! (1) tvair, (2) pmb

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
       dim_len = (/-1,soil_dim,-1/)  ! (1) mp, (2) soil, (3) time [re-set]



  !local only
  INTEGER :: ncok      !ncdf return status

  ! END header
#ifndef UM_BUILD
  dim_len(1)        = mp
  dim_len(num_dims) = NF90_unlimited

  num_vars = 14

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    num_vars=num_vars+1
  ENDIF

  allocate(var_name(num_vars))
  allocate(varID(num_vars))

  var_name =  (/"lat          ", &
                "lon          ", &
                "casamet_tairk", &
                "tsoil        ", &
                "moist        ", &
                "cgpp         ", &
                "crmplant     ", &
                "phenphase    ", &
                "phendoyphase1", &
                "phendoyphase2", &
                "phendoyphase3", &
                "phendoyphase4", &
                "Ndep         ", &
                "Pdep         "/)

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    var_name(num_vars)="mtemp"
  ENDIF

  IF (n_call == 1) THEN

     ! create netCDF dataset: enter define mode
     ncok = nf90_create(path = TRIM(ncfile), cmode = nf90_clobber, ncid = ncid)
     IF (ncok /= nf90_noerr) CALL stderr_nc(ncok,'ncdf creating ', ncfile)

     !ncok = nf90_redef(ncid)
     !if (ncok /= nf90_noerr) call stderr_nc(ncok,'enter def mode', ncfile)

     ! define dimensions: from name and length
     CALL def_dims(num_dims, ncid, dimID, dim_len, dim_name )

     ! define variables: from name, type, dims
     CALL def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

     ! define variable attributes
     !CLN LATER!             CALL def_var_atts( ncfile, ncid, varID )

     ncok = nf90_enddef(ncid)
     if (ncok /= nf90_noerr) call stderr_nc(ncok,'end def mode', ncfile)


     CALL put_var_ncr1(ncid, var_name(1), REAL(casamet%lat)  )
     CALL put_var_ncr1(ncid, var_name(2), REAL(casamet%lon)  )


  ENDIF


  CALL put_var_ncr2(ncid, var_name(3), casamet%tairk    ,n_call )
  CALL put_var_ncr3(ncid, var_name(4), casamet%tsoil    ,n_call, ms )
  CALL put_var_ncr3(ncid, var_name(5), casamet%moist    ,n_call, ms )
  CALL put_var_ncr2(ncid, var_name(6), casaflux%cgpp    ,n_call )
  CALL put_var_ncr3(ncid, var_name(7), casaflux%crmplant,n_call, mplant )
  CALL put_var_ncr2(ncid, var_name(8), real(phen%phase , r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(9), real(phen%doyphase(:,1), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(10), real(phen%doyphase(:,2), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(11), real(phen%doyphase(:,3), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(12), real(phen%doyphase(:,4), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(13), real(casaflux%Nmindep,r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(14), real(casaflux%Pdep,r_2)    ,n_call )
  if (cable_user%CALL_climate) then
     CALL put_var_ncr2(ncid, var_name(15), real(climate%mtemp_max,r_2)    ,n_call )
  endif

  deallocate(var_name)

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

  ! first initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  npleafx(:) = casabiome%ratioNPplantmin(veg%iveg(:),leaf)

  DO np=1,mp
    ivt=veg%iveg(np)
    IF (casamet%iveg2(np)/=icewater &
        .AND. casamet%glai(np)>casabiome%glaimin(ivt)  &
        .AND. casapool%cplant(np,leaf)>0.0) THEN
      ncleafx(np) = MIN(casabiome%ratioNCplantmax(ivt,leaf), &
                        MAX(casabiome%ratioNCplantmin(ivt,leaf), &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf)))
!      ! new equation for ncleafx, added by YPW Apr 2013
!      ncleafx(np) = (casapool%nplant(np,leaf)/casapool%cplant(np,leaf)) &
!                  * veg%extkn(np) * casamet%glai(np) &
!                  / (1.0 - exp(-veg%extkn(np)*casamet%glai(np)))

      IF (icycle>2 .AND. casapool%pplant(np,leaf)>0.0) THEN
        npleafx(np) = MIN(30.0,MAX(8.0,casapool%nplant(np,leaf) &
                                      /casapool%pplant(np,leaf)))
      ENDIF
    ENDIF

    IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
      IF (ivt/=2) THEN
        veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                        + casabiome%nslope(ivt)     &
                        * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
      ELSE
        IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
          veg%vcmax(np) = ( casabiome%nintercept(ivt)  &
                          + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ELSE
          veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                          + casabiome%nslope(ivt)     &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ENDIF
      ENDIF
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

!    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
!!   For prognostic Vcmax, and NEE should include clabloss under nutrient limitation
!!   Q.Zhang 12/09/2011
!    canopy%fnee(:) = canopy%fnee(:) + casaflux%clabloss(:)/86400.0
!    ! Q.Zhang 08/06/2011. return NEE from casaflux when N/NP mode is activated.
!    ! NPP of CABLE's output is "potential" NPP, not "real" C input to casacnp
!    ! To derive nutrient limited NPP from CABLE's standard output, use NEE+Crsoil
!!    if (icycle>1) then
!!     canopy%fnee = (casaflux%Crsoil-casaflux%cnpp)/86400.0
!!    else
!!     canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
!!    end if

!    write(*,101) ktau,casaflux%Crsoil(:)
!    write(*,101) ktau,dels,veg%vlai,veg%vcmax*1.0e6,casaflux%cgpp,canopy%fpn*1.0e6/12.0,canopy%frp*1.0e6/12.0,canopy%frs*1.0e6/12.0,canopy%fnee*1.0e6/12.0
!101  format(i6,2x,100(f12.5,2x))
!    if(ktau==kend) then
!       PRINT *, 'carbon fluxes'
!       PRINT *, 'sumpn', sum_flux%sumpn
!       PRINT *, 'sumrd', sum_flux%sumrd
!       PRINT *, 'sumrp', sum_flux%sumrp
!       PRINT *, 'sumrs', sum_flux%sumrs
!       PRINT *, 'npp', sum_flux%sumpn+sum_flux%sumrp
!       PRINT *, 'nee', sum_flux%sumpn+sum_flux%sumrp+sum_flux%sumrs
!     !  PRINT *, 'carbon pools', leaf,wood,froot
!     !  PRINT *,  casapool%cplant(1,2),casaflux%crmplant(1,wood),casaflux%Crsoil(1)
!     !  PRINT *, 'respiration rate'
!     !  PRINT *,  casabiome%rmplant(1,2)*365.0
!    endif

END SUBROUTINE sumcflux

call totcnppools(1,veg,casamet,casapool, &
                bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter,  &
                bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb, &
                bmpsoilocc,bmarea)
!  call spinanalytic(fcnpspin,dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
!                       casaflux,casamet,casabal,phen)

  nloop1= max(1,mloop-4)

  DO nloop=2,mloop

  OPEN(91,file=fcnpspin)
  read(91,*) myearspin
  DO nyear=1,myearspin
     read(91,901) ncfile
     call read_casa_dump(ncfile,casamet,casaflux,ktau,kend)

    DO idoy=1,mdyear
      ktauy=idoy*ktauday
       casamet%tairk(:)       = casamet%Tairkspin(:,idoy)
       casamet%tsoil(:,1)     = casamet%Tsoilspin_1(:,idoy)
       casamet%tsoil(:,2)     = casamet%Tsoilspin_2(:,idoy)
       casamet%tsoil(:,3)     = casamet%Tsoilspin_3(:,idoy)
       casamet%tsoil(:,4)     = casamet%Tsoilspin_4(:,idoy)
       casamet%tsoil(:,5)     = casamet%Tsoilspin_5(:,idoy)
       casamet%tsoil(:,6)     = casamet%Tsoilspin_6(:,idoy)
       casamet%moist(:,1)     = casamet%moistspin_1(:,idoy)
       casamet%moist(:,2)     = casamet%moistspin_2(:,idoy)
       casamet%moist(:,3)     = casamet%moistspin_3(:,idoy)
       casamet%moist(:,4)     = casamet%moistspin_4(:,idoy)
       casamet%moist(:,5)     = casamet%moistspin_5(:,idoy)
       casamet%moist(:,6)     = casamet%moistspin_6(:,idoy)
       casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
       casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
       casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
       casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
       call biogeochem(ktauy,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                      casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
                      cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                      nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
!    write(89,891) idoy,nyear,nloop,nptx,&
!                  casaflux%cgpp(nptx),casaflux%cnpp(nptx),casaflux%crmplant(nptx,:),casaflux%Crgplant(nptx),casaflux%Crsoil(nptx), &
!                  casaflux%clabloss(nptx),casapool%dClabiledt(nptx),casaflux%fracClabile(nptx),                                  &
!                  casapool%Cplant(nptx,:), casapool%clitter(nptx,:), casapool%csoil(nptx,:),casabal%cbalance(nptx)
    ENDDO   ! end of idoy

  ENDDO   ! end of nyear
  close(91)

891 format(4(i6,2x),100(f9.2,1x))

    if(nloop>=nloop1) then
       print *, 'kloop =', 2+nloop-nloop1, nloop,nloop1
       call totcnppools(2+nloop-nloop1,veg,casamet,casapool, &
            bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
            bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb, &
            bmpsoilocc,bmarea)
    endif

  ENDDO     ! end of nloop

  ! write the last five loop pool size by PFT type
  open(92,file='cnpspinlast5.txt', position='append')
  write(92,*) 'myearspin =', myearspin
  write(92,921)
921 format('PFT total area in 10**12 m2', f12.4)
  do nvt=1,mvtype
     write(92,*) bmarea(nvt)
  enddo

  do nvt=1,mvtype
  if(bmarea(nvt) >0.0) then
     do kloop=1,5
        write(92,922) nvt, bmcplant(kloop,nvt,:),bmclitter(kloop,nvt,:),bmcsoil(kloop,nvt,:)
     enddo
     if (icycle >1) then
        do kloop=1,5
           write(92,922) nvt, bmnplant(kloop,nvt,:),bmnlitter(kloop,nvt,:),bmnsoil(kloop,nvt,:), bmnsoilmin(kloop,nvt)
        enddo
     endif

     if(icycle >2) then
        do kloop=1,5
           write(92,922) nvt, bmpplant(kloop,nvt,:),bmplitter(kloop,nvt,:),bmpsoil(kloop,nvt,:),  &
                              bmpsoillab(kloop,nvt), bmpsoilsorb(kloop,nvt), bmpsoilocc(kloop,nvt)
        enddo
     endif
  endif
  enddo

901  format(A99)
922 format(i4,20(f10.4,2x))
    CLOSE(92)
151 FORMAT(i6,100(f12.5,2x))
END SUBROUTINE spincasacnp




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

  real(r_2), dimension(mp)             :: cuemet, cuestr,cuecwd
  real, parameter                      :: cnmic=10.0                ! microbial biomass C:N ratio

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

!  WHERE(casamet%iveg2/=icewater)
!    cuemet(:) = (cnmic)*(casapool%nlitter(:,mic)/(casapool%clitter(:,mic)+1.0e-3))**0.76
!    cuestr(:) = (cnmic)*(casapool%nlitter(:,str)/(casapool%clitter(:,str)+1.0e-3))**0.76
!    WHERE(casapool%clitter(:,cwd) >1.0e-3)
!      cuecwd(:) = (0.45/0.52)*(cnmic)*(casapool%nlitter(:,cwd)/(casapool%clitter(:,cwd)+1.0e-3))**0.76
!    ELSEWHERE
!      cuecwd(:) = 0.4
!    ENDWHERE
!  ENDWHERE

  cuemet(:) = 0.45
  cuestr(:) = 0.45
  cuecwd(:) = 0.4

!  write(517,*),'before analyticpool:casaflux%klitter(39:40,1:3),casa%ksoil(39:40,1:3),avgxkNlimiting(39:40),avgxklitter(39:40),casabiome%fracLigninplant(veg%iveg(39:40),1:3),casapool%clitter(39:40,1:3),casa%csoil%(39:40,1:3)',casaflux%klitter(39:40,1:3),casaflux%ksoil(39:40,1:3),avgxkNlimiting(39:40),avgxklitter(39:40),casabiome%fracLigninplant(veg%iveg(39:40),1:3),casapool%clitter(39:40,:),casapool%csoil(39:40,:)
  do npt=1,mp
  if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) THEN
    casaflux%fromLtoS(npt,mic,metb)   = cuemet(npt)
                                          ! metb -> mic
    casaflux%fromLtoS(npt,mic,str)   = (1.0-casabiome%fracLigninplant(veg%iveg(npt),leaf)) * cuestr(npt)
                                          ! str -> mic
    casaflux%fromLtoS(npt,slow,str)  = casabiome%fracLigninplant(veg%iveg(npt),leaf) *cuestr(npt)
                                          ! str -> slow
    casaflux%fromLtoS(npt,mic,cwd)   = (1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood)) *cuecwd(npt)
                                          ! CWD -> fmic
    casaflux%fromLtoS(npt,slow,cwd)  =  casabiome%fracLigninplant(veg%iveg(npt),wood) * cuecwd(npt)
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
!                                   * casapool%ratioNCsoil(npt,mic)/casaflux%ksoil(npt,mic)
!                                  * casapool%ratioNCsoilnew(npt,mic)/casaflux%ksoil(npt,mic)
                                   * avgratioNCsoilmic(npt)/casaflux%ksoil(npt,mic)
         casapool%nsoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
!                                   * casapool%ratioNCsoil(npt,slow)/casaflux%ksoil(npt,slow)
!                                  * casapool%ratioNCsoilnew(npt,slow)/casaflux%ksoil(npt,slow)
                                   * avgratioNCsoilslow(npt)/casaflux%ksoil(npt,slow)
         casapool%nsoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
!                                   * casapool%ratioNCsoil(npt,pass)/casaflux%ksoil(npt,pass)
!                                   * casapool%ratioNCsoilnew(npt,pass)/casaflux%ksoil(npt,pass)
                                   * avgratioNCsoilpass(npt)/casaflux%ksoil(npt,pass)
          casapool%Nsoilmin(npt)    = avgnsoilmin(npt)

        ENDIF

        IF (icycle<=2) THEN
            totpsoil(npt)          = psorder(casamet%isorder(npt)) *xpsoil50(casamet%isorder(npt))
            casapool%plitter(npt,:)= casapool%Nlitter(npt,:) / casapool%ratioNPlitter(npt,:)
            casapool%psoil(npt,:)  = casapool%Nsoil(npt,:)   / casapool%ratioNPsoil(npt,:)
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
                                   * (casapool%ratioNCsoil(npt,pass)/casapool%ratioNPsoil(npt,pass))/casaflux%ksoil(npt,pass)
          ! assign the mineral pools
          casapool%psoillab(npt)      = avgpsoillab(npt)
          casapool%psoilsorb(npt)     = avgPsoilsorb(npt)
          casapool%psoilocc(npt)      = avgPsoilocc(npt)
        ENDIF
  ENDIF
  ENDDO

!  write(517,*),'before analyticpool:casaflux%klitter(39:40,1:3),casa%ksoil(39:40,1:3),avgxkNlimiting(39:40),avgxklitter(39:40),casabiome%fracLigninplant(veg%iveg(39:40),1:3),casapool%clitter(39:40,1:3),casa%csoil%(39:40,1:3)',casaflux%klitter(39:40,1:3),casaflux%ksoil(39:40,1:3),avgxkNlimiting(39:40),avgxklitter(39:40),casabiome%fracLigninplant(veg%iveg(39:40),1:3),casapool%clitter(39:40,1:3),casapool%csoil(39:40,:)
  END SUBROUTINE analyticpool



  ! added by ypwang following Chris Lu 5/nov/2012
  subroutine ncdf_dump(casamet, n_call, kend, ncfile)
        use netcdf
        use cable_def_types_mod
        use casadimension, only : mdyear, mplant
        USE casavariable
        use casa_dump_module, only : def_dims, def_vars, def_var_atts, &
                                     put_var_nc, stderr_nc
        USE cable_io_vars_module, only : patch

        implicit none
        !var (type) to write
        TYPE (casa_met),            INTENT(INOUT) :: casamet

        integer , intent(in) :: &
           n_call, &         ! this timestep #
           kend              ! final timestep of run

        !number of instances. dummied here and so=1
        !integer :: inst =1

        !netcdf IDs/ names
        character(len=*),intent(in):: ncfile
        integer, parameter :: num_vars=7
        integer, parameter :: num_dims=4

        integer, dimension(ms)     :: &
           soil

        integer, dimension(mvtype) :: &
           tile

        integer, dimension(mplant) :: &
           plant

        integer, save :: ncid       ! netcdf file ID

        !vars
        character(len=*), dimension(num_vars), parameter :: &
              var_name =  (/  "latitude     ", &
                             "longitude    ", &
                              "casamet_tairk", &
                              "tsoil        ", &
                              "moist        ", &
                              "cgpp         ", &
                              "crmplant     " /)

        integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb

        !dims
        character(len=*), dimension(num_dims), parameter :: &
              dim_name =  (/ "tile", &
                             "soil", &
                             "plant",&
                             "time" /)

        integer, dimension(num_dims)  :: &
              dimID   ! (1) mp, (2) ms, (3) mplant (4) time

        integer, dimension(num_dims)  :: &
              !x,y generally lat/lon BUT for single site = 1,1
              dim_len
        !local only
       integer :: ncok      !ncdf return status
       integer :: i,j


        real(r_2), dimension(1:mp) :: &
           cgpp,   &
           tairk

        real(r_2), dimension(1:mp,1:ms) :: &
           tsoil, &
           moist

        real(r_2), dimension(1:mp,1:mplant) :: &
           crmplant

        ! END header
  !      tairk = 0
  !      cgpp  = 0
  !      tsoil = 0
  !      moist = 0
  !      crmplant = 0
  !      write(89,*) 'tsoil,gpp',casamet%Tsoilspin_1(1,:),casamet%cgppspin(1,:),mplant,ms
  !      write(89,*) 'creating dump file'
  !      n_call = 1
  !      kend   = mdyear


        print *, 'yp wang: calling ncdf_dump'
  !      print *, 'latitude= ', patch(:)%latitude
  !      print *, 'longitude= ', patch(:)%longitude
        print *, 'filename= ', ncfile
        print *, 'constants= ', mp,ms,mplant,mdyear,ncid

        dim_len(1) = mp
        dim_len(2) = ms
        dim_len(3) = mplant
        dim_len(4) = mdyear
           ! create netCDF dataset: enter define mode
        ncok = nf90_create(path = ncfile, cmode = nf90_clobber, ncid = ncid)

  !      print *, 'ncok =', ncok

  !      ncok = nf90_create(path = 'dump_casamet.nc', cmode = nf90_noclobber, ncid = ncid)
        if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile)
  !      if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', 'dump_casamet.nc')

  !      print *, 'here 1' ,ncid

        ! define dimensions: from name and length
  !      write(89,*) 'defining dims'
        call def_dims(num_dims, ncid, dimID, dim_len, dim_name )

  !      print *, 'here 2',varID,num_dims
        ! define variables: from name, type, dims
  !      write(89,*) 'defining vars'
        call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

  !      print *, 'here 3',varID,num_vars
        ! define variable attributes
  !      write(89,*) 'defining attribution'
        call def_var_atts(ncfile, ncid, varID )
  !      call def_var_atts('dump_casamet.nc', ncid, varID )

  !      print *, 'here 4', varID
        ncok = nf90_enddef(ncid)

  !      print *, 'here 5', var_name(1), size(patch(:)%latitude)
  !      write(89,*) 'writing latitude'
        call put_var_nc(ncid, var_name(1), patch(:)%latitude )

  !      print *, 'here 6',var_name(2) , size(patch(:)%longitude)
  !      write(89,*) 'writing longitude'
        call put_var_nc(ncid, var_name(2), patch(:)%longitude )

        write(*,901)  mdyear
  901   format(' yp wang at ncdf_dump', I6)
  !      write(*,*) casamet%cgppspin(10,:)

        do i=1,mdyear
           tairk(:)      = casamet%Tairkspin(:,i)
           tsoil(:,1)    = casamet%Tsoilspin_1(:,i)
  !      tairk           = casamet%Tsoilspin_1
           tsoil(:,2)    = casamet%Tsoilspin_2(:,i)
           tsoil(:,3)    = casamet%Tsoilspin_3(:,i)
           tsoil(:,4)    = casamet%Tsoilspin_4(:,i)
           tsoil(:,5)    = casamet%Tsoilspin_5(:,i)
           tsoil(:,6)    = casamet%Tsoilspin_6(:,i)
           moist(:,1)    = casamet%moistspin_1(:,i)
           moist(:,2)    = casamet%moistspin_2(:,i)
           moist(:,3)    = casamet%moistspin_3(:,i)
           moist(:,4)    = casamet%moistspin_4(:,i)
           moist(:,5)    = casamet%moistspin_5(:,i)
           moist(:,6)    = casamet%moistspin_6(:,i)
           cgpp (:)      = casamet%cgppspin   (:,i)
           crmplant(:,1) = casamet%crmplantspin_1(:,i)
           crmplant(:,2) = casamet%crmplantspin_2(:,i)
           crmplant(:,3) = casamet%crmplantspin_3(:,i)

  !         write(89,*) 'writing ',var_name(3)
           call put_var_nc(ncid, var_name(3), tairk, i,kend )
  !         write(89,*) 'writing '//trim(var_name(4))
           call put_var_nc(ncid, var_name(4), tsoil, i,kend ,ms)
  !         write(89,*) 'writing '//trim(var_name(5))
           call put_var_nc(ncid, var_name(5), moist, i,kend ,ms)
  !         write(89,*) 'writing '//trim(var_name(6))
           call put_var_nc(ncid, var_name(6), cgpp , i,kend )
  !         write(89,*) 'writing '//trim(var_name(7))
           call put_var_nc(ncid, var_name(7), crmplant, i,kend,mplant )
       end do

       !      if (n_call == kend ) &
       ncok = nf90_close(ncid)            ! close: save new netCDF dataset

     end subroutine ncdf_dump


End module casa_cable
