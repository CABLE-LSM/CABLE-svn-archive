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

!#define UM_BUILD YES
SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssnow,canopy,veg,soil, &
                     climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
                     pop, spinConv, spinup, ktauday, idoy,loy, dump_read,   &
                     dump_write, LALLOC)

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
   USE cable_phenology_module, ONLY: cable_phenology_clim

   IMPLICIT NONE

   INTEGER,      INTENT(IN) :: ktau ! integration step number
   INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
   INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

   INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
   INTEGER,      INTENT(IN)                  :: ktauday
   logical,      INTENT(IN) :: spinConv, spinup
   logical,      INTENT(IN) :: dump_read, dump_write
   INTEGER,      INTENT(IN)                  :: LALLOC

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
   TYPE(POP_TYPE),             INTENT(INOUT) :: POP
   TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

   ! local variables added ypwang 5/nov/2012
   real,      dimension(mp)  :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
   real,      dimension(mp)  :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
   real,      dimension(mp)  :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
   real(r_2), dimension(mp)  :: xnplimit,  xkNlimiting, xklitter, xksoil ,xkleaf,xkleafcold,xkleafdry

   INTEGER                                   :: it, nit
   REAL(dp)                               :: StemNPP(mp,2)
   CHARACTER                                 :: cyear*4
   CHARACTER                                 :: ncfile*99
 

  ! INTEGER, INTENT(IN) :: wlogn
   INTEGER , parameter :: wlogn=6

  
   IF ( .NOT. dump_read ) THEN  ! construct casa met and flux inputs from current CABLE run
      IF ( TRIM(cable_user%MetType) .EQ. 'cru' ) THEN
         casaflux%Nmindep = met%Ndep
      ENDIF

      IF(ktau == kstart) THEN
         casamet%tairk  = 0.0
         casamet%tsoil  = 0.0
         casamet%moist  = 0.0
         if(Ticket200) then
           casaflux%cgpp  = 0.0
           !add initializations (BP jul2010)
           casaflux%Crsoil   = 0.0
           casaflux%crgplant = 0.0
           casaflux%crmplant = 0.0
           casaflux%clabloss = 0.0
           !casaflux%crmplant(:,leaf) = 0.0
         End 
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
         casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + &
                                       canopy%frday*dels
      ENDIF

      IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
         casamet%tairk  =casamet%tairk/FLOAT(ktauday)
         casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
         casamet%moist=casamet%moist/FLOAT(ktauday)
  
         IF ( icycle .GT. 0 ) THEN
            IF (trim(cable_user%PHENOLOGY_SWITCH)=='climate') &
             ! get climate_dependent phenology
             call cable_phenology_clim(veg, climate, phen)
         
           if(Ticket200) then           
         
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
            CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
                casamet,casabal,phen,POP,climate, xnplimit,xkNlimiting,xklitter,xksoil, &
                xkleaf,xkleafcold,xkleafdry,&
                cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

!Ticket200 - YP's CALL - check and consolidate
!         CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
!                         casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,  &
!                         xksoil,xkleaf,xkleafcold,xkleafdry,&
!                         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,   &
!                         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,   &
!                         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

            IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP
               ! accumulate annual variables for use in POP
               IF(MOD(ktau/ktauday,LOY)==1 ) THEN
                  casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 
                  ! (assumes 70% of wood NPP is allocated above ground)
                  casabal%LAImax = casamet%glai
                  casabal%Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
                  casabal%Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
               ELSE
                  casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
                  casabal%LAImax = max(casamet%glai, casabal%LAImax)
                  casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
                  casabal%Crootmean = casabal%Crootmean + casapool%cplant(:,3)/real(LOY)/1000.
               ENDIF
            ELSE
               casaflux%stemnpp = 0.
            ENDIF ! CALL_POP

         ENDIF  ! icycle .gt. 0

      ENDIF  ! end of day

   ELSE ! dump_read: ! use casa met and flux inputs from dumpfile

      IF( MOD((ktau-kstart+1),ktauday) == 0 ) THEN  ! end of day

!Ticket200 - again consolidate call. VH also has POP loop etc
!       CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
!                       casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,  &
!                       xksoil,xkleaf,xkleafcold,xkleafdry,&
!                       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,   &
!                       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,   &
!                       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

         CALL biogeochem(ktau,dels,idoy,LALLOC,veg,soil,casabiome,casapool,casaflux, &
              casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf, &
              xkleafcold,xkleafdry,&
              cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
              nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
              pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

         IF (cable_user%CALL_POP) THEN ! accumulate input variables for POP

            ! accumulate annual variables for use in POP
            IF(MOD(ktau/ktauday,LOY)==1) THEN
               casaflux%stemnpp =  casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7 
               ! (assumes 70% of wood NPP is allocated above ground)
               casabal%LAImax = casamet%glai
               casabal%Cleafmean = casapool%cplant(:,1)/real(LOY)/1000.
               casabal%Crootmean = casapool%cplant(:,3)/real(LOY)/1000.
            ELSE
               casaflux%stemnpp = casaflux%stemnpp + casaflux%cnpp * casaflux%fracCalloc(:,2) * 0.7
               casabal%LAImax = max(casamet%glai, casabal%LAImax)
               casabal%Cleafmean = casabal%Cleafmean + casapool%cplant(:,1)/real(LOY)/1000.
               casabal%Crootmean = casabal%Crootmean +casapool%cplant(:,3)/real(LOY)/1000.
            ENDIF


         ENDIF ! CALL_POP

      ENDIF ! end of day

   ENDIF ! dump_read

END SUBROUTINE bgcdriver



