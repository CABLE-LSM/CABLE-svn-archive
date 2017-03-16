!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
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



