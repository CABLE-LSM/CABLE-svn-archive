
MODULE casavariable
  USE casadimension
  IMPLICIT NONE

  SAVE

  TYPE casa_biome
    INTEGER,   DIMENSION(:),POINTER :: ivt2
    REAL(r_2), DIMENSION(:),POINTER :: xkleafcoldmax,  &
                                       xkleafcoldexp,  &
                                       xkleafdrymax,   &
                                       xkleafdryexp,   &
                                       glaimax,        &
                                       glaimin,        &
                                       sla,            &
                                       ratiofrootleaf, &
                                       kroot,          &
                                       krootlen,       &
                                       rootdepth,      &
                                       kuptake,        &
                                       kminN,          &
                                       kuplabP,        &
                                       kclabrate,      &
                                       xnpmax,         &
                                       q10soil,        &
                                       xkoptlitter,    &
                                       xkoptsoil,      &
                                       xkplab,         &
                                       xkpsorb,        &
                                       xkpocc,         &
                                       prodptase,      &
                                       costnpup,       &
                                       maxfinelitter,  &
                                       maxcwd,         &
                                       nintercept,     &
                                       nslope

    REAL(r_2), DIMENSION(:,:),POINTER :: plantrate,     &
                                       rmplant,         &
                                       fracnpptoP,      &
                                       fraclignin,      &
                                       fraclabile,      &
                                       ratioNCplantmin, &
                                       ratioNCplantmax, &
                                       ratioNPplantmin, &
                                       ratioNPplantmax, &
                                       fracLigninplant, &
                                       ftransNPtoL,     &
                                       ftransPPtoL,     &
                                       litterrate,      &
                                       ratioPcplantmin, &
                                       ratioPcplantmax
    REAL(r_2), DIMENSION(:,:),POINTER :: soilrate
  END TYPE casa_biome

  TYPE casa_pool
    REAL(r_2), DIMENSION(:),POINTER :: Clabile,       &
                                       dClabiledt,    &
                                       Ctot ,         &          !! vh_js !!
                                       Ctot_0
    REAL(r_2), DIMENSION(:,:),POINTER :: Cplant,      &
                                       Nplant,        &
                                       Pplant,        &
                                       dCplantdt,     &
                                       dNplantdt,     &
                                       dPplantdt,     &
                                       ratioNCplant,  &
                                       ratioNPplant
    REAL(r_2), DIMENSION(:),POINTER :: Nsoilmin,      &
                                       Psoillab,      &
                                       Psoilsorb,     &
                                       Psoilocc,      &
                                       dNsoilmindt,   &
                                       dPsoillabdt,   &
                                       dPsoilsorbdt,  &
                                       dPsoiloccdt
    REAL(r_2), DIMENSION(:,:), POINTER :: Clitter,    &
                                       Nlitter,       &
                                       Plitter,       &
                                       dClitterdt,    &
                                       dNlitterdt,    &
                                       dPlitterdt,    &
                                       ratioNClitter, &
                                       ratioNPlitter
    REAL(r_2), DIMENSION(:,:),POINTER :: Csoil,       &
                                       Nsoil,         &
                                       Psoil,         &
                                       dCsoildt,      &
                                       dNsoildt,      &
                                       dPsoildt,      &
                                       ratioNCsoil,   &
                                       ratioNCsoilnew,&
                                       ratioNPsoil,   &
                                       ratioNCsoilmin,&
                                       ratioNCsoilmax,&
                                       ratioPcsoil,   &
                                       ratioPcplant,  &
                                       ratioPclitter
 END TYPE casa_pool

  TYPE casa_flux
    REAL(r_2), DIMENSION(:),POINTER :: Cgpp,          &
                                       Cnpp,          &
                                       Crp,           &
                                       Crgplant,      &
                                       Nminfix,       &
                                       Nminuptake,    &
                                       Plabuptake,    &
                                       Clabloss,      &
                                       fracClabile, &
!! vh_js !! the 3 variables below are needed for POP coupling to CASA
                                       stemnpp, &
                                       frac_sapwood, &
                                       sapwood_area
    REAL(r_2), DIMENSION(:,:),POINTER :: fracCalloc,  &
                                       fracNalloc,    &
                                       fracPalloc,    &
                                       Crmplant,      &
                                       kplant, &
!! vh_js !! additional diagnostic
                                       Cplant_turnover
    REAL(r_2), DIMENSION(:,:,:),POINTER :: fromPtoL
    REAL(r_2), DIMENSION(:),POINTER :: Cnep,        &
                                       Crsoil,      &
                                       Nmindep,     &
                                       Nminloss,    &
                                       Nminleach,   &
                                       Nupland,     &
                                       Nlittermin,  &
                                       Nsmin,       &
                                       Nsimm,       &
                                       Nsnet,       &
                                       fNminloss,   &
                                       fNminleach,  &
                                       Pdep,        &
                                       Pwea,        &
                                       Pleach,      &
                                       Ploss,       &
                                       Pupland,     &
                                       Plittermin,  &
                                       Psmin,       &
                                       Psimm,       &
                                       Psnet,       &
                                       fPleach,     &
                                       kplab,       &
                                       kpsorb,      &
                                       kpocc,       &
                                       kmlabp,      &
                                       Psorbmax,    &
!! additional diagnostics for partitioning biomass turnover
                                       Cplant_turnover_disturbance, &
                                       Cplant_turnover_crowding , &
                                       Cplant_turnover_resource_limitation

    REAL(r_2), DIMENSION(:,:),POINTER    :: klitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: ksoil
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromLtoS
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromStoS
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromLtoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromStoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtosoil
    REAL(r_2), DIMENSION(:),POINTER      :: FluxCtoCO2
    REAL(r_2), DIMENSION(:),POINTER    :: FluxCtohwp
    REAL(r_2), DIMENSION(:),POINTER    :: FluxNtohwp
    REAL(r_2), DIMENSION(:),POINTER    :: FluxPtohwp
    REAL(r_2), DIMENSION(:),POINTER    :: FluxCtoclear
    REAL(r_2), DIMENSION(:),POINTER    :: FluxNtoclear
    REAL(r_2), DIMENSION(:),POINTER    :: FluxPtoclear
    REAL(r_2), DIMENSION(:),POINTER    :: CtransferLUC

  
  END TYPE casa_flux

  TYPE casa_met
    REAL(r_2), DIMENSION(:),POINTER    :: glai,     &
                                          Tairk,    &
                                          precip,   &
                                          tsoilavg, &
                                          moistavg, &
                                          btran
    INTEGER, DIMENSION(:), POINTER     :: lnonwood
    REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil,    &
                                          moist
    INTEGER, DIMENSION(:), POINTER     :: iveg2,    &
                                          ijgcm,    &
                                          isorder
    REAL(r_2), DIMENSION(:), POINTER   :: lat,      &
                                          lon,      &
                                          areacell
    ! added yp wang 5/nov/2012
    REAL(r_2), DIMENSION(:,:), POINTER :: Tairkspin,&
                                          cgppspin,&
                                          crmplantspin_1,&
                                          crmplantspin_2,&
                                          crmplantspin_3,&
                                          Tsoilspin_1,&
                                          Tsoilspin_2,&
                                          Tsoilspin_3,&
                                          Tsoilspin_4,&
                                          Tsoilspin_5,&
                                          Tsoilspin_6,&
                                          moistspin_1,&
                                          moistspin_2,&
                                          moistspin_3,&
                                          moistspin_4,&
                                          moistspin_5,&
                                          moistspin_6, &
                                          mtempspin

  END TYPE casa_met

  TYPE casa_balance
    REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear,FCnppyear,                 &
            FCrmleafyear,FCrmwoodyear,FCrmrootyear,FCrgrowyear,               &
            FCrpyear, FCrsyear,FCneeyear,  dCdtyear ,                         &
            LAImax, Cleafmean, Crootmean,                          &
            FNdepyear,FNfixyear, FNsnetyear,FNupyear, FNleachyear,FNlossyear, &
            FPweayear,FPdustyear,FPsnetyear,FPupyear, FPleachyear,FPlossyear
 
    REAL(r_2), DIMENSION(:,:),POINTER :: glaimon,glaimonx
    REAL(r_2), DIMENSION(:,:),POINTER :: cplantlast,nplantlast,pplantlast
    REAL(r_2), DIMENSION(:,:),POINTER :: clitterlast,nlitterlast,plitterlast
    REAL(r_2), DIMENSION(:,:),POINTER :: csoillast,nsoillast,psoillast
    REAL(r_2), DIMENSION(:),  POINTER :: nsoilminlast,psoillablast,  &
                                         psoilsorblast,psoilocclast, &
                                         cbalance,nbalance,pbalance, &
                                         sumcbal,sumnbal,sumpbal
    REAL(r_2), DIMENSION(:),POINTER   :: clabilelast
 END TYPE casa_balance

! The following declarations are removed and have to be passed using
! parameter list for each subroutine (BP apr2010)
!  TYPE (casa_biome)              :: casabiome
!  TYPE (casa_pool)               :: casapool
!  TYPE (casa_flux)               :: casaflux
!  TYPE (casa_met)                :: casamet
!  TYPE (casa_balance)            :: casabal

! Added filename type for casaCNP (BP apr2010)
  TYPE casafiles_type
    CHARACTER(LEN=99) :: cnpbiome    ! file for biome-specific BGC parameters
    CHARACTER(LEN=99) :: cnppoint    ! file for point-specific BGC inputs
    CHARACTER(LEN=99) :: cnpepool    ! file for end-of-run pool sizes
    CHARACTER(LEN=99) :: cnpipool=''    ! file for inital pool sizes
    CHARACTER(LEN=99) :: cnpmetin      ! met file for spin up
    CHARACTER(LEN=99) :: cnpmetout     ! met file for spin up
    CHARACTER(LEN=99) :: ndep          ! N deposition input file
! added yp wang
    CHARACTER(LEN=99) :: cnpspin       ! input file for spin up
    CHARACTER(LEN=99) :: dump_cnpspin  ! name of dump file for spinning casa-cnp

    CHARACTER(LEN=99) :: phen        ! leaf phenology datafile
    CHARACTER(LEN=99) :: cnpflux     ! modelled mean yearly CNP fluxes
    LOGICAL           :: l_ndep
! added vh
    CHARACTER(LEN=99) :: c2cdumppath='' ! cable2casa dump for casa spinup
  END TYPE casafiles_type
  TYPE(casafiles_type) :: casafile

Contains


END MODULE casavariable


