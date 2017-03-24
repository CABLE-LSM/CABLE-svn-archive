
SUBROUTINE alloc_casavariable(casabiome,casapool,casaflux, &
      casamet,casabal,arraysize)

  use casaparm, ONLY : leaf
  IMPLICIT NONE
  TYPE (casa_biome)  , INTENT(INOUT) :: casabiome
  TYPE (casa_pool)   , INTENT(INOUT) :: casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: casaflux
  TYPE (casa_met)    , INTENT(INOUT) :: casamet
  TYPE (casa_balance), INTENT(INOUT) :: casabal
  INTEGER,             INTENT(IN) :: arraysize

  ALLOCATE(casabiome%ivt2(mvtype),                   &
           casabiome%xkleafcoldmax(mvtype),          &
           casabiome%xkleafcoldexp(mvtype),          &
           casabiome%xkleafdrymax(mvtype),           &
           casabiome%xkleafdryexp(mvtype),           &
           casabiome%glaimax(mvtype),                &
           casabiome%glaimin(mvtype),                &
           casabiome%sla(mvtype),                    &
           casabiome%ratiofrootleaf(mvtype),         &
           casabiome%kroot(mvtype),                  &
           casabiome%krootlen(mvtype),               &
           casabiome%rootdepth(mvtype),              &
           casabiome%kuptake(mvtype),                &
           casabiome%kminN(mvtype),                  &
           casabiome%KuplabP(mvtype),                &
           casabiome%kclabrate(mvtype),              &
           casabiome%xnpmax(mvtype),                 &
           casabiome%q10soil(mvtype),                &
           casabiome%xkoptlitter(mvtype),            &
           casabiome%xkoptsoil(mvtype),              &
           casabiome%xkplab(mso),                    &
           casabiome%xkpsorb(mso),                   &
           casabiome%xkpocc(mso),                    &
           casabiome%prodptase(mvtype),              &
           casabiome%costnpup(mvtype),               &
           casabiome%maxfinelitter(mvtype),          &
           casabiome%maxcwd(mvtype),                 &
           casabiome%nintercept(mvtype),             &
           casabiome%nslope(mvtype),                 &
           casabiome%plantrate(mvtype,mplant),       &
           casabiome%rmplant(mvtype,mplant),         &
           casabiome%fracnpptoP(mvtype,mplant),      &
           casabiome%fraclignin(mvtype,mplant),      &
           casabiome%fraclabile(mvtype,mplant),      &
           casabiome%ratioNCplantmin(mvtype,mplant), &
           casabiome%ratioNCplantmax(mvtype,mplant), &
           casabiome%ratioNPplantmin(mvtype,mplant), &
           casabiome%ratioNPplantmax(mvtype,mplant), &
           casabiome%fracLigninplant(mvtype,mplant), &
           casabiome%ftransNPtoL(mvtype,mplant),     &
           casabiome%ftransPPtoL(mvtype,mplant),     &
           casabiome%litterrate(mvtype,mlitter),     &
           casabiome%soilrate(mvtype,msoil),         &
         !  casabiome%ratioPcplantmax(mvtype,leaf),   &
         !  casabiome%ratioPcplantmin(mvtype,leaf)    &
         !! vh_js !!
           casabiome%ratioPcplantmax(mvtype,mplant),   &
           casabiome%ratioPcplantmin(mvtype,mplant)    &
          )

  ALLOCATE(casapool%Clabile(arraysize),               &
           casapool%dClabiledt(arraysize),            &
           casapool%Cplant(arraysize,mplant),         &
           casapool%Nplant(arraysize,mplant),         &
           casapool%Pplant(arraysize,mplant),         &
           casapool%dCplantdt(arraysize,mplant),      &
           casapool%dNplantdt(arraysize,mplant),      &
           casapool%dPplantdt(arraysize,mplant),      &
           casapool%ratioNCplant(arraysize,mplant),   &
           casapool%ratioNPplant(arraysize,mplant),   &
           casapool%Nsoilmin(arraysize),              &
           casapool%Psoillab(arraysize),              &
           casapool%Psoilsorb(arraysize),             &
           casapool%Psoilocc(arraysize),              &
           casapool%dNsoilmindt(arraysize),           &
           casapool%dPsoillabdt(arraysize),           &
           casapool%dPsoilsorbdt(arraysize),          &
           casapool%dPsoiloccdt(arraysize),           &
           casapool%Clitter(arraysize,mlitter),       &
           casapool%Nlitter(arraysize,mlitter),       &
           casapool%Plitter(arraysize,mlitter),       &
           casapool%dClitterdt(arraysize,mlitter),    &
           casapool%dNlitterdt(arraysize,mlitter),    &
           casapool%dPlitterdt(arraysize,mlitter),    &
           casapool%ratioNClitter(arraysize,mlitter), &
           casapool%ratioNPlitter(arraysize,mlitter), &
           casapool%Csoil(arraysize,msoil),           &
           casapool%Nsoil(arraysize,msoil),           &
           casapool%Psoil(arraysize,msoil),           &
           casapool%dCsoildt(arraysize,msoil),        &
           casapool%dNsoildt(arraysize,msoil),        &
           casapool%dPsoildt(arraysize,msoil),        &
           casapool%ratioNCsoil(arraysize,msoil),     &
           casapool%ratioNPsoil(arraysize,msoil),     &
           casapool%ratioNCsoilnew(arraysize,msoil),  &
           casapool%ratioNCsoilmin(arraysize,msoil),  &
           casapool%ratioNCsoilmax(arraysize,msoil),  &
           casapool%ratioPcsoil(arraysize,msoil),     &
           casapool%ratioPcplant(arraysize,mplant),   &
           casapool%ratioPclitter(arraysize,mlitter), &
           casapool%Ctot_0(arraysize),                &
           casapool%Ctot(arraysize)   )               
    
  ALLOCATE(casaflux%Cgpp(arraysize),                     &
           casaflux%Cnpp(arraysize),                     &
           casaflux%Crp(arraysize),                      &
           casaflux%Crgplant(arraysize),                 &
           casaflux%Nminfix(arraysize),                  &
           casaflux%Nminuptake(arraysize),               &
           casaflux%Plabuptake(arraysize),               &
           casaflux%Clabloss(arraysize),                 &
           casaflux%fracClabile(arraysize),              &
           casaflux%fracCalloc(arraysize,mplant),        &
           casaflux%fracNalloc(arraysize,mplant),        &
           casaflux%fracPalloc(arraysize,mplant),        &
           casaflux%kplant(arraysize,mplant),            &
           casaflux%Crmplant(arraysize,mplant),          &
           casaflux%fromPtoL(arraysize,mlitter,mplant),  &
           casaflux%Cnep(arraysize),                     &
           casaflux%Crsoil(arraysize),                   &
           casaflux%Nmindep(arraysize),                  &
           casaflux%Nminloss(arraysize),                 &
           casaflux%Nminleach(arraysize),                &
           casaflux%Nupland(arraysize),                  &
           casaflux%Nlittermin(arraysize),               &
           casaflux%Nsmin(arraysize),                    &
           casaflux%Nsimm(arraysize),                    &
           casaflux%Nsnet(arraysize),                    &
           casaflux%fNminloss(arraysize),                &
           casaflux%fNminleach(arraysize),               &
           casaflux%Pdep(arraysize),                     &
           casaflux%Pwea(arraysize),                     &
           casaflux%Pleach(arraysize),                   &
           casaflux%Ploss(arraysize),                    &
           casaflux%Pupland(arraysize),                  &
           casaflux%Plittermin(arraysize),               &
           casaflux%Psmin(arraysize),                    &
           casaflux%Psimm(arraysize),                    &
           casaflux%Psnet(arraysize),                    &
           casaflux%fPleach(arraysize),                  &
           casaflux%kplab(arraysize),                    &
           casaflux%kpsorb(arraysize),                   &
           casaflux%kpocc(arraysize),                    &
           casaflux%kmlabP(arraysize),                   &
           casaflux%Psorbmax(arraysize),                 &
           casaflux%klitter(arraysize,mlitter),          &
           casaflux%ksoil(arraysize,msoil),              &
           casaflux%fromLtoS(arraysize,msoil,mlitter),   &
           casaflux%fromStoS(arraysize,msoil,msoil),     &
           casaflux%fromLtoCO2(arraysize,mlitter),       &
           casaflux%fromStoCO2(arraysize,msoil),         &
           casaflux%stemnpp(arraysize),                  &
           casaflux%frac_sapwood(arraysize),             &
           casaflux%sapwood_area(arraysize), &
           casaflux%Cplant_turnover(arraysize,mplant) , &
           casaflux%Cplant_turnover_disturbance(arraysize) , &
           casaflux%Cplant_turnover_crowding(arraysize) , &
           casaflux%Cplant_turnover_resource_limitation(arraysize))

  ALLOCATE(casaflux%FluxCtolitter(arraysize,mlitter),    &
           casaflux%FluxNtolitter(arraysize,mlitter),    &
           casaflux%FluxPtolitter(arraysize,mlitter))

  ALLOCATE(casaflux%FluxCtosoil(arraysize,msoil),        &
           casaflux%FluxNtosoil(arraysize,msoil),        &
           casaflux%FluxPtosoil(arraysize,msoil))

  ALLOCATE(casaflux%FluxCtohwp(arraysize),    &
           casaflux%FluxNtohwp(arraysize),    &
           casaflux%FluxPtohwp(arraysize))

  ALLOCATE(casaflux%FluxCtoclear(arraysize),    &
           casaflux%FluxNtoclear(arraysize),    &
           casaflux%FluxPtoclear(arraysize))

  ALLOCATE(casaflux%CtransferLUC(arraysize))

  ALLOCATE(casaflux%FluxCtoco2(arraysize))

  ALLOCATE(casamet%glai(arraysize),                &
           casamet%lnonwood(arraysize),            &
           casamet%Tairk(arraysize),               &
           casamet%precip(arraysize),              &
           casamet%tsoilavg(arraysize),            &
           casamet%moistavg(arraysize),            &
           casamet%btran(arraysize),               &
           casamet%Tsoil(arraysize,ms),            &
           casamet%moist(arraysize,ms),            &
           casamet%iveg2(arraysize),               &
           casamet%ijgcm(arraysize),               &
           casamet%isorder(arraysize),             &
           casamet%lat(arraysize),                 &
           casamet%lon(arraysize),                 &
           casamet%areacell(arraysize),             &

           casamet%Tairkspin(arraysize,mdyear),     &
           casamet%cgppspin(arraysize,mdyear),      &
           casamet%crmplantspin_1(arraysize,mdyear),&
           casamet%crmplantspin_2(arraysize,mdyear),&
           casamet%crmplantspin_3(arraysize,mdyear),&
           casamet%Tsoilspin_1(arraysize,mdyear),   &
           casamet%Tsoilspin_2(arraysize,mdyear),   &
           casamet%Tsoilspin_3(arraysize,mdyear),   &
           casamet%Tsoilspin_4(arraysize,mdyear),   &
           casamet%Tsoilspin_5(arraysize,mdyear),   &
           casamet%Tsoilspin_6(arraysize,mdyear),   &
           casamet%moistspin_1(arraysize,mdyear),   &
           casamet%moistspin_2(arraysize,mdyear),   &
           casamet%moistspin_3(arraysize,mdyear),   &
           casamet%moistspin_4(arraysize,mdyear),   &
           casamet%moistspin_5(arraysize,mdyear),   &
           casamet%moistspin_6(arraysize,mdyear),  &
           casamet%mtempspin(arraysize,mdyear))     

  ALLOCATE(casabal%FCgppyear(arraysize),           &
           casabal%FCnppyear(arraysize),           &
           casabal%FCrpyear(arraysize),            &
           casabal%FCrmleafyear(arraysize),        &
           casabal%FCrmwoodyear(arraysize),        &
           casabal%FCrmrootyear(arraysize),        &
           casabal%FCrgrowyear(arraysize),         &
           casabal%FCrsyear(arraysize),            &
           casabal%FCneeyear(arraysize),           &
           casabal%FNdepyear(arraysize),           &
           casabal%FNfixyear(arraysize),           &
           casabal%FNsnetyear(arraysize),          &
           casabal%FNupyear(arraysize),            &
           casabal%FNleachyear(arraysize),         &
           casabal%FNlossyear(arraysize),          &
           casabal%FPweayear(arraysize),           &
           casabal%FPdustyear(arraysize),          &
           casabal%FPsnetyear(arraysize),          &
           casabal%FPupyear(arraysize),            &
           casabal%FPleachyear(arraysize),         &
           casabal%FPlossyear(arraysize),          &
           casabal%dCdtyear(arraysize),            & 
           casabal%LAImax(arraysize),              &  
           casabal%Cleafmean(arraysize),           & 
           casabal%Crootmean(arraysize)            ) 
  
     
  ALLOCATE(casabal%glaimon(arraysize,12),          &
           casabal%glaimonx(arraysize,12))

  ALLOCATE(casabal%cplantlast(arraysize,mplant),   &
           casabal%nplantlast(arraysize,mplant),   &
           casabal%pplantlast(arraysize,mplant))

  ALLOCATE(casabal%clitterlast(arraysize,mlitter), &
           casabal%nlitterlast(arraysize,mlitter), &
           casabal%plitterlast(arraysize,mlitter))

  ALLOCATE(casabal%csoillast(arraysize,msoil),     &
           casabal%nsoillast(arraysize,msoil),     &
           casabal%psoillast(arraysize,msoil))

  ALLOCATE(casabal%nsoilminlast(arraysize),        &
           casabal%psoillablast(arraysize),        &
           casabal%psoilsorblast(arraysize),       &
           casabal%psoilocclast(arraysize),        &
           casabal%cbalance(arraysize),            &
           casabal%nbalance(arraysize),            &
           casabal%pbalance(arraysize),            &
           casabal%sumcbal(arraysize),             &
           casabal%sumnbal(arraysize),             &
           casabal%sumpbal(arraysize),             &
           casabal%clabilelast(arraysize))       
END SUBROUTINE alloc_casavariable


