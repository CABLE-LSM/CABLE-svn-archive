SUBROUTINE alloc_sum_casavariable(  sum_casapool, sum_casaflux &
     ,arraysize)

  use casaparm, ONLY : leaf
  IMPLICIT NONE
  TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux
  INTEGER,             INTENT(IN) :: arraysize


 ALLOCATE(sum_casapool%Clabile(arraysize),               &
           sum_casapool%dClabiledt(arraysize),            &
           sum_casapool%Cplant(arraysize,mplant),         &
           sum_casapool%Nplant(arraysize,mplant),         &
           sum_casapool%Pplant(arraysize,mplant),         &
           sum_casapool%dCplantdt(arraysize,mplant),      &
           sum_casapool%dNplantdt(arraysize,mplant),      &
           sum_casapool%dPplantdt(arraysize,mplant),      &
           sum_casapool%ratioNCplant(arraysize,mplant),   &
           sum_casapool%ratioNPplant(arraysize,mplant),   &
           sum_casapool%Nsoilmin(arraysize),              &
           sum_casapool%Psoillab(arraysize),              &
           sum_casapool%Psoilsorb(arraysize),             &
           sum_casapool%Psoilocc(arraysize),              &
           sum_casapool%dNsoilmindt(arraysize),           &
           sum_casapool%dPsoillabdt(arraysize),           &
           sum_casapool%dPsoilsorbdt(arraysize),          &
           sum_casapool%dPsoiloccdt(arraysize),           &
           sum_casapool%Clitter(arraysize,mlitter),       &
           sum_casapool%Nlitter(arraysize,mlitter),       &
           sum_casapool%Plitter(arraysize,mlitter),       &
           sum_casapool%dClitterdt(arraysize,mlitter),    &
           sum_casapool%dNlitterdt(arraysize,mlitter),    &
           sum_casapool%dPlitterdt(arraysize,mlitter),    &
           sum_casapool%ratioNClitter(arraysize,mlitter), &
           sum_casapool%ratioNPlitter(arraysize,mlitter), &
           sum_casapool%Csoil(arraysize,msoil),           &
           sum_casapool%Nsoil(arraysize,msoil),           &
           sum_casapool%Psoil(arraysize,msoil),           &
           sum_casapool%dCsoildt(arraysize,msoil),        &
           sum_casapool%dNsoildt(arraysize,msoil),        &
           sum_casapool%dPsoildt(arraysize,msoil),        &
           sum_casapool%ratioNCsoil(arraysize,msoil),     &
           sum_casapool%ratioNPsoil(arraysize,msoil),     &
           sum_casapool%ratioNCsoilnew(arraysize,msoil),  &
           sum_casapool%ratioNCsoilmin(arraysize,msoil),  &
           sum_casapool%ratioNCsoilmax(arraysize,msoil),  &
           sum_casapool%ratioPcsoil(arraysize,msoil),     &
           sum_casapool%ratioPcplant(arraysize,mplant),   &
           sum_casapool%ratioPclitter(arraysize,mlitter)  &
          )

 ALLOCATE(sum_casaflux%Cgpp(arraysize),                     &
           sum_casaflux%Cnpp(arraysize),                     &
           sum_casaflux%Crp(arraysize),                      &
           sum_casaflux%Crgplant(arraysize),                 &
           sum_casaflux%Nminfix(arraysize),                  &
           sum_casaflux%Nminuptake(arraysize),               &
           sum_casaflux%Plabuptake(arraysize),               &
           sum_casaflux%Clabloss(arraysize),                 &
           sum_casaflux%fracClabile(arraysize),              &
           sum_casaflux%fracCalloc(arraysize,mplant),        &
           sum_casaflux%fracNalloc(arraysize,mplant),        &
           sum_casaflux%fracPalloc(arraysize,mplant),        &
           sum_casaflux%kplant(arraysize,mplant),            &
           sum_casaflux%Crmplant(arraysize,mplant),          &
           sum_casaflux%fromPtoL(arraysize,mlitter,mplant),  &
           sum_casaflux%Cnep(arraysize),                     &
           sum_casaflux%Crsoil(arraysize),                   &
           sum_casaflux%Nmindep(arraysize),                  &
           sum_casaflux%Nminloss(arraysize),                 &
           sum_casaflux%Nminleach(arraysize),                &
           sum_casaflux%Nupland(arraysize),                  &
           sum_casaflux%Nlittermin(arraysize),               &
           sum_casaflux%Nsmin(arraysize),                    &
           sum_casaflux%Nsimm(arraysize),                    &
           sum_casaflux%Nsnet(arraysize),                    &
           sum_casaflux%fNminloss(arraysize),                &
           sum_casaflux%fNminleach(arraysize),               &
           sum_casaflux%Pdep(arraysize),                     &
           sum_casaflux%Pwea(arraysize),                     &
           sum_casaflux%Pleach(arraysize),                   &
           sum_casaflux%Ploss(arraysize),                    &
           sum_casaflux%Pupland(arraysize),                  &
           sum_casaflux%Plittermin(arraysize),               &
           sum_casaflux%Psmin(arraysize),                    &
           sum_casaflux%Psimm(arraysize),                    &
           sum_casaflux%Psnet(arraysize),                    &
           sum_casaflux%fPleach(arraysize),                  &
           sum_casaflux%kplab(arraysize),                    &
           sum_casaflux%kpsorb(arraysize),                   &
           sum_casaflux%kpocc(arraysize),                    &
           sum_casaflux%kmlabP(arraysize),                   &
           sum_casaflux%Psorbmax(arraysize),                 &
           sum_casaflux%klitter(arraysize,mlitter),          &
           sum_casaflux%ksoil(arraysize,msoil),              &
           sum_casaflux%fromLtoS(arraysize,msoil,mlitter),   &
           sum_casaflux%fromStoS(arraysize,msoil,msoil),     &
           sum_casaflux%fromLtoCO2(arraysize,mlitter),       &
           sum_casaflux%fromStoCO2(arraysize,msoil),         &
           sum_casaflux%stemnpp(arraysize),                  &
           sum_casaflux%frac_sapwood(arraysize),             &
           sum_casaflux%sapwood_area(arraysize), &
           sum_casaflux%Cplant_turnover(arraysize,mplant) , &
           sum_casaflux%Cplant_turnover_disturbance(arraysize) , &
           sum_casaflux%Cplant_turnover_crowding(arraysize) , &
           sum_casaflux%Cplant_turnover_resource_limitation(arraysize))

  ALLOCATE(sum_casaflux%FluxCtolitter(arraysize,mlitter),    &
           sum_casaflux%FluxNtolitter(arraysize,mlitter),    &
           sum_casaflux%FluxPtolitter(arraysize,mlitter))

  ALLOCATE(sum_casaflux%FluxCtosoil(arraysize,msoil),        &
           sum_casaflux%FluxNtosoil(arraysize,msoil),        &
           sum_casaflux%FluxPtosoil(arraysize,msoil))

  ALLOCATE(sum_casaflux%FluxCtoco2(arraysize))

END SUBROUTINE alloc_sum_casavariable


