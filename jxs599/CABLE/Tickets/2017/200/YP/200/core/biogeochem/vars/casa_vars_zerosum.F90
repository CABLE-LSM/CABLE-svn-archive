SUBROUTINE zero_sum_casa(sum_casapool, sum_casaflux)

  IMPLICIT NONE
  TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux

           sum_casapool%Clabile = 0
           sum_casapool%dClabiledt = 0
           sum_casapool%Cplant = 0
           sum_casapool%Nplant = 0
           sum_casapool%Pplant = 0
           sum_casapool%dCplantdt = 0
           sum_casapool%dNplantdt = 0
           sum_casapool%dPplantdt = 0
           sum_casapool%ratioNCplant = 0
           sum_casapool%ratioNPplant = 0
           sum_casapool%Nsoilmin = 0
           sum_casapool%Psoillab = 0
           sum_casapool%Psoilsorb  = 0
           sum_casapool%Psoilocc = 0
           sum_casapool%dNsoilmindt = 0
           sum_casapool%dPsoillabdt = 0
           sum_casapool%dPsoilsorbdt = 0
           sum_casapool%dPsoiloccdt = 0
           sum_casapool%Clitter = 0
           sum_casapool%Nlitter = 0
           sum_casapool%Plitter = 0
           sum_casapool%dClitterdt = 0
           sum_casapool%dNlitterdt = 0
           sum_casapool%dPlitterdt = 0
           sum_casapool%ratioNClitter = 0
           sum_casapool%ratioNPlitter = 0
           sum_casapool%Csoil = 0
           sum_casapool%Nsoil = 0
           sum_casapool%Psoil = 0
           sum_casapool%dCsoildt = 0
           sum_casapool%dNsoildt = 0
           sum_casapool%dPsoildt = 0
           sum_casapool%ratioNCsoil = 0
           sum_casapool%ratioNPsoil = 0
           sum_casapool%ratioNCsoilnew = 0
           sum_casapool%ratioNCsoilmin = 0
           sum_casapool%ratioNCsoilmax = 0
           sum_casapool%ratioPcsoil = 0
           sum_casapool%ratioPcplant = 0
           sum_casapool%ratioPclitter = 0

           sum_casaflux%Cgpp = 0
           sum_casaflux%Cnpp = 0
           sum_casaflux%Crp = 0
           sum_casaflux%Crgplant = 0
           sum_casaflux%Nminfix = 0
           sum_casaflux%Nminuptake = 0
           sum_casaflux%Plabuptake = 0
           sum_casaflux%Clabloss = 0
           sum_casaflux%fracClabile = 0
           sum_casaflux%fracCalloc = 0
           sum_casaflux%fracNalloc = 0
           sum_casaflux%fracPalloc = 0
           sum_casaflux%kplant = 0
           sum_casaflux%Crmplant = 0
           sum_casaflux%fromPtoL = 0
           sum_casaflux%Cnep = 0
           sum_casaflux%Crsoil = 0
           sum_casaflux%Nmindep = 0
           sum_casaflux%Nminloss = 0
           sum_casaflux%Nminleach = 0
           sum_casaflux%Nupland = 0
           sum_casaflux%Nlittermin = 0
           sum_casaflux%Nsmin = 0
           sum_casaflux%Nsimm = 0
           sum_casaflux%Nsnet = 0
           sum_casaflux%fNminloss = 0
           sum_casaflux%fNminleach = 0
           sum_casaflux%Pdep = 0
           sum_casaflux%Pwea = 0
           sum_casaflux%Pleach = 0
           sum_casaflux%Ploss = 0
           sum_casaflux%Pupland = 0
           sum_casaflux%Plittermin = 0
           sum_casaflux%Psmin = 0
           sum_casaflux%Psimm = 0
           sum_casaflux%Psnet = 0
           sum_casaflux%fPleach = 0
           sum_casaflux%kplab = 0
           sum_casaflux%kpsorb = 0
           sum_casaflux%kpocc = 0
           sum_casaflux%kmlabP = 0
           sum_casaflux%Psorbmax = 0
           sum_casaflux%klitter = 0
           sum_casaflux%ksoil = 0
           sum_casaflux%fromLtoS = 0
           sum_casaflux%fromStoS = 0
           sum_casaflux%fromLtoCO2 = 0
           sum_casaflux%fromStoCO2 = 0
           sum_casaflux%stemnpp = 0
           sum_casaflux%frac_sapwood = 0
           sum_casaflux%sapwood_area = 0
           sum_casaflux%Cplant_turnover = 0
           sum_casaflux%Cplant_turnover_disturbance = 0
           sum_casaflux% Cplant_turnover_crowding = 0
           sum_casaflux%Cplant_turnover_resource_limitation = 0

           sum_casaflux%FluxCtolitter = 0
           sum_casaflux%FluxNtolitter = 0
           sum_casaflux%FluxPtolitter = 0

           sum_casaflux%FluxCtosoil = 0
           sum_casaflux%FluxNtosoil = 0
           sum_casaflux%FluxPtosoil = 0

           sum_casaflux%FluxCtoco2 = 0




END SUBROUTINE zero_sum_casa


