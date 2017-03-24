SUBROUTINE update_sum_casa(sum_casapool, sum_casaflux, casapool, casaflux, &
   sum_now, average_now, nsteps)

  IMPLICIT NONE
  TYPE (casa_pool)   , INTENT(INOUT) :: sum_casapool
  TYPE (casa_flux)   , INTENT(INOUT) :: sum_casaflux
  TYPE (casa_pool)   , INTENT(IN) :: casapool
  TYPE (casa_flux)   , INTENT(IN) :: casaflux
  LOGICAL, INTENT(IN) :: sum_now, average_now
  INTEGER, INTENT(IN) :: nsteps

        IF (sum_now) then

           sum_casapool%Clabile = sum_casapool%Clabile + casapool%Clabile
           sum_casapool%dClabiledt = sum_casapool%Clabile + casapool%Clabile
           sum_casapool%Cplant =sum_casapool%Cplant + casapool%Cplant
           sum_casapool%Nplant =  sum_casapool%Nplant + casapool%Nplant
           sum_casapool%Pplant =  sum_casapool%Pplant + casapool%Pplant
           sum_casapool%dCplantdt = sum_casapool%dCplantdt + casapool%dCplantdt
           sum_casapool%dNplantdt = sum_casapool%dNplantdt + casapool%dNplantdt
           sum_casapool%dPplantdt = sum_casapool%dPplantdt + casapool%dPplantdt
           sum_casapool%ratioNCplant = sum_casapool%ratioNCplant + casapool%ratioNCplant
           sum_casapool%ratioNPplant = sum_casapool%ratioNPplant + casapool%ratioNPplant
           sum_casapool%Nsoilmin =  sum_casapool%Nsoilmin + casapool%Nsoilmin
           sum_casapool%Psoillab = sum_casapool%Psoillab + casapool%Psoillab
           sum_casapool%Psoilsorb  = sum_casapool%Psoilsorb + casapool%Psoilsorb
           sum_casapool%Psoilocc = sum_casapool%Psoilocc + casapool%Psoilocc
           sum_casapool%dNsoilmindt =  sum_casapool%dNsoilmindt + casapool%dNsoilmindt
           sum_casapool%dPsoillabdt = sum_casapool%dPsoillabdt + casapool%dPsoillabdt
           sum_casapool%dPsoilsorbdt =sum_casapool%dPsoilsorbdt + casapool%dPsoilsorbdt
           sum_casapool%dPsoiloccdt = sum_casapool%dPsoiloccdt +casapool%dPsoiloccdt
           sum_casapool%Clitter = sum_casapool%Clitter + casapool%Clitter
           sum_casapool%Nlitter = sum_casapool%Nlitter  + casapool%Nlitter
           sum_casapool%Plitter =  sum_casapool%Plitter + casapool%Plitter
           sum_casapool%dClitterdt = sum_casapool%dClitterdt  + casapool%dClitterdt
           sum_casapool%dNlitterdt = sum_casapool%dNlitterdt  + casapool%dNlitterdt
           sum_casapool%dPlitterdt = sum_casapool%dPlitterdt + casapool%dPlitterdt
           sum_casapool%ratioNClitter = sum_casapool%ratioNClitter + casapool%ratioNClitter
           sum_casapool%ratioNPlitter = sum_casapool%ratioNPlitter + casapool%ratioNPlitter
           sum_casapool%Csoil =  sum_casapool%Csoil + casapool%Csoil
           sum_casapool%Nsoil =  sum_casapool%Nsoil + casapool%Nsoil
           sum_casapool%Psoil =  sum_casapool%Psoil + casapool%Psoil
           sum_casapool%dCsoildt = sum_casapool%dCsoildt + casapool%dCsoildt
           sum_casapool%dNsoildt = sum_casapool%dNsoildt + casapool%dNsoildt
           sum_casapool%dPsoildt = sum_casapool%dPsoildt + casapool%dPsoildt
           sum_casapool%ratioNCsoil = sum_casapool%ratioNCsoil + casapool%ratioNCsoil
           sum_casapool%ratioNPsoil = sum_casapool%ratioNPsoil + casapool%ratioNPsoil
           sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew + casapool%ratioNCsoilnew
           sum_casapool%ratioNCsoilmin =  sum_casapool%ratioNCsoilmin + casapool%ratioNCsoilmin
           sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax + casapool%ratioNCsoilmax
           sum_casapool%ratioPcsoil =  sum_casapool%ratioPcsoil  + casapool%ratioPcsoil
           sum_casapool%ratioPcplant =  sum_casapool%ratioPcplant + casapool%ratioPcplant
           sum_casapool%ratioPclitter =   sum_casapool%ratioPclitter + casapool%ratioPclitter

           sum_casaflux%Cgpp = sum_casaflux%Cgpp  + casaflux%Cgpp
           sum_casaflux%Cnpp = sum_casaflux%Cnpp + casaflux%Cnpp
           sum_casaflux%Crp = sum_casaflux%Crp  + casaflux%Crp
           sum_casaflux%Crgplant = sum_casaflux%Crgplant + casaflux%Crgplant
           sum_casaflux%Nminfix =  sum_casaflux%Nminfix + casaflux%Nminfix
           sum_casaflux%Nminuptake =  sum_casaflux%Nminuptake + casaflux%Nminuptake
           sum_casaflux%Plabuptake =  sum_casaflux%Plabuptake + casaflux%Plabuptake
           sum_casaflux%Clabloss =  sum_casaflux%Clabloss + casaflux%Clabloss
           sum_casaflux%fracClabile = sum_casaflux%fracClabile +  casaflux%fracClabile
           sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc + casaflux%fracCalloc*casapool%cplant
           sum_casaflux%fracNalloc = sum_casaflux%fracNalloc + casaflux%fracNalloc
           sum_casaflux%fracPalloc =  sum_casaflux%fracPalloc + casaflux%fracPalloc
           sum_casaflux%kplant =  sum_casaflux%kplant + casaflux%kplant*casapool%cplant
           sum_casaflux%Crmplant =   sum_casaflux%Crmplant + casaflux%Crmplant
           sum_casaflux%fromPtoL =  sum_casaflux%fromPtoL + casaflux%fromPtoL
           sum_casaflux%Cnep =  sum_casaflux%Cnep +  casaflux%Cnep
           sum_casaflux%Crsoil = sum_casaflux%Crsoil + casaflux%Crsoil
           sum_casaflux%Nmindep = sum_casaflux%Nmindep + casaflux%Nmindep
           sum_casaflux%Nminloss = sum_casaflux%Nminloss + casaflux%Nminloss
           sum_casaflux%Nminleach =  sum_casaflux%Nminleach + casaflux%Nminleach
           sum_casaflux%Nupland = sum_casaflux%Nupland + casaflux%Nupland
           sum_casaflux%Nlittermin =  sum_casaflux%Nlittermin +  casaflux%Nlittermin
           sum_casaflux%Nsmin =  sum_casaflux%Nsmin +  casaflux%Nsmin
           sum_casaflux%Nsimm =  sum_casaflux%Nsimm + casaflux%Nsimm
           sum_casaflux%Nsnet = sum_casaflux%Nsnet + casaflux%Nsnet
           sum_casaflux%fNminloss =  sum_casaflux%fNminloss + casaflux%fNminloss
           sum_casaflux%fNminleach =  sum_casaflux%fNminleach + casaflux%fNminleach
           sum_casaflux%Pdep = sum_casaflux%Pdep +  casaflux%Pdep
           sum_casaflux%Pwea = sum_casaflux%Pwea + casaflux%Pwea
           sum_casaflux%Pleach =  sum_casaflux%Pleach + casaflux%Pleach
           sum_casaflux%Ploss =  sum_casaflux%Ploss +  casaflux%Ploss
           sum_casaflux%Pupland =  sum_casaflux%Pupland  + casaflux%Pupland
           sum_casaflux%Plittermin =  sum_casaflux%Plittermin + casaflux%Plittermin
           sum_casaflux%Psmin = sum_casaflux%Psmin +  casaflux%Psmin
           sum_casaflux%Psimm = sum_casaflux%Psimm + casaflux%Psimm
           sum_casaflux%Psnet =  sum_casaflux%Psnet + casaflux%Psnet
           sum_casaflux%fPleach =  sum_casaflux%fPleach + casaflux%fPleach
           sum_casaflux%kplab =  sum_casaflux%kplab + casaflux%kplab
           sum_casaflux%kpsorb = sum_casaflux%kpsorb + casaflux%kpsorb
           sum_casaflux%kpocc =  sum_casaflux%kpocc + casaflux%kpocc
           sum_casaflux%kmlabP =  sum_casaflux%kmlabP + casaflux%kmlabP
           sum_casaflux%Psorbmax =  sum_casaflux%Psorbmax + casaflux%Psorbmax
           sum_casaflux%klitter =  sum_casaflux%klitter + casaflux%klitter
           sum_casaflux%ksoil =  sum_casaflux%ksoil + casaflux%ksoil
           sum_casaflux%fromLtoS =  sum_casaflux%fromLtoS + casaflux%fromLtoS
           sum_casaflux%fromStoS =  sum_casaflux%fromStoS + casaflux%fromStoS
           sum_casaflux%fromLtoCO2 =  sum_casaflux%fromLtoCO2 + casaflux%fromLtoCO2
           sum_casaflux%fromStoCO2 = sum_casaflux%fromStoCO2 + casaflux%fromStoCO2
           sum_casaflux%stemnpp =  sum_casaflux%stemnpp + casaflux%stemnpp
           sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood + casaflux%frac_sapwood
           sum_casaflux%sapwood_area = sum_casaflux%sapwood_area + casaflux%sapwood_area
           sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover + casaflux%Cplant_turnover
           sum_casaflux%Cplant_turnover_disturbance = sum_casaflux%Cplant_turnover_disturbance + &
                casaflux%Cplant_turnover_disturbance
           sum_casaflux%Cplant_turnover_crowding = sum_casaflux%Cplant_turnover_crowding + &
                casaflux%Cplant_turnover_crowding
           sum_casaflux%Cplant_turnover_resource_limitation =  sum_casaflux%Cplant_turnover_resource_limitation +  &
                casaflux%Cplant_turnover_resource_limitation


           sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter + casaflux%FluxCtolitter
           sum_casaflux%FluxNtolitter =  sum_casaflux%FluxNtolitter + casaflux%FluxNtolitter
           sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter + casaflux%FluxPtolitter

           sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil + casaflux%FluxCtosoil
           sum_casaflux%FluxNtosoil =  sum_casaflux%FluxNtosoil + casaflux%FluxNtosoil
           sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil  + casaflux%FluxPtosoil

           sum_casaflux%FluxCtoco2 =  sum_casaflux%FluxCtoco2 + casaflux%FluxCtoco2
        endif

           if (average_now) then
           sum_casapool%Clabile = sum_casapool%Clabile/real(nsteps)
           sum_casapool%dClabiledt = sum_casapool%Clabile/real(nsteps)
           where (sum_casapool%Cplant.gt.1.e-12) 
              sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc/sum_casapool%Cplant
           elsewhere
              sum_casaflux%fracCalloc = 0.0
           endwhere

           where (sum_casapool%Cplant.gt.1.e-12) 
              sum_casaflux%kplant =  sum_casaflux%kplant/sum_casapool%Cplant
           elsewhere
              sum_casaflux%kplant = 0.0
           endwhere
           sum_casapool%Cplant =sum_casapool%Cplant/real(nsteps)
           sum_casapool%Nplant =  sum_casapool%Nplant/real(nsteps)
           sum_casapool%Pplant =  sum_casapool%Pplant/real(nsteps)
           sum_casapool%dCplantdt = sum_casapool%dCplantdt/real(nsteps)
           sum_casapool%dNplantdt = sum_casapool%dNplantdt/real(nsteps)
           sum_casapool%dPplantdt = sum_casapool%dPplantdt/real(nsteps)
           sum_casapool%ratioNCplant = sum_casapool%ratioNCplant/real(nsteps)
           sum_casapool%ratioNPplant = sum_casapool%ratioNPplant/real(nsteps)
           sum_casapool%Nsoilmin =  sum_casapool%Nsoilmin/real(nsteps)
           sum_casapool%Psoillab = sum_casapool%Psoillab/real(nsteps)
           sum_casapool%Psoilsorb  = sum_casapool%Psoilsorb/real(nsteps)
           sum_casapool%Psoilocc = sum_casapool%Psoilocc/real(nsteps)
           sum_casapool%dNsoilmindt =  sum_casapool%dNsoilmindt/real(nsteps)
           sum_casapool%dPsoillabdt = sum_casapool%dPsoillabdt/real(nsteps)
           sum_casapool%dPsoilsorbdt =sum_casapool%dPsoilsorbdt/real(nsteps)
           sum_casapool%dPsoiloccdt = sum_casapool%dPsoiloccdt/real(nsteps)
           sum_casapool%Clitter = sum_casapool%Clitter/real(nsteps)
           sum_casapool%Nlitter = sum_casapool%Nlitter /real(nsteps)
           sum_casapool%Plitter =  sum_casapool%Plitter/real(nsteps)
           sum_casapool%dClitterdt = sum_casapool%dClitterdt /real(nsteps)
           sum_casapool%dNlitterdt = sum_casapool%dNlitterdt /real(nsteps)
           sum_casapool%dPlitterdt = sum_casapool%dPlitterdt/real(nsteps)
           sum_casapool%ratioNClitter = sum_casapool%ratioNClitter/real(nsteps)
           sum_casapool%ratioNPlitter = sum_casapool%ratioNPlitter/real(nsteps)
           sum_casapool%Csoil =  sum_casapool%Csoil/real(nsteps)
           sum_casapool%Nsoil =  sum_casapool%Nsoil/real(nsteps)
           sum_casapool%Psoil =  sum_casapool%Psoil/real(nsteps)
           sum_casapool%dCsoildt = sum_casapool%dCsoildt/real(nsteps)
           sum_casapool%dNsoildt = sum_casapool%dNsoildt/real(nsteps)
           sum_casapool%dPsoildt = sum_casapool%dPsoildt/real(nsteps)
           sum_casapool%ratioNCsoil = sum_casapool%ratioNCsoil/real(nsteps)
           sum_casapool%ratioNPsoil = sum_casapool%ratioNPsoil/real(nsteps)
           sum_casapool%ratioNCsoilnew = sum_casapool%ratioNCsoilnew/real(nsteps)
           sum_casapool%ratioNCsoilmin =  sum_casapool%ratioNCsoilmin/real(nsteps)
           sum_casapool%ratioNCsoilmax = sum_casapool%ratioNCsoilmax/real(nsteps)
           sum_casapool%ratioPcsoil =  sum_casapool%ratioPcsoil /real(nsteps)
           sum_casapool%ratioPcplant =  sum_casapool%ratioPcplant/real(nsteps)
           sum_casapool%ratioPclitter =   sum_casapool%ratioPclitter/real(nsteps)

           sum_casaflux%Cgpp = sum_casaflux%Cgpp /real(nsteps)
           sum_casaflux%Cnpp = sum_casaflux%Cnpp/real(nsteps)
           sum_casaflux%Crp = sum_casaflux%Crp /real(nsteps)
           sum_casaflux%Crgplant = sum_casaflux%Crgplant/real(nsteps)
           sum_casaflux%Nminfix =  sum_casaflux%Nminfix/real(nsteps)
           sum_casaflux%Nminuptake =  sum_casaflux%Nminuptake/real(nsteps)
           sum_casaflux%Plabuptake =  sum_casaflux%Plabuptake/real(nsteps)
           sum_casaflux%Clabloss =  sum_casaflux%Clabloss/real(nsteps)
           sum_casaflux%fracClabile = sum_casaflux%fracClabile/real(nsteps)
         !  sum_casaflux%fracCalloc =  sum_casaflux%fracCalloc/real(nsteps)
           sum_casaflux%fracNalloc = sum_casaflux%fracNalloc/real(nsteps)
           sum_casaflux%fracPalloc =  sum_casaflux%fracPalloc/real(nsteps)
          ! sum_casaflux%kplant =  sum_casaflux%kplant/real(nsteps)

          
           sum_casaflux%Crmplant =   sum_casaflux%Crmplant/real(nsteps)
           sum_casaflux%fromPtoL =  sum_casaflux%fromPtoL/real(nsteps)
           sum_casaflux%Cnep =  sum_casaflux%Cnep/real(nsteps)
           sum_casaflux%Crsoil = sum_casaflux%Crsoil/real(nsteps)
           sum_casaflux%Nmindep = sum_casaflux%Nmindep/real(nsteps)
           sum_casaflux%Nminloss = sum_casaflux%Nminloss/real(nsteps)
           sum_casaflux%Nminleach =  sum_casaflux%Nminleach/real(nsteps)
           sum_casaflux%Nupland = sum_casaflux%Nupland/real(nsteps)
           sum_casaflux%Nlittermin =  sum_casaflux%Nlittermin/real(nsteps)
           sum_casaflux%Nsmin =  sum_casaflux%Nsmin/real(nsteps)
           sum_casaflux%Nsimm =  sum_casaflux%Nsimm/real(nsteps)
           sum_casaflux%Nsnet = sum_casaflux%Nsnet/real(nsteps)
           sum_casaflux%fNminloss =  sum_casaflux%fNminloss/real(nsteps)
           sum_casaflux%fNminleach =  sum_casaflux%fNminleach/real(nsteps)
           sum_casaflux%Pdep = sum_casaflux%Pdep/real(nsteps)
           sum_casaflux%Pwea = sum_casaflux%Pwea/real(nsteps)
           sum_casaflux%Pleach =  sum_casaflux%Pleach/real(nsteps)
           sum_casaflux%Ploss =  sum_casaflux%Ploss/real(nsteps)
           sum_casaflux%Pupland =  sum_casaflux%Pupland /real(nsteps)
           sum_casaflux%Plittermin =  sum_casaflux%Plittermin/real(nsteps)
           sum_casaflux%Psmin = sum_casaflux%Psmin/real(nsteps)
           sum_casaflux%Psimm = sum_casaflux%Psimm/real(nsteps)
           sum_casaflux%Psnet =  sum_casaflux%Psnet/real(nsteps)
           sum_casaflux%fPleach =  sum_casaflux%fPleach/real(nsteps)
           sum_casaflux%kplab =  sum_casaflux%kplab/real(nsteps)
           sum_casaflux%kpsorb = sum_casaflux%kpsorb/real(nsteps)
           sum_casaflux%kpocc =  sum_casaflux%kpocc/real(nsteps)
           sum_casaflux%kmlabP =  sum_casaflux%kmlabP/real(nsteps)
           sum_casaflux%Psorbmax =  sum_casaflux%Psorbmax/real(nsteps)
           sum_casaflux%klitter =  sum_casaflux%klitter/real(nsteps)
           sum_casaflux%ksoil =  sum_casaflux%ksoil/real(nsteps)
           sum_casaflux%fromLtoS =  sum_casaflux%fromLtoS/real(nsteps)
           sum_casaflux%fromStoS =  sum_casaflux%fromStoS/real(nsteps)
           sum_casaflux%fromLtoCO2 =  sum_casaflux%fromLtoCO2/real(nsteps)
           sum_casaflux%fromStoCO2 = sum_casaflux%fromStoCO2/real(nsteps)
           sum_casaflux%stemnpp =  sum_casaflux%stemnpp/real(nsteps)
           sum_casaflux%frac_sapwood = sum_casaflux%frac_sapwood/real(nsteps)
           sum_casaflux%sapwood_area = sum_casaflux%sapwood_area/real(nsteps)
           sum_casaflux%Cplant_turnover = sum_casaflux%Cplant_turnover/real(nsteps)
           sum_casaflux%Cplant_turnover_disturbance = casaflux%Cplant_turnover_disturbance/real(nsteps)
           sum_casaflux% Cplant_turnover_crowding = sum_casaflux%Cplant_turnover_crowding/real(nsteps)
           sum_casaflux%Cplant_turnover_resource_limitation = &
                sum_casaflux%Cplant_turnover_resource_limitation/real(nsteps)
           sum_casaflux%FluxCtolitter = sum_casaflux%FluxCtolitter/real(nsteps)
           sum_casaflux%FluxNtolitter =  sum_casaflux%FluxNtolitter/real(nsteps)
           sum_casaflux%FluxPtolitter = sum_casaflux%FluxPtolitter/real(nsteps)

           sum_casaflux%FluxCtosoil = sum_casaflux%FluxCtosoil/real(nsteps)
           sum_casaflux%FluxNtosoil =  sum_casaflux%FluxNtosoil/real(nsteps)
           sum_casaflux%FluxPtosoil = sum_casaflux%FluxPtosoil /real(nsteps)

           sum_casaflux%FluxCtoco2 =  sum_casaflux%FluxCtoco2/real(nsteps)
        endif



END SUBROUTINE update_sum_casa


