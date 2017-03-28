SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,fzeroflux)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
  LOGICAL, optional :: fzeroflux
  LOGICAL :: zeroflux
  LOGICAL :: Ticket200 = .false.
  INTEGER n
  
  IF( present(fzeroflux) ) THEN
    zeroflux = fzeroflux
  else
    zeroflux = .false.
  endif
    
  IF(zeroflux) THEN
     casabal%FCgppyear    = 0.0
     casabal%FCrpyear     = 0.0   
     casabal%FCrmleafyear = 0.0
     casabal%FCrmwoodyear = 0.0
     casabal%FCrmrootyear = 0.0
     casabal%FCrgrowyear  = 0.0
     casabal%FCnppyear    = 0.0
     casabal%FCrsyear     = 0.0
     casabal%FCneeyear    = 0.0
     casabal%dCdtyear    = 0.0
    

     casabal%FNdepyear    = 0.0
     casabal%FNfixyear    = 0.0
     casabal%FNsnetyear   = 0.0
     casabal%FNupyear     = 0.0
     casabal%FNleachyear  = 0.0
     casabal%FNlossyear   = 0.0

     casabal%FPweayear   = 0.0
     casabal%FPdustyear  = 0.0
     casabal%FPsnetyear  = 0.0
     casabal%FPupyear    = 0.0
     casabal%FPleachyear = 0.0
     casabal%FPlossyear  = 0.0

     casaflux%FluxCtohwp = 0.0
     casaflux%FluxNtohwp = 0.0
     casaflux%FluxPtohwp = 0.0
     casaflux%FluxCtoclear = 0.0
     casaflux%FluxNtoclear = 0.0
     casaflux%FluxPtoclear = 0.0
     casaflux%CtransferLUC = 0.02
  ELSE

     casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
     casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
     casabal%FCrmleafyear(:)  = casabal%FCrmleafyear(:)  + casaflux%Crmplant(:,leaf)    * deltpool
     casabal%FCrmwoodyear(:)  = casabal%FCrmwoodyear(:)  + casaflux%Crmplant(:,wood)    * deltpool
     casabal%FCrmrootyear(:)  = casabal%FCrmrootyear(:)  + casaflux%Crmplant(:,froot)    * deltpool
     casabal%FCrgrowyear      = casabal%FCrgrowyear  + casaflux%Crgplant              * deltpool
     casabal%FCnppyear        = casabal%FCnppyear + (casaflux%Cnpp+casapool%dClabiledt)   * deltpool
     casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
     casabal%FCneeyear = casabal%FCneeyear &
          + (casaflux%Cnpp+casapool%dClabiledt-casaflux%Crsoil) * deltpool
     if( .NOT. Ticket200 ) &     
       casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool
   
     IF (icycle >1) THEN
        casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
        casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
        casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
        casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
        casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
        casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
     ENDIF

     IF (icycle >2) THEN
        casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
        casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
        casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
        casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
        casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool
        casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool
     ENDIF

  ENDIF

END SUBROUTINE casa_cnpflux

