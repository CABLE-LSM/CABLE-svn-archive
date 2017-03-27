SUBROUTINE casa_cnpbal(casapool,casaflux,casabal)

  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_balance),          INTENT(INOUT) :: casabal

  ! local variables
  INTEGER :: npt
  REAL(r_2), DIMENSION(mp) :: cbalplant,  nbalplant,  pbalplant
  REAL(r_2), DIMENSION(mp) :: cbalsoil,   nbalsoil,   pbalsoil
  REAL(r_2), DIMENSION(mp) :: cbalplantx, nbalplantx, pbalplantx
  logical :: Ticket200 = .false.

  cbalplant(:) = 0.0
  cbalsoil(:)  = 0.0
  nbalplant(:) = 0.0
  nbalsoil(:)  = 0.0
  pbalplant(:) = 0.0
  pbalsoil(:)  = 0.0

  casabal%cbalance(:)  = 0.0
  casabal%nbalance(:)  = 0.0
  casabal%pbalance(:)  = 0.0

!C balance
   Cbalplant(:)  = sum(casabal%cplantlast,2) -sum(casapool%cplant,2)            &
                 + casabal%Clabilelast(:)-casapool%clabile(:)                   &
                 +(casaflux%Cnpp(:) - SUM((casaflux%kplant*casabal%cplantlast),2))*deltpool &
                 + casapool%dClabiledt(:)* deltpool
   Cbalsoil(:)   = sum(casabal%clitterlast,2) - sum(casapool%clitter,2)         &
                 + sum(casabal%csoillast,2)   - sum(casapool%csoil,2)           &
                 +(SUM((casaflux%kplant*casabal%cplantlast),2)-casaflux%Crsoil(:))*deltpool

   casabal%cbalance(:) = Cbalplant(:) + Cbalsoil(:)

if(.NOT. Ticket200) then
 do npt=1,mp
    IF(abs(casabal%cbalance(npt))>1e-10) THEN
      write(*,*) 'cbalance',  npt, Cbalplant(npt), Cbalsoil(npt)
      write(*,*) 'cplant', casapool%cplant(npt,:)
      write(*,*) 'gpp, npp',casaflux%Cgpp(npt) , &
           casaflux%Cnpp(npt)
      write(*,*) 'dcplandt',  casapool%dcplantdt(npt,:), sum(casapool%dcplantdt(npt,:))
      write(*,*) 'rmplant, rgplant',  casaflux%crmplant(npt,:) , casaflux%crgplant(npt)
      write(*,*), 'dclabile',  casapool%dClabiledt(npt)* deltpool
       
     ENDIF
 ENDDO

   casapool%ctot_0 = sum(casabal%cplantlast,2)+sum(casabal%clitterlast,2) &
        + sum(casabal%csoillast,2)+ casabal%clabilelast
   casapool%ctot = sum(casapool%cplant,2)+sum(casapool%clitter,2) &
        + sum(casapool%csoil,2)+ casapool%clabile
endif   
   
   casabal%cplantlast  = casapool%cplant
   casabal%clabilelast = casapool%clabile
   casabal%clitterlast = casapool%clitter
   casabal%csoillast   = casapool%csoil
   casabal%sumcbal     = casabal%sumcbal + casabal%cbalance

   IF(icycle >1) THEN
      Nbalplant(:) = sum(casabal%nplantlast,2) -sum(casapool%nplant,2)                  &
                    +casaflux%Nminuptake(:) *deltpool
      Nbalsoil(:)  = -sum(casapool%nlitter,2)-sum(casapool%nsoil,2)                     &
                     -casapool%nsoilmin(:)+ casabal%nsoilminlast(:)                     &
                     + sum(casabal%nlitterlast,2)    + sum(casabal%nsoillast,2)         &
                     +(casaflux%Nmindep(:) + casaflux%Nminfix(:)- casaflux%Nminloss(:)                        &
                       -casaflux%Nminleach(:)-casaflux%Nupland(:)) * deltpool

      casabal%nbalance(:) = Nbalplant(:) + Nbalsoil(:)

      casabal%nplantlast  = casapool%nplant
      casabal%nlitterlast = casapool%nlitter
      casabal%nsoillast   = casapool%nsoil
      casabal%nsoilminlast= casapool%nsoilmin
      casabal%sumnbal     = casabal%sumnbal + casabal%nbalance

   ENDIF

   IF(icycle >2) THEN
      Pbalplant(:) = sum(casabal%Pplantlast,2) -sum(casapool%Pplant,2)                      &
                   + casaflux%Plabuptake(:) *deltpool
      Pbalsoil(:)  = -sum(casapool%Plitter,2)        - sum(casapool%Psoil,2)                      &
                     + sum(casabal%Plitterlast,2)    + sum(casabal%Psoillast,2)                   &
                   -casapool%psoillab(:)-casapool%psoilsorb(:)-casapool%psoilocc(:)               &
                   + casabal%psoillablast(:) + casabal%psoilsorblast(:) + casabal%psoilocclast(:) &
                   +(casaflux%Pdep(:) + casaflux%Pwea(:)                                          &
                     -casaflux%Pleach(:)-casaflux%Pupland(:)                                      &
                     -casaflux%Ploss(:)) * deltpool

      casabal%pbalance(:) = pbalplant(:) + pbalsoil(:)

      casabal%pplantlast   = casapool%pplant
      casabal%plitterlast  = casapool%plitter
      casabal%psoillast    = casapool%psoil
      casabal%psoillablast = casapool%psoillab
      casabal%psoilsorblast= casapool%psoilsorb
      casabal%psoilocclast = casapool%psoilocc
      casabal%sumpbal  = casabal%sumpbal + casabal%pbalance
   ENDIF

END SUBROUTINE casa_cnpbal