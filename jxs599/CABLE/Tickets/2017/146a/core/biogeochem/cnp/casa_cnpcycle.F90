SUBROUTINE casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet, LALLOC)
! update all pool sizes
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  INTEGER , INTENT(IN) :: LALLOC
  ! local variables
  REAL(r_2), DIMENSION(mp)   :: plabsorb,deltap
  INTEGER i,j,k,np,nland

!  PRINT *, 'mp here is ', mp
!  print*,'cplant',casapool%cplant(:,leaf)

  DO np=1,mp
  IF(casamet%iveg2(np) == icewater) THEN
    casamet%glai(np)   = 0.0
  ELSE

!if (np==2) write(912,91) casapool%cplant(np,:),  casapool%dcplantdt(np,:)  * deltpool
91 format (100(e12.4,2x))
    casapool%cplant(np,:)  = casapool%cplant(np,:)  &
                           + casapool%dcplantdt(np,:)  * deltpool
    casapool%clabile(np)   = casapool%clabile(np)   &
                           + casapool%dclabiledt(np)   * deltpool



    IF(casapool%cplant(np,leaf) > 0.0) THEN
      IF(icycle >1) casapool%Nplant(np,:) = casapool%Nplant(np,:) &
                                 +casapool%dNplantdt(np,:)*deltpool
      IF(icycle >2) casapool%Pplant(np,:) = casapool%Pplant(np,:) &
                                 +casapool%dPplantdt(np,:)*deltpool
    ENDIF
!    casamet%glai(np)   = MIN(0.0, casabiome%sla(veg%iveg(np))  &
!                                  * casapool%cplant(np,leaf))
    casamet%glai(np)   = MAX(casabiome%glaimin(veg%iveg(np)), &
                               casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf))
   ! vh !
    IF (LALLOC.ne.3) THEN
       casamet%glai(np)   = MIN(casabiome%glaimax(veg%iveg(np)), casamet%glai(np))
    ENDIF
    casapool%clitter(np,:) = casapool%clitter(np,:) &
                           + casapool%dClitterdt(np,:) * deltpool
    casapool%csoil(np,:)   = casapool%csoil(np,:)   &
                           + casapool%dCsoildt(np,:)   * deltpool

    IF(icycle >1) THEN
      casapool%Nlitter(np,:) = casapool%Nlitter(np,:) &
                             + casapool%dNlitterdt(np,:)* deltpool
      casapool%Nsoil(np,:)   = casapool%Nsoil(np,:)   &
                             + casapool%dNsoildt(np,:)  * deltpool
      ! vh ! put lower bound of 1.e-3 to prevent Nsoilmin from going negative
        ! Ticket #108
      casapool%Nsoilmin(np)  = max(casapool%Nsoilmin(np)  &
                             + casapool%dNsoilmindt(np) * deltpool,1.e-3)
    ENDIF

    IF(icycle >2) THEN
      casapool%Plitter(np,:) = casapool%Plitter(np,:) &
                             + casapool%dPlitterdt(np,:)* deltpool
      casapool%Psoil(np,:)   = casapool%Psoil(np,:)   &
                             + casapool%dPsoildt(np,:)  * deltpool
      casapool%Psoillab(np)  = casapool%Psoillab(np)  &
                             + casapool%dPsoillabdt(np) * deltpool
      casapool%Psoilsorb(np) = casaflux%Psorbmax(np)*casapool%Psoillab(np) &
                             /(casaflux%kmlabp(np)+casapool%Psoillab(np))
!      casapool%Psoilsorb(np) = casapool%Psoilsorb(np)  &
!                             + casapool%dPsoilsorbdt(np) * deltpool
      casapool%Psoilocc(np)   = casapool%Psoilocc(np)  &
                              + casapool%dPsoiloccdt(np)  * deltpool
    ENDIF

    DO i=1,mplant
      IF(casapool%cplant(np,i) < 0.0)  THEN
        WRITE(57,*)  'Cpool: np,ivt',np,casamet%lat(np),casamet%lon(np), &
             casamet%iveg2(np),casapool%cplant(np,:)
        call casa_poolzero(np,1,casapool)
!stop
        casapool%cplant(np,i) = max(0.0, casapool%cplant(np,i))
      ENDIF
    ENDDO
    IF(icycle >1) THEN
      DO i=1,mplant
        IF(casapool%nplant(np,i) < 0.0) THEN
          WRITE(57,*) 'Npool:', 'np,ivt,ipool',np,casamet%iveg2(np),casapool%nplant(np,:)
          call casa_poolzero(np,2,casapool)
          casapool%nplant(np,i) = max(0.0, casapool%nplant(np,i))
        ENDIF
      ENDDO
    ENDIF ! end of "icycle >1"

    DO j=1,mlitter
      IF(casapool%clitter(np,j) < 0.0)  THEN
        WRITE(57,*)  'Clitter: np,ivt2',np,casamet%iveg2(np),casapool%clitter(np,:)
        call casa_poolzero(np,3,casapool)
        casapool%clitter(np,j) = max(0.0, casapool%clitter(np,j))
      ENDIF
    ENDDO

    DO k=1,msoil
      IF(casapool%csoil(np,k) < 0.0)    THEN
        WRITE(57,*)  'Csoil: np,ivt2',np,casamet%iveg2(np),casapool%csoil(np,:)
        call casa_poolzero(np,5,casapool)
        casapool%csoil(np,k) = max(0.0, casapool%csoil(np,k))
      ENDIF
    ENDDO

!  check if any pool size, and terminate model run if any pool size is negative!!
    IF(icycle >1) THEN
      DO j=1,mlitter
        IF(casapool%nlitter(np,j) < 0.0)  THEN
          WRITE(57,*)  'Nlitter: np,ivt2',np,casamet%iveg2(np),casapool%Nlitter(np,:)
          call casa_poolzero(np,4,casapool)
          casapool%nlitter(np,j) = max(0.0, casapool%nlitter(np,j))
        ENDIF
      ENDDO
      DO k=1,msoil
        IF(casapool%nsoil(np,k) < 0.0) THEN
          WRITE(57,*)  'Nsoil: np,ivt2',np,casamet%iveg2(np),casapool%nsoil(np,:)
          call casa_poolzero(np,6,casapool)
          casapool%nsoil(np,k) = max(0.0, casapool%nsoil(np,k))
        ENDIF
      ENDDO
    ENDIF  !end of "icycle >1"
  ENDIF
  ENDDO !end of "np"



END SUBROUTINE casa_cnpcycle


