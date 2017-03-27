SUBROUTINE casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)
! update all pool sizes
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  REAL(r_2), DIMENSION(mp)   :: plabsorb,deltap
  INTEGER i,j,k,np,nland

  DO np=1,mp

  IF(casamet%iveg2(np) == icewater) THEN
    casamet%glai(np)   = 0.0
  ELSE  
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
    casamet%glai(np)   = casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf)

!! changes made by ypw on 09-08-2016 
      if(casamet%glai(np) < casabiome%glaimin(veg%iveg(np)) &
        .or. casamet%glai(np) > casabiome%glaimax(veg%iveg(np) ) ) then

     casamet%glai(np)         = MIN(casabiome%glaimax(veg%iveg(np)),MAX(casabiome%glaimin(veg%iveg(np)), &
                                casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf)))
     casapool%cplant(np,leaf) = casamet%glai(np)/ casabiome%sla(veg%iveg(np))
     casapool%nplant(np,leaf) = casapool%cplant(np,leaf) * casabiome%ratioNCplantmin(veg%iveg(np),leaf) 
     casapool%pplant(np,leaf) = casapool%nplant(np,leaf) / casabiome%ratioNPplantmax(veg%iveg(np),leaf)
    endif
    casapool%clitter(np,:) = casapool%clitter(np,:) &
                           + casapool%dClitterdt(np,:) * deltpool 
    casapool%csoil(np,:)   = casapool%csoil(np,:)   &
                           + casapool%dCsoildt(np,:)   * deltpool

    IF(icycle >1) THEN
      casapool%Nlitter(np,:) = casapool%Nlitter(np,:) &
                             + casapool%dNlitterdt(np,:)* deltpool
      casapool%Nsoil(np,:)   = casapool%Nsoil(np,:)   &
                             + casapool%dNsoildt(np,:)  * deltpool
      casapool%Nsoilmin(np)  = casapool%Nsoilmin(np)  &
                             + casapool%dNsoilmindt(np) * deltpool
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
        if(veg%iveg(np)/=14) then
           WRITE(57,*)  'Cpool: np,ivt',np,veg%iveg(np),casapool%cplant(np,:)
           call casa_poolzero(np,1,casapool)
        endif
        casapool%cplant(np,i) = max(0.0, casapool%cplant(np,i))
      ENDIF
    ENDDO
    IF(icycle >1) THEN
      DO i=1,mplant
        IF(casapool%nplant(np,i) < 0.0) THEN
          if(veg%iveg(np)/=14) then
             WRITE(57,*) 'Npool:', 'np,ivt,ipool',np,veg%iveg(np),casapool%nplant(np,:)
             call casa_poolzero(np,2,casapool)
          endif
          casapool%nplant(np,i) = max(0.0, casapool%nplant(np,i))
        ENDIF
      ENDDO
    ENDIF ! end of "icycle >1"

    DO j=1,mlitter
      IF(casapool%clitter(np,j) < 0.0)  THEN
        if(veg%iveg(np)/=14) then
           WRITE(57,*)  'Clitter: np,ivt',np,veg%iveg(np),casapool%clitter(np,:)
           call casa_poolzero(np,3,casapool) 
        endif
        casapool%clitter(np,j) = max(0.0, casapool%clitter(np,j))
      ENDIF
    ENDDO

    DO k=1,msoil
      IF(casapool%csoil(np,k) < 0.0)    THEN
        if(veg%iveg(np)/=14) then
           WRITE(57,*)  'Csoil: np,ivt',np,veg%iveg(np),casapool%csoil(np,:)
           call casa_poolzero(np,5,casapool) 
        endif
        casapool%csoil(np,k) = max(0.0, casapool%csoil(np,k))
      ENDIF
    ENDDO

!  check if any pool size, and terminate model run if any pool size is negative!!
    IF(icycle >1) THEN
      DO j=1,mlitter
        IF(casapool%nlitter(np,j) < 0.0)  THEN
        if(veg%iveg(np)/=14)  then
            WRITE(57,*)  'Nlitter: np,ivt',np,veg%iveg(np),casapool%Nlitter(np,:)
            call casa_poolzero(np,4,casapool) 
        endif
          casapool%nlitter(np,j) = max(0.0, casapool%nlitter(np,j))
        ENDIF
      ENDDO
      DO k=1,msoil
        IF(casapool%nsoil(np,k) < 0.0) THEN
        if(veg%iveg(np)/=14)   then
            WRITE(57,*)  'Nsoil: np,ivt',np,veg%iveg(np),casapool%nsoil(np,:)
            call casa_poolzero(np,6,casapool) 
        endif
          casapool%nsoil(np,k) = max(0.0, casapool%nsoil(np,k))
        ENDIF
      ENDDO
    ENDIF  !end of "icycle >1"
  ENDIF
  ENDDO !end of "np"

END SUBROUTINE casa_cnpcycle


