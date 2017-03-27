SUBROUTINE casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
! calculate changes in litter and soil pools

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome

  ! local variables
  REAL(r_2), DIMENSION(mp)    :: xdplabsorb, fluxptase
  INTEGER i,j,jj,k,kk,kkk,n,iv,npt,nL,nS,nSS,nland
  logical :: Ticket200 = .false.
  
  casaflux%fluxCtoCO2    = 0.0
  casaflux%fluxCtosoil   = 0.0
  casaflux%fluxNtosoil   = 0.0
  casaflux%fluxPtosoil   = 0.0
  casaflux%Crsoil        = 0.0  ! initialization added by BP jul2010

  casapool%dClitterdt = 0.0
  casapool%dCsoildt   = 0.0

  casapool%dNlitterdt = 0.0
  casapool%dNsoildt   = 0.0
  casapool%dNsoilmindt= 0.0
  casaflux%Nsmin = 0.0
  casaflux%Nsimm = 0.0
  casaflux%Nsnet = 0.0
  casaflux%Nminloss = 0.0
  casaflux%Nminleach = 0.0
  casaflux%Nlittermin=0.0

  casapool%dPlitterdt   = 0.0
  casapool%dPsoildt     = 0.0
  casapool%dPsoillabdt  = 0.0
  casapool%dPsoilsorbdt = 0.0
  casapool%dPsoiloccdt  = 0.0
  casaflux%Psmin = 0.0
  casaflux%Psimm = 0.0
  casaflux%Psnet = 0.0
  casaflux%Pleach = 0.0
  casaflux%Ploss  = 0.0
  casaflux%Plittermin=0.0
  fluxptase = 0.0

DO nland=1,mp
IF(casamet%iveg2(nland)/=icewater) THEN
   !Ticket200
   IF(icycle > 1) THEN
      !vh! set klitter to zero where Nlitter will go -ve 
      !(occurs occasionally for metabolic litter pool) Ticket#108
      where (casaflux%klitter(nland,:) * max(0.0,casapool%Nlitter(nland,:)).gt. &
           casapool%Nlitter(nland,:)+casaflux%fluxNtolitter(nland,:)) casaflux%klitter(nland,:) = 0.0
   endif

   DO nL=1,mlitter
      casaflux%fluxCtoCO2(nland) = casaflux%fluxCtoCO2(nland)  &
                        + casaflux%fromLtoCO2(nland,nL)  &
                        * casaflux%klitter(nland,nL) &
                        * casapool%clitter(nland,nL)
   ENDDO

   DO nS=1,msoil
      DO nL=1,mlitter
         casaflux%fluxCtosoil(nland,nS) = casaflux%fluxCtosoil(nland,nS) &
                               + casaflux%fromLtoS(nland,nS,nL) &
                               * casaflux%klitter(nland,nL) &
                               * casapool%clitter(nland,nL)
      ENDDO
      DO nSS=1,msoil
         IF(nSS/=nS) THEN
            casaflux%fluxCtosoil(nland,nS) = casaflux%fluxCtosoil(nland,nS) &
                                  + casaflux%fromStoS(nland,nS,nSS) &
                                  * casaflux%ksoil(nland,nSS) &
                                  * casapool%csoil(nland,nSS)
         ENDIF
      ENDDO
      casaflux%fluxCtoCO2(nland) = casaflux%fluxCtoCO2(nland)  &
                        + casaflux%fromStoCO2(nland,nS) &
                        * casaflux%ksoil(nland,nS) &
                        * casapool%csoil(nland,nS)
   ENDDO

   IF(icycle>1) THEN
      DO j=1,mlitter
         casaflux%Nlittermin(nland) = casaflux%Nlittermin(nland) &
                                    + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
      ENDDO
      DO k=1,msoil
         casaflux%Nsmin(nland)   = casaflux%Nsmin(nland)   &
                                    + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
      ENDDO    !gross mineralisation

      DO kk=1,msoil
         DO jj=1,mlitter    ! immobilisation from litter to soil
            casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                                     - casaflux%fromLtoS(nland,kk,jj) &
                                     * casaflux%klitter(nland,jj)     &
                                     * casapool%Clitter(nland,jj)     &
                                     * casapool%ratioNCsoilnew(nland,kk)
         ENDDO
         DO kkk=1,msoil      ! immobilisation from soil to soil
            IF(kkk.ne.kk) THEN
               casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                                        - casaflux%fromStoS(nland,kk,kkk)  &
                                        * casaflux%ksoil(nland,kkk) &
                                        * casapool%Csoil(nland,kkk) &
                                        * casapool%ratioNCsoilnew(nland,kk)
            ENDIF
         ENDDO
      ENDDO  ! immobilization

      casaflux%Nsnet(nland)=casaflux%Nlittermin(nland) &
                                 +casaflux%Nsmin(nland)   &
                                 +casaflux%Nsimm(nland)
                                      ! net mineralization
      IF(casapool%Nsoilmin(nland)>2.0.AND.casamet%tsoilavg(nland)>273.12) THEN
        casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                                   * MAX(0.0,casaflux%Nsnet(nland))
        casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                                   * MAX(0.0,casapool%Nsoilmin(nland))
      ELSE
        casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                                   * MAX(0.0,casaflux%Nsnet(nland)) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
        casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
      ENDIF
      DO k=1,msoil
         DO j=1,mlitter
            casaflux%FluxNtosoil(nland,k) =  casaflux%FluxNtosoil(nland,k)  &
                                 + casaflux%fromLtoS(nland,k,j) &
                                 * casaflux%klitter(nland,j)    &
                                 * casapool%Clitter(nland,j)    &
                                 * casapool%ratioNCsoilnew(nland,k)
         ENDDO  ! end of "j"
         DO kk=1,msoil
            IF(kk.ne.k) THEN
               casaflux%FluxNtosoil(nland,k) = casaflux%FluxNtosoil(nland,k)  &
                                    + casaflux%fromStoS(nland,k,kk) &
                                    * casaflux%ksoil(nland,kk)      &
                                    * casapool%Csoil(nland,kk)      &
                                    * casapool%ratioNCsoilnew(nland,k)
            ENDIF
         ENDDO ! end of "kk"
      ENDDO    ! end of "k"

   ENDIF !end of icycle >1

   IF(icycle >2) THEN
      DO j=1,mlitter
         casaflux%Plittermin(nland) = casaflux%Plittermin(nland) &
                                    + casaflux%klitter(nland,j) * casapool%Plitter(nland,j)
      ENDDO
      DO k=1,msoil
         casaflux%Psmin(nland)   = casaflux%Psmin(nland)   &
                                    + casaflux%ksoil(nland,k)   * casapool%Psoil(nland,k)
      ENDDO    !gross mineralisation

      DO kk=1,msoil
         DO jj=1,mlitter    ! immobilisation from litter to soil
           casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                                     - casaflux%fromLtoS(nland,kk,jj) &
                                     * casaflux%klitter(nland,jj)     &
                                     * casapool%Nlitter(nland,jj)     &
                                     / casapool%ratioNPsoil(nland,kk)
         ENDDO
         DO kkk=1,msoil      ! immobilisation from soil to soil
            IF(kkk.ne.kk) THEN
               casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                                        - casaflux%fromStoS(nland,kk,kkk)  &
                                        * casaflux%ksoil(nland,kkk) &
                                        * casapool%Nsoil(nland,kkk) &
                                        / casapool%ratioNPsoil(nland,kk)
            ENDIF
         ENDDO
      ENDDO  ! immobilization

      casaflux%Psnet(nland)=casaflux%Plittermin(nland) &
                                 +casaflux%Psmin(nland)   &
                                 +casaflux%Psimm(nland)
                                      ! net mineralization


      casaflux%Pleach(nland)  =  casaflux%fPleach(nland) &
                                 * max(0.0,casapool%Psoillab(nland))

      DO k=1,msoil
         DO j=1,mlitter
            casaflux%FluxPtosoil(nland,k) =  casaflux%FluxPtosoil(nland,k)  &
                                 + casaflux%fromLtoS(nland,k,j) &
                                 * casaflux%klitter(nland,j)    &
                                 * casapool%Nlitter(nland,j)    &
                                 / casapool%ratioNPsoil(nland,k)
         ENDDO  ! end of "j"
         DO kk=1,msoil
            IF(kk.ne.k) THEN
               casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                                    + casaflux%fromStoS(nland,k,kk) &
                                    * casaflux%ksoil(nland,kk)      &
                                    * casapool%Nsoil(nland,kk)      &
                                    / casapool%ratioNPsoil(nland,k)
            ENDIF
         ENDDO ! end of "kk"
      ENDDO    ! end of "k"
! need to account for flow from sorbed to occluded pool
   ENDIF
ENDIF  ! end of /=icewater
ENDDO  ! end of nland

DO nland=1,mp
IF(casamet%iveg2(nland)/=icewater) THEN
   casapool%dClitterdt(nland,:) =  casaflux%fluxCtolitter(nland,:) - casaflux%klitter(nland,:) * casapool%clitter(nland,:)
   casapool%dCsoildt(nland,:)   =  casaflux%fluxCtosoil(nland,:)   - casaflux%ksoil(nland,:)   * casapool%csoil(nland,:)
   casaflux%Crsoil(nland)       =  casaflux%fluxCtoCO2(nland)
   if( .NOT. Ticket200 ) &
     casaflux%cnep(nland)         =  casaflux%cnpp(nland) - casaflux%Crsoil(nland)

   IF(icycle > 1) THEN
      casapool%dNlitterdt(nland,:) =  casaflux%fluxNtolitter(nland,:)  &
                                   - casaflux%klitter(nland,:) &
                                   * max(0.0,casapool%Nlitter(nland,:))

      casapool%dNsoildt(nland,:) = casaflux%FluxNtosoil(nland,:) &
                                 - casaflux%ksoil(nland,:) * casapool%Nsoil(nland,:)
      casapool%dNsoilmindt(nland)= casaflux%Nsnet(nland)&
                                 + casaflux%Nmindep(nland) + casaflux%Nminfix(nland)   &
                                 - casaflux%Nminloss(nland)   &
                                 - casaflux%Nminleach(nland)   &
                                 - casaflux%Nupland(nland)

   ENDIF

   IF(icycle >2) THEN
      fluxptase(nland) =  casabiome%prodptase( veg%iveg(nland) ) * deltcasa    &
                       * max( 0.0_r_2, ( casapool%Psoil(nland,2)                   &
                                      * casaflux%ksoil(nland,2)                &
                                      + casapool%Psoil(nland,3)                &
                                      * casaflux%ksoil(nland,3) )              &
                             )                                                 &
                        * max( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )    &
                                      - 15.0 )                                 &
                             )                                                 &
                        / ( max( 0.0_r_2, ( casabiome%costNPup( veg%iveg(nland) )  &
                                        - 15.0 )                               &
                               ) + 150.0                                       &
                          )

      xdplabsorb(nland) = 1.0+ casaflux%Psorbmax(nland)*casaflux%kmlabp(nland) &
                        /((casaflux%kmlabp(nland)+casapool%Psoillab(nland))**2)
      casapool%dPlitterdt(nland,:) = casaflux%fluxPtolitter(nland,:)  &
                                   - casaflux%klitter(nland,:)                 &
                                   * max(zero,casapool%Plitter(nland,:))

      casapool%dPsoildt(nland,1) = casaflux%FluxPtosoil(nland,1)                        &
                                 - casaflux%ksoil(nland,1) * casapool%Psoil(nland,1)

      casapool%dPsoildt(nland,2) = casaflux%FluxPtosoil(nland,2)                        &
                                 - casaflux%ksoil(nland,2) * casapool%Psoil(nland,2)    &
                                 - fluxptase(nland) * casaflux%ksoil(nland,2)*casapool%Psoil(nland,2) &
                                  /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

      casapool%dPsoildt(nland,3) = casaflux%FluxPtosoil(nland,3)                        &
                                 - casaflux%ksoil(nland,3) * casapool%Psoil(nland,3)    &
                                 - fluxptase(nland) * casaflux%ksoil(nland,3)*casapool%Psoil(nland,3) &
                                  /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

      casapool%dPsoillabdt(nland)= casaflux%Psnet(nland) + fluxptase(nland)         &
                                 + casaflux%Pdep(nland) + casaflux%Pwea(nland)      &
                                 - casaflux%Pleach(nland)-casaflux%pupland(nland)   &
                                 - casaflux%kpsorb(nland)*casapool%Psoilsorb(nland) &
                                 + casaflux%kpocc(nland) * casapool%Psoilocc(nland)
      ! here the dPsoillabdt =(dPsoillabdt+dPsoilsorbdt)
      ! dPsoilsorbdt  = xdplabsorb
      casapool%dPsoillabdt(nland)  = casapool%dPsoillabdt(nland)/xdplabsorb(nland)
      casapool%dPsoilsorbdt(nland) = 0.0

      casapool%dPsoiloccdt(nland)  = casaflux%kpsorb(nland)* casapool%Psoilsorb(nland) &
                                   - casaflux%kpocc(nland) * casapool%Psoilocc(nland)
      ! P loss to non-available P pools
!      casaflux%Ploss(nland)        = casaflux%kpocc(nland) * casapool%Psoilocc(nland)

!      casaflux%Ploss(nland)       = casaflux%fPleach(nland) &
!                                 * max(zero,casapool%Psoillab(nland))
      casaflux%Ploss(nland)       = 0.0
   ENDIF
ENDIF
ENDDO

END SUBROUTINE casa_delsoil


