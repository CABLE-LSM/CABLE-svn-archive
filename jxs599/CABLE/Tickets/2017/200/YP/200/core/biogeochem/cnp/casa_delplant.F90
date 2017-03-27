! modified by ypw following Chris Lu 5/nov/2012
SUBROUTINE casa_delplant(veg,casabiome,casapool,casaflux,casamet,            &
     cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

  !  calculate the chnage in plant C, N and P pools
  !  uptake of N and P will be computed in casa_uptake
  !  labile C pool will be computed casa_labile

  IMPLICIT NONE
  TYPE (veg_parameter_type), INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),         INTENT(INOUT) :: casabiome
  TYPE (casa_pool),          INTENT(INOUT) :: casapool
  TYPE (casa_flux),          INTENT(INOUT) :: casaflux
  TYPE (casa_met),           INTENT(INOUT) :: casamet

  ! added by ypwang following Chris Lu 5/nov/2012
  real, dimension(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

  INTEGER  npt,nL,nP,nland, ivt
  real(r_2)      :: Ygrow, ratioPNplant
  logical :: Ticket200 = .false.

  casaflux%FluxCtolitter = 0.0
  casaflux%FluxNtolitter = 0.0
  casaflux%FluxPtolitter = 0.0

  ! added by ypwang following Chris Lu 5/nov/2012

  cleaf2met = 0.0
  cleaf2str = 0.0
  croot2met = 0.0
  croot2str = 0.0
  cwood2cwd = 0.0

  nleaf2met = 0.0
  nleaf2str = 0.0
  nroot2met = 0.0
  nroot2str = 0.0
  nwood2cwd = 0.0

  pleaf2met = 0.0
  pleaf2str = 0.0
  proot2met = 0.0
  proot2str = 0.0
  pwood2cwd = 0.0

  !MPI
  DO npt=1,mp
     IF(casamet%iveg2(npt)/=icewater) THEN
        !    PRINT *, 'npt = ', npt
        !    PRINT *, 'casapool%cplant(npt,:) = ', casapool%cplant(npt,:)
        casapool%dcplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
             - casaflux%kplant(npt,:)  * casapool%cplant(npt,:)

        if( .NOT. Ticket200 ) then
        !! adjust turnover and autotrophic respiration to avoid negative stores.
        !! Ticket#108

        where (((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)).lt. 0.0) &
                  .OR. ((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
                  .lt. 0.5 * casapool%cplant(npt,2:3) ))
           casaflux%kplant(npt,2:3) = 0.0
           casaflux%crmplant(npt,2:3)= 0.0
        endwhere
        IF(ANY((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)).lt. 0.0)) THEN
           casaflux%kplant(npt,1) = 0.0
           casaflux%crmplant(npt,1)= min(casaflux%crmplant(npt,1),0.5*casaflux%Cgpp(npt))
        ENDIF

        !! revise turnover and NPP and dcplantdt to reflect above adjustments
       
        casaflux%Cplant_turnover(npt,:) = casaflux%kplant(npt,:)  * casapool%cplant(npt,:)
        if (any((casapool%dcplantdt(npt,:)*deltpool + casapool%cplant(npt,:)).lt. 0.0) &

          .OR. any((casapool%dcplantdt(npt,2:3)*deltpool + casapool%cplant(npt,2:3)) &
                  .lt. 0.5 * casapool%cplant(npt,2:3) )) then
        
           ratioPNplant = 0.0  
           if (casapool%Nplant(npt,leaf)>0.0) &
                ratioPNplant = casapool%Pplant(npt,leaf)/(casapool%Nplant(npt,leaf)+ 1.0e-10)  

           Ygrow = 0.65+0.2*ratioPNplant/(ratioPNplant+1.0/15.0)
           IF ((casaflux%Cgpp(npt)-SUM(casaflux%crmplant(npt,:)))>0.0) THEN
              ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
              casaflux%crgplant(npt)  = (1.0-Ygrow)* max(0.0,casaflux%Cgpp(npt)- &
                   SUM(casaflux%crmplant(npt,:)))
           ELSE
              casaflux%crgplant(npt) = 0.0
           ENDIF
           casaflux%Cnpp(npt) = casaflux%Cgpp(npt)-SUM(casaflux%crmplant(npt,:)) &
                - casaflux%crgplant(npt) - casaflux%fracClabile(npt) * casaflux%cgpp(npt)

           casapool%dcplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
                - casaflux%kplant(npt,:)  * casapool%cplant(npt,:)

        endif
        endif !Ticket200 
        
        ! change here made by ypw on 26august 2011
        ! calculate fraction c to labile pool as a fraction of gpp, not npp
        ! casapool%dClabiledt(npt)   = casaflux%Cnpp(npt)    * casaflux%fracClabile(npt)
        casapool%dClabiledt(npt)   =  casaflux%Cgpp(npt)  * casaflux%fracClabile(npt) &
             - casaflux%clabloss(npt)
        ! added by ypwang 5/nov/2012
        cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
        cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
        croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
        croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
        cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

        !    PRINT *, 'npt, mp, iveg', npt, mp, veg%iveg(npt)
        IF(icycle > 1) THEN
           !    PRINT *, 'casapool%Nplant(npt,:) = ', casapool%Nplant(npt,:)
           IF(casaflux%fracNalloc(npt,leaf)==0.0) THEN
              casapool%dNplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf)
           else
              casapool%dNplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf) &
                   * casabiome%ftransNPtoL(veg%iveg(npt),leaf)
           ENDIF
           !Ticket200: if
           IF (casamet%lnonwood(npt)==0) THEN
              casapool%dNplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Nplant(npt,wood) &
                                        * casabiome%ftransNPtoL(veg%iveg(npt),wood)
           ELSE
              casapool%dNplantdt(npt,wood) = 0.0
           ENDIF
          
           casapool%dNplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Nplant(npt,froot) &
                * casabiome%ftransNPtoL(veg%iveg(npt),froot)
           ! added by ypwang 5/nov/2012

           nleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix
           nroot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix

           nleaf2met(npt) = - casapool%dNplantdt(npt,leaf)  - nleaf2str(npt)
           nroot2met(npt) = - casapool%dNplantdt(npt,froot) - nroot2str(npt)

           nwood2cwd(npt) = -casapool%dNplantdt(npt,wood)

        ENDIF


        IF(icycle >2) THEN

           IF(casaflux%fracPalloc(npt,leaf)==0.0) THEN
              casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf)
           else
              casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf) &
                   * casabiome%ftransPPtoL(veg%iveg(npt),leaf)
           ENDIF

           casapool%dPplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Pplant(npt,wood) &
                * casabiome%ftransPPtoL(veg%iveg(npt),wood)

          
           casapool%dPplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Pplant(npt,froot) &
                * casabiome%ftransPPtoL(veg%iveg(npt),froot)
           ! added by ypwang 5/nov/2012

           pleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix
           proot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
           pleaf2met(npt) = -casapool%dPplantdt(npt,leaf)  - pleaf2str(npt)
           proot2met(npt) = -casapool%dPplantdt(npt,froot) - proot2str(npt)
           pwood2cwd(npt) = -casapool%dPplantdt(npt,wood)


        ENDIF

        DO nL=1,mlitter
           DO nP=1,mplant
              casaflux%FluxCtolitter(npt,nL) = casaflux%FluxCtolitter(npt,nL) &
                   + casaflux%fromPtoL(npt,nL,nP) &
                   * casaflux%kplant(npt,nP) &
                   * casapool%cplant(npt,nP)
           ENDDO
        ENDDO

        !    PRINT *, 'before 2nd icycle >1; npt, mp', npt, mp
        IF(icycle > 1) THEN
           !Ticket#108 to avoid -ve Nitrogen pools 
           casaflux%FluxNtolitter(npt,str) = min( &
                 casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)   &
               * casapool%cplant(npt,leaf) * ratioNCstrfix ,                   &
                   (-)casapool%dNplantdt(npt,leaf))               &
                  + min( casaflux%fromPtoL(npt,str,froot)         &
                  * casaflux%kplant(npt,froot)                    &
                  * casapool%cplant(npt,froot) * ratioNCstrfix  , &
                    (-)casapool%dNplantdt(npt,froot) )

           casaflux%FluxNtolitter(npt,metb) = - casapool%dNplantdt(npt,leaf)-casapool%dNplantdt(npt,froot) &
                - casaflux%FluxNtolitter(npt,str)
           casaflux%FluxNtolitter(npt,CWD) = -casapool%dNplantdt(npt,wood)

           ! adding N uptake
           casapool%dNplantdt(npt,:) = casapool%dNplantdt(npt,:) &
                + casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,:)
! now accounted for in delsoil
           !       casapool%Nsoilmin(npt)    = casapool%Nsoilmin(npt) - casaflux%Nminuptake(npt) *deltpool
        ENDIF !end "icycle >1"

       
        IF(icycle>2) THEN
           casaflux%FluxPtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix        &
                + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
           casaflux%FluxPtolitter(npt,metb) = -casapool%dPplantdt(npt,leaf)-casapool%dPplantdt(npt,froot) &
                - casaflux%FluxPtolitter(npt,str)
           casaflux%FluxPtolitter(npt,CWD) = -casapool%dPplantdt(npt,wood)
           ! add P uptake
           casapool%dPplantdt(npt,:) = casapool%dPplantdt(npt,:) &
                + casaflux%Plabuptake(npt)*casaflux%fracPalloc(npt,:)
           !       casapool%Psoillab(npt)    = casapool%Psoillab(npt) - casaflux%Plabuptake(npt) * deltpool
        ENDIF  !of "icycle >2"
    

     ENDIF
  ENDDO

END SUBROUTINE casa_delplant


