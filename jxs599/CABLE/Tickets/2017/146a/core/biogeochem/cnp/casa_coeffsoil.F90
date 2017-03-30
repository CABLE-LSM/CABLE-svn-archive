SUBROUTINE casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
!Ticket 146 SUBROUTINE casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casapool,casaflux,casamet)
!  calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
!  and the transfer coefficients between different pools
!
! inputs:
!     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
!
! outputs:
!     klitter(mp,mlitter):      decomposition rate of litter pool (1/day)
!     ksoil(mp,msoil):          decomposition rate of soil pool (1/day)
!     fromLtoS(mp,mlitter,msoil):  fraction of decomposed litter to soil (fraction)
!     fromStoS(mp,msoil,msoil):    fraction of decomposed soil C to another soil pool (fraction)
!     fromLtoCO2(mp,mlitter):      fraction of decomposed litter emitted as CO2 (fraction)
!     fromStoCO2(mp,msoil):        fraction of decomposed soil C emitted as Co2 (fraction)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(IN) :: xklitter,xksoil
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
!Ticket146
  TYPE (casa_pool),             INTENT(IN)    :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  INTEGER j,k,kk,nland             !i: for plant pool, j for litter, k for soil
!Ticket146
  real(r_2), dimension(mp)             :: cuemet, cuestr,cuecwd
  real, parameter                      :: cnmic=10.0               ! microbial biomass C:N ratio 
  logical :: Ticket146 = .false.
   
  casaflux%fromLtoS(:,:,:)      = 0.0
  casaflux%fromStoS(:,:,:)      = 0.0
                                          ! flow from soil to soil
  DO k = 1, msoil
     casaflux%fromStoS(:,k,k)   = -1.0
  ENDDO  ! "k"
  casaflux%fromLtoCO2(:,:) = 0.0             ! flow from L or S to CO2
  casaflux%fromStoCO2(:,:) = 0.0

  casaflux%klitter(:,:) = 0.0        !initialize klitter (Q.Zhang 03/03/2011)

  WHERE(casamet%iveg2/=icewater)

    casaflux%klitter(:,metb)   = xklitter(:) * casabiome%litterrate(veg%iveg(:),metb)
    casaflux%klitter(:,str)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),str) &
                                 * exp(-3.0*casabiome%fracLigninplant(veg%iveg(:),leaf))
    casaflux%klitter(:,cwd)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),cwd)

    casaflux%ksoil(:,mic)      = xksoil(:) * casabiome%soilrate(veg%iveg(:),mic)   &
                               * (1.0 - 0.75 *(soil%silt(:)+soil%clay(:)))
    casaflux%ksoil(:,slow)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),slow)
    casaflux%ksoil(:,pass)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),pass)
    casaflux%kplab(:)          = xksoil(:) * casabiome%xkplab(casamet%isorder(:))
    casaflux%kpsorb(:)         = xksoil(:) * casabiome%xkpsorb(casamet%isorder(:))
    casaflux%kpocc(:)          = xksoil(:) * casabiome%xkpocc(casamet%isorder(:))


    WHERE(veg%iveg==cropland)      ! for cultivated land type
       casaflux%ksoil(:,mic)  = casaflux%ksoil(:,mic) * 1.25
       casaflux%ksoil(:,slow) = casaflux%ksoil(:,slow)* 1.5
       casaflux%ksoil(:,pass) = casaflux%ksoil(:,pass)* 1.5
    ENDWHERE  !

                                          ! flow from litter to soil
    if(Ticket146) then
      cuemet(:) = 0.45
      cuestr(:) = 0.45
      cuecwd(:) = 0.4  
    else
      cuemet(:) = 0.45
      cuestr(:) = 0.7
      cuecwd(:) = 0.4
    endif
    
    casaflux%fromLtoS(:,mic,metb)  = cuemet(:)                                  
                                          ! metb -> mic
    casaflux%fromLtoS(:,mic,str)   = cuestr(:)*(1.0-casabiome%fracLigninplant(veg%iveg(:),leaf))  
                                          ! str -> mic
    casaflux%fromLtoS(:,slow,str)  = cuestr(:) * casabiome%fracLigninplant(veg%iveg(:),leaf)       
                                          ! str -> slow
    casaflux%fromLtoS(:,mic,cwd)   = cuecwd(:) *(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood)) 
                                          ! CWD -> fmic
    if(Ticket146) then
      casaflux%fromLtoS(:,slow,cwd)  = cuecwd(:) * casabiome%fracLigninplant(veg%iveg(:),wood)        
                                          ! CWD -> slow
    else
      casaflux%fromLtoS(:,slow,cwd)  = cuestr(:) * casabiome%fracLigninplant(veg%iveg(:),wood)        
    endif
    
!! set the following two backflow to set (see Bolker 199x)
!    casaflux%fromStoS(:,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(:))
!    casaflux%fromStoS(:,mic,pass)  = 0.45

    casaflux%fromStoS(:,slow,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
                                     * (0.997 - 0.032*soil%clay(:))
    casaflux%fromStoS(:,pass,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
                                     * (0.003 + 0.032*soil%clay(:))
    casaflux%fromStoS(:,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(:) )

  ENDWHERE

  DO nland=1,mp
    IF(casamet%iveg2(nland)/=icewater) THEN
      DO j=1,mlitter
        DO k=1,msoil
          casaflux%fromLtoCO2(nland,j) = casaflux%fromLtoCO2(nland,j)  &
                                       + casaflux%fromLtoS(nland,k,j)
        ENDDO  !"k"
        casaflux%fromLtoCO2(nland,j) = 1.0 - casaflux%fromLtoCO2(nland,j)
      ENDDO !"j"
      DO k=1,msoil
        DO kk=1,msoil
          casaflux%fromStoCO2(nland,k) = casaflux%fromStoCO2(nland,k) &
                                       + casaflux%fromStoS(nland,kk,k)
        ENDDO  !"kk"
      ENDDO   !"k"
      casaflux%fromStoCO2(nland,:) = -casaflux%fromStoCO2(nland,:)
    ENDIF
  ENDDO   ! "nland"

END SUBROUTINE casa_coeffsoil


