MODULE cbl_snow_albedo_module

  IMPLICIT NONE

  PUBLIC surface_albedosn
  PRIVATE

CONTAINS

  SUBROUTINE surface_albedosn( AlbSnow, AlbSoil, mp, nrb, jls_radiation, surface_type, soil_type, &
                            SnowDepth, SnowODepth, SnowFlag_3L, SnowDensity, &
                            SoilTemp, SnowTemp, SnowAge, & 
                            metTk, coszen )
    !H!jhan:Eliminate these USE data statements  
    !USE cable_common_module, ONLY : kwidth_gl     !not needed anymore
    !USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ  !not needed anymore

    IMPLICIT NONE

    !re-decl input args
    INTEGER, INTENT(IN) :: mp
    INTEGER, INTENT(IN) :: nrb
    LOGICAL, INTENT(IN) :: jls_radiation             !runtime switch def. in cable_*main routines - not needed anymore
    REAL, INTENT(OUT)   :: AlbSnow(mp,nrb) 
    REAL, INTENT(IN)    :: AlbSoil(mp,nrb)           !NB to become IN OUT because of soil colour parameterization
    REAL, INTENT(IN)    :: MetTk(mp)                 !not needed anymore
    REAL, INTENT(IN)    :: coszen(mp) 
    REAL, INTENT(IN)    :: SnowDepth(mp)
    REAL, INTENT(IN)    :: SnowODepth(mp)            !not needed anymore
    REAL, INTENT(IN)    :: SnowDensity(mp)
    REAL, INTENT(IN)    :: SoilTemp(mp)                
    REAL, INTENT(IN)    :: SnowTemp(mp)              !not needed anymore
    REAL, INTENT(IN)    :: SnowAge(mp)
    INTEGER, INTENT(IN) :: SnowFlag_3L(mp)           !not needed anymore
    INTEGER, INTENT(IN) :: surface_type(mp)          
    INTEGER, INTENT(IN) :: soil_type(mp) 

    !working variables  
    REAL, DIMENSION(mp) ::                                                      &
         alv,     &  ! Snow albedo for visible
         alir,    &  ! Snow albedo for near infra-red
         !ar1,     &  ! crystal growth  (-ve)            !not needed anymore
         !ar2,     &  ! freezing of melt water           !not needed anymore
         !ar3,     &  !                                  !not needed anymore
         !dnsnow,  &  ! new snow albedo                  !not needed anymore
         !dtau,    &  !                                  !not needed anymore
         fage,    &  ! age factor
         fzenm,   &  !
         sfact,   &  !
         snr,     &  !
         snrat,   &  !
         !talb,    &  ! snow albedo                      !legacy - not used
         tmp,     &  ! temporary value
         SoilAlbsoilF

    REAL, PARAMETER ::                                                          &
         alvo  = 0.95,  &  ! albedo for vis. on a new snow
         aliro = 0.70      ! albedo for near-infr. on a new snow

    !hard wired indexes to be substituted with arg list or USEd from module
    INTEGER, PARAMETER :: perm_ice = 9
    INTEGER, PARAMETER :: lake = 16
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    INTEGER :: i    !looping variable

    !initialise to the no-snow value for albedo for all land points
    SoilAlbsoilF = Albsoil(:,1)

    ! lakes - with/without snow cover
    WHERE( surface_type == lake )                                                     
       SoilAlbsoilF = -0.022*( MIN( 275.0, MAX( 260.0, SoilTemp) ) - 260.0 ) + 0.45
    END WHERE
    WHERE(SnowDepth > snow_depth_thresh .and. surface_type == lake )
       SoilAlbsoilF = 0.85
    END WHERE

    sfact(:) = 0.68
    WHERE (SoilAlbsoilF <= 0.14)
       sfact = 0.5
    ELSEWHERE (SoilAlbsoilF > 0.14 .and. SoilAlbsoilF <= 0.20)
       sfact = 0.62
    END WHERE

    !first estimate of snow-affected albedos
    AlbSnow(:,2) = 2.0 * SoilAlbsoilF / (1.0 + sfact)
    AlbSnow(:,1) = sfact * AlbSnow(:,2)

    ! calc soil albedo based on colour - Ticket #27
    !H!IF (calcsoilalbedo) THEN
    !H!   CALL soilcol_albedo(ssnow, soil)
    !H!END IF

    !no snow values for working variables.
    snrat(:)=0.0
    alir(:) =0.0
    alv(:)  =0.0

    !Ticket 331 - snow aging evaluation moved to soilsnow routine
    !NB since soil_type==perm_ice is overwritten at the end
    !these conditions/that code can be removed; 
    !the outer DO-IF can then be collapsed back to a WHERE loop
    ! with the remaining IF replaced with a MERGE function
    !for fzenm as per the original code
    DO i = 1,mp
       IF (SnowDepth(i) > snow_depth_thresh) THEN

          snr(i) = SnowDepth(i) / MAX (SnowDensity(i), 200.0)

          IF (soil_type(i) == perm_ice) THEN
             snrat(i) = 1.0
          ELSE
             snrat(i) = MIN (1.0, snr(i) / (snr(i) + 0.1) )
          END IF
        
          !snow age and zenith angle factors
          fage(i) = 1.0 - 1.0 / (1.0 + SnowAge(i) ) !age factor
          tmp(i) = MAX (0.17365, coszen(i) )
          !Ticket 331 - originally as MAX(0.0,MERGE(0.0,function,tmp>0.5))
          IF (tmp(i)>0.5) THEN
             fzenm(i) = 0.0
          ELSE
             !Ticket 331 - why isn't this written as
             !fzenm(i) = 1.5/(1.0+4.0*tmp(i)) - 0.5
             fzenm(i) = ( 1.0 + 1.0/2.0 ) / ( 1.0 + 2.0* 2.0 * tmp(i) ) - 1.0/2.0
          END IF
          fzenm(i) = MAX( 0.0, fzenm(i))
           
          
          tmp(i) = alvo * (1.0 - 0.2 * fage(i))
          alv(i) = 0.4 * fzenm(i) * (1.0 - tmp(i)) + tmp(i)
          tmp(i) = aliro * (1.0 - 0.5 * fage(i))

          ! use dry snow albedo on permanent ice points
          IF (soil_type(i) == perm_ice) THEN
             tmp(i) = 0.95 * (1.0 - 0.2 * fage(i))
             alv(i) = 0.4 * fzenm(i) * (1.0 - tmp(i)) + tmp(i)
             tmp(i) = 0.75 * (1.0 - 0.5 * fage(i))
          END IF

          alir(i) = 0.4 * fzenm(i) * (1.0 - tmp(i)) + tmp(i)
          !talb(i) = 0.5 * (alv(i) + alir(i)) ! snow albedo - not needed as not used!!

       END IF
    END DO        ! snowd > snow_deth_thresh
    
    !H!jhan:SLI currently not available
    !H!IF(cable_user%SOIL_STRUC=='sli') THEN
    !H!   WHERE (SnowDepth.GT.snow_depth_thresh)
    !H!      snrat = 1.0   ! using default parameterisation, albedo is too low,
    !H!      ! inhibiting snowpack initiation
    !H!   ENDWHERE
    !H!ENDIF

    !final values of soil-snow albedos - 1=vis, 2=nir
    AlbSnow(:,2) = MIN( aliro,                                          &
                          ( 1.0 - snrat ) * AlbSnow(:,2) + snrat * alir)

    AlbSnow(:,1) = MIN( alvo,                                           &
                          ( 1.0 - snrat ) * AlbSnow(:,1) + snrat * alv )

    !except for ice regions
    WHERE (soil_type == perm_ice) ! use dry snow albedo: 1=vis, 2=nir
       AlbSnow(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow 
       AlbSnow(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
    END WHERE

    RETURN

  END SUBROUTINE surface_albedosn

END MODULE cbl_snow_albedo_module
