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
USE cable_common_module, ONLY : kwidth_gl     !not needed anymore
use cable_phys_constants_mod, ONLY : CTFRZ => TFRZ

implicit none

!re-decl input args
integer :: mp
integer :: nrb
LOGICAL :: jls_radiation            !runtime switch def. in cable_*main routines 
REAL :: AlbSnow(mp,nrb) 
REAL :: AlbSoil(mp,nrb) 
REAL :: MetTk(mp)                   !not needed anymore
REAL :: coszen(mp) 
REAL :: SnowDepth(mp)
REAL :: SnowODepth(mp)              !not needed anymore
REAL :: SnowDensity(mp)
REAL :: SoilTemp(mp)                !not needed anymore
REAL :: SnowTemp(mp)                !not needed anymore
REAL :: SnowAge(mp)
integer:: SnowFlag_3L(mp)           !not needed anymore
integer:: surface_type(mp)          
integer:: soil_type(mp) 


    REAL, DIMENSION(mp) ::                                                      &
         alv,     &  ! Snow albedo for visible
         alir,    &  ! Snow albedo for near infra-red
         ar1,     &  ! crystal growth  (-ve)            !not needed anymore
         ar2,     &  ! freezing of melt water           !not needed anymore
         ar3,     &  !                                  !not needed anymore
         dnsnow,  &  ! new snow albedo                  !not needed anymore
         dtau,    &  !                                  !not needed anymore
         fage,    &  ! age factor
         fzenm,   &  !
         sfact,   &  !
         snr,     &  !
         snrat,   &  !
         talb,    &  ! snow albedo                      !legacy - not used
         tmp         ! temporary value

    REAL, PARAMETER ::                                                          &
         alvo  = 0.95,  &  ! albedo for vis. on a new snow
         aliro = 0.70      ! albedo for near-infr. on a new snow

    !hard wired indexes to be subistituted with arg list or USE from moduled
    INTEGER, PARAMETER :: perm_ice = 9
    INTEGER, PARAMETER :: lake = 16
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    ! local vars    
    REAL :: SoilAlbsoilf(mp) 

    SoilAlbsoilF = Albsoil(:,1)

    ! lakes: hard-wired number to be removed in future
    WHERE( surface_type == lake )                                                     
       SoilAlbsoilF = -0.022*( MIN( 275., MAX( 260., SoilTemp) ) - 260. ) + 0.45
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

    !first estimate of albedos
    AlbSnow(:,2) = 2. * SoilAlbsoilF / (1. + sfact)
    AlbSnow(:,1) = sfact * AlbSnow(:,2)

    ! calc soil albedo based on colour - Ticket #27
    !H!IF (calcsoilalbedo) THEN
    !H!   CALL soilcol_albedo(ssnow, soil)
    !H!END IF

    !no snow values for working variables.
    snrat(:)=0.0
    alir(:) =0.0
    alv(:)  =0.0

    !Snow aging evaluation moved to soilsnow routine
    WHERE (SnowDepth > snow_depth_thresh) ! .and. jls_radiation )

       snr = SnowDepth / MAX (SnowDensity, 200.)

       WHERE (soil_type == perm_ice)
          ! permanent ice: hard-wired number to be removed
          snrat = 1.
       ELSEWHERE
          snrat = MIN (1., snr / (snr + .1) )
       END WHERE

       fage = 1. - 1. / (1. + SnowAge ) !age factor
       tmp = MAX (.17365, coszen )
       fzenm = MAX( 0., MERGE( 0.0,                                             &
            ( 1. + 1./2. ) / ( 1. + 2. * 2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)

       ! use dry snow albedo
       WHERE (soil_type == perm_ice)
          ! permanent ice: hard-wired number to be removed

          tmp = 0.95 * (1.0 - 0.2 * fage)
          alv = .4 * fzenm * (1. - tmp) + tmp
          tmp = 0.75 * (1. - .5 * fage)

       END WHERE

       alir = .4 * fzenm * (1.0 - tmp) + tmp
       talb = .5 * (alv + alir) ! snow albedo - not needed as not used!!

    END WHERE        ! snowd > snow_deth_thresh=
    
    !H!jhan:SLI currently not available
    !H!IF(cable_user%SOIL_STRUC=='sli') THEN
    !H!   WHERE (SnowDepth.GT.1.0)
    !H!      snrat = 1.0   ! using default parameterisation, albedo is too low,
    !H!      ! inhibiting snowpack initiation
    !H!   ENDWHERE
    !H!ENDIF

    !final values of soil-snow albedos
    AlbSnow(:,2) = MIN( aliro,                                          &
                          ( 1. - snrat ) * AlbSnow(:,2) + snrat * alir)

    AlbSnow(:,1) = MIN( alvo,                                           &
                          ( 1. - snrat ) * AlbSnow(:,1) + snrat * alv )

    !ecept for ice regions
    WHERE (soil_type == perm_ice) ! use dry snow albedo: 1=vis, 2=nir
       AlbSnow(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow 
       AlbSnow(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
    END WHERE

    RETURN

END SUBROUTINE surface_albedosn

END MODULE cbl_snow_albedo_module
