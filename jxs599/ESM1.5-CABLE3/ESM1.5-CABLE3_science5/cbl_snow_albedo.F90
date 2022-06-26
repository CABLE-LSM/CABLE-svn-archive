MODULE cbl_snow_albedo_module

  IMPLICIT NONE

  PUBLIC surface_albedosn
  PRIVATE

CONTAINS

  SUBROUTINE surface_albedosn( AlbSnow, AlbSoil, mp, nrb, surface_type, soil_type, &
                            SnowDepth, SnowDensity, SoilTemp, SnowAge, coszen, metTk )

    USE cable_common_module, ONLY: cable_runtime   
    IMPLICIT NONE
   
    !re-decl input args
    INTEGER, INTENT(IN) :: mp
    INTEGER, INTENT(IN) :: nrb
    REAL, INTENT(OUT)   :: AlbSnow(mp,nrb) 
    REAL, INTENT(IN)    :: AlbSoil(mp,nrb) !becomes INOUT w soilColour param^n
    REAL, INTENT(IN)    :: coszen(mp) 
    REAL, INTENT(IN)    :: SnowDepth(mp)
    REAL, INTENT(IN)    :: SnowDensity(mp)
    REAL, INTENT(IN)    :: SoilTemp(mp)                
    REAL, INTENT(IN)    :: SnowAge(mp)
    INTEGER, INTENT(IN) :: surface_type(mp)          
    INTEGER, INTENT(IN) :: soil_type(mp) 
REAL :: MetTk(mp)                   !Air Temperture at surface - atmospheric forcing (met%tk)

    !working variables  
   REAL, DIMENSION(mp) ::                                                      &
      alv,     &  ! Snow albedo for visible
      alir,    &  ! Snow albedo for near infra-red
      fage,    &  ! age factor
      fzenm,   &  ! zenith factor
      sfact,   &  ! soil factor
      snrat,   &  ! (1-) fraction of soil 'seen' when evaluating surface albedo
      snr,     &  ! when evaluating surface albedo
      tmp         ! temporary value
   
   REAL, PARAMETER ::                                                          &
      alvo  = 0.95,  &  ! albedo for vis. on a new snow
      aliro = 0.70      ! albedo for near-infr. on a new snow
   
    ! local vars    
    REAL :: SoilAlbsoilf(mp) 
    REAL :: Albsoilf_min(mp) 

    !hard wired indexes to be substituted with arg list or USEd from module
    INTEGER, PARAMETER :: perm_ice = 9
    INTEGER, PARAMETER :: lake = 16
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    INTEGER :: i    !looping variable

    !initialise to the no-snow value for albedo for all land points
   SoilAlbsoilF = Albsoil(:,1)

    IF( cable_runtime%esm15_albedo ) THEN
   Albsoilf_min = MetTk   
    ELSE  
   Albsoilf_min = SoilTemp
    ENDIF

    ! lakes - with/without snow cover
    WHERE( surface_type == lake )                                                     
      SoilAlbsoilf = -0.022*( MIN( 275., MAX( 260., Albsoilf_min ) ) - 260. ) + 0.45
    END WHERE
    WHERE(SnowDepth > snow_depth_thresh .and. surface_type == lake )
      SoilAlbsoilF = 0.85
    END WHERE

    sfact(:) = 0.68
   WHERE (SoilAlbsoilf <= 0.14)
      sfact = 0.5
   ELSEWHERE (SoilAlbsoilf > 0.14 .and. SoilAlbsoilf <= 0.20)
      sfact = 0.62
   END WHERE

    !first estimate of snow-affected surface albedos
    AlbSnow(:,2) = 2.0 * SoilAlbsoilF / (1.0 + sfact)
   AlbSnow(:,1) = sfact * AlbSnow(:,2)

    ! calc soil albedo based on colour - Ticket #27
   !H!IF (calcsoilalbedo) THEN
   !H!   CALL soilcol_albedo(ssnow, soil)
   !H!END IF
  
   snrat=0.
   alir =0.
   alv  =0.

   WHERE ( SnowDepth > 1. )
       
      ! Snow age depends on snow crystal growth, freezing of melt water,
      ! accumulation of dirt and amount of new snow.
      snr = SnowDepth / max (SnowDensity, 200.)
      
      !snrat is how little (as fraction) of the underlying soil 'seen'
      snrat = MIN(1.0, snr/ (snr + 0.1))
       
      !IF( cable_runtime%esm15_albedo ) THEN
       WHERE (soil_type == 9)
         !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
         !dnsnow = max (dnsnow, .5) !increase refreshing of snow in Antarctic
         snrat = 1.
      END WHERE
      !END IF

      !331!!snow age and zenith angle factors
      !331!fage = 1.0 - 1.0 / (1.0 + SnowAge )
      !331!tmp = MAX (0.17365, coszen )
      !331!fzenm = MAX(0.0, MERGE(0.0, 1.5/(1.0+4.0*tmp) - 0.5,tmp>0.5) )
      fage = 1. - 1. / (1. + SnowAge ) !age factor
      tmp = MAX( .17365, coszen )
      !Share fzenm = MAX(0.0, MERGE(0.0, 1.5/(1.0+4.0*tmp) - 0.5,tmp>0.5) )
      fzenm = MAX( 0.0, MERGE( 0.0,                                           &
              ( 1. + 1./2. ) / ( 1. + 2.*2. * tmp ) - 1./2., tmp > 0.5 ) )

      !alv and alir: aged-snow albedo
      tmp = alvo * (1.0 - 0.2 * fage)
      alv = 0.4 * fzenm * (1.0 - tmp) + tmp
      tmp = aliro * (1.0 - 0.5 * fage )

      !IF( cable_runtime%esm15_albedo ) THEN
      ! use dry snow albedo for pernament land ice: hard-wired no to be removed
       WHERE (soil_type == 9)
         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
         tmp = 0.75 * (1. - .5 * fage)
      END WHERE
      !END IF
      
      alir = .4 * fzenm * (1.0 - tmp) + tmp
      !talb = .5 * (alv + alir) ! snow albedo
    
   ENDWHERE        ! snowd > 0
   
   !H!! when it is called from cable_rad_driver (UM)
   !H!! no need to recalculate snage
WHERE (SnowDepth > 1 )
      
   snr = SnowDepth / MAX (SnowDensity, 200.)
      
       WHERE (soil_type == 9)
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
       WHERE (soil_type == 9)
         ! permanent ice: hard-wired number to be removed

         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
         tmp = 0.75 * (1. - .5 * fage)
      
      END WHERE
      
      alir = .4 * fzenm * (1.0 - tmp) + tmp

ENDWHERE        ! snowd > 0
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
    IF( cable_runtime%esm15_albedo ) THEN
      WHERE (soil_type == perm_ice) ! use dry snow albedo: 1=vis, 2=nir
    AlbSnow(:,1) = 0.82
    AlbSnow(:,2) = 0.82
  ENDWHERE
    ELSE  
      WHERE (soil_type == perm_ice) ! use dry snow albedo: 1=vis, 2=nir
    AlbSnow(:,1) = alvo - 0.05 ! al*o = albedo appropriate for new snow 
    AlbSnow(:,2) = aliro - 0.05 ! => here al*o LESS arbitrary aging 0.05
  ENDWHERE
    ENDIF

    RETURN

  END SUBROUTINE surface_albedosn

END MODULE cbl_snow_albedo_module
