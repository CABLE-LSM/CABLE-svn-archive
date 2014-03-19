! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine FTSA -------------------------------------------------
!
!    It calculates (true) surface albedos for sea and sea ice
!    Release 2.8 of the UM allows for separate surface
!    albedos for direct and diffuse light over sea.
!
!    Suitable for single column model use.
!
!    Author: William Ingram
!
!    Offline documentation is in UMDP 23, sections "True surface albedo
!    specification" and "Modifications to the radiation scheme to
!    accommodate the leads model"
!
!    There are two options for the sea ice albedo scheme:
!    1) the default albedo parameterisation of the Los Alamos sea ice
!    model CICE vn4.1 (Hunke and Lipscomb, 2010) which was used in CCSM3.
!    This calculates the albedo for 4 bands, visible/near-IR
!    diffuse/direct.
!    2) the original UM scheme which simply uses 1 band.
!
SUBROUTINE ftsa (                                                              &
  flandg, aice, aice_cat, tstar, tstar_sice_cat, cosz, ws_10m, chloro,         &
  s, s_sea_cat, di_cat,                                                        &
  alpham,alphac,alphab,dtice,                                                  &
  dt_bare,dalb_bare_wet,pen_rad_frac,beta,                                     &
  n_points, n_band_dim, n_bands, nice_use, spec_band_bot, spec_band_top,       &
  sa_sice, saos)

USE switches, ONLY: i_sea_alb_method, l_spec_sea_alb, l_sea_alb_var_chl,       &
                    l_ssice_albedo, l_sice_meltponds, l_sice_scattering,       &
                    l_sice_hadgem1a, l_cice_alb
USE c_0_dg_c
USE fluxes, ONLY: alb_sice
USE ancil_info, ONLY: sea_pts, sea_index, sice_pts, ssi_index, sice_pts_ncat,  &
                      sice_index_ncat, sice_frac_ncat
USE c_kappai, ONLY : kappai,kappai_snow,rhosnow
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE Jin11_osa_mod, ONLY: Jin11_osa

#if defined(UM_JULES)
USE ereport_mod, ONLY : ereport
#endif

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
     n_points,                                                                 &
!      Full field dimension
     n_band_dim,                                                               &
!      Dimension of the sea albedo in bands (to fit all SW calls)
     n_bands,                                                                  &
!      Number of SW radiation bands in this call
     nice_use
!      Number of sea ice categories used in radiation scheme

REAL, INTENT(IN) ::                                                            &
     flandg(n_points),                                                         &
                                    ! Land fraction
     aice(n_points),                                                           &
                                    ! Sea-ice fraction
     aice_cat(n_points,nice_use),                                              &
                                    ! Category sea-ice fraction
     tstar(n_points),                                                          &
                                    ! Surface temperature
     tstar_sice_cat(n_points,nice_use),                                        &
                                    ! Category seaice surface temperature
     cosz(n_points),                                                           &
                                    ! cos(solar zenith angle)
     ws_10m(n_points),                                                         &
                                    ! 10m wind speed
     chloro(n_points),                                                         &
                                    ! nr surface chlorophyll content
     s(n_points),                                                              &
                                    ! Snow amount (mass/area)
     s_sea_cat(n_points,nice_use),                                             &
                     ! Category snow amount on sea ice (mass/area of ice)

     di_cat(n_points, nice_use),                                               &
                     ! Effective thickness of sea ice categories
     spec_band_bot(n_bands), spec_band_top(n_bands)
                     ! boundaries of the spectral bands


REAL, INTENT(IN) ::                                               &
!
!     Constants used to determine the albedo of sea-ice:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
!
     alpham                                                       &
                       ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
!
    ,alphac                                                       &
                       ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
    ,alphab                                                       &
                       ! "B" for "bare"
!
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
    ,dtice                                                        &
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
    ,dt_bare                                                      &
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
    ,dalb_bare_wet                                                &
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
    ,pen_rad_frac                                                 &
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
   ,beta

REAL, INTENT(OUT) ::                                                           &
     sa_sice(n_points, 4),                                                     &
!      Surface Albedo for sea ice
!             (*,1) - direct beam visible
!             (*,2) - diffuse visible
!             (*,3) - direct beam near-ir
!             (*,4) - diffuse near-ir
     saos(n_points, 2, n_band_dim)
!      Surface albedo for Open Sea (direct and diffuse components, for each
!        band, with zeros for safety where no value applies)

! Local variables
! reduce some inputs on n_points, onto sea-points:
REAL, ALLOCATABLE :: cosz_s(:)
REAL, ALLOCATABLE :: ws_10m_s(:)
REAL, ALLOCATABLE :: chloro_s(:)
REAL, ALLOCATABLE :: saos_s(:,:,:)

INTEGER band, n, j, l, ll           ! Loops over points
REAL dsa                            ! Deep-snow albedo (alphasubD)
REAL bx                                                           &
                                    ! 1 - COSZ
! Temperature at which (snow on) sea-ice reaches its "cold" value
    ,tcice                                                        &
! Slope and intercept of temperature dependence of the albedo of
! (snow on) sea-ice
    ,ice1,ice2

! Local parameters

REAL, PARAMETER :: adifc = 0.06
!      Surface albedo of ice-free sea for the diffuse beam
REAL, PARAMETER :: chloro_constant = 5.0E-7
!      Chlorophyll content of the sea near-surface, if constant (0.5mg m-3)

REAL, PARAMETER :: maskd = 0.2
!      Masking depth (S in 3.6.1)

! Parameters for the CICE scheme:
REAL, PARAMETER :: albicev = 0.78
!      Visible ice albedo for h > ahmax and T < -1C
REAL, PARAMETER :: albicei = 0.36
!      Near-ir ice albedo for h > ahmax and T < -1C
REAL, PARAMETER :: albsnowv = 0.98
!      Visible snow albedo for cold snow (T < -1C)
REAL, PARAMETER :: albsnowi = 0.70
!      Near-ir snow albedo for cold snow (T < -1C)
REAL, PARAMETER :: ahmax = 0.3
!      Albedo is constant above this thickness (metres)
REAL, PARAMETER :: dalb_mlt  = -0.075
!      Ice albedo change per 1 degree change in temp for ice (for -1C to 0C)
REAL :: dalb_mlts(4)
DATA dalb_mlts  /-0.1,-0.1,-0.15,-0.15/
!      Snow albedo change per 1 degree change in temp (for -1C to 0C)
!       (values for each radiation band)
REAL, PARAMETER :: snowpatch = 0.02
!      Length scale for parameterizing non uniform snow coverage (m)


REAL :: fh, albo, albice(4), albsnow(4), fT, albs(4), area_snow
REAL :: fhtan

REAL :: snow_albedo
!      Snow albedo

REAL :: hice(n_points,nice_use)
!      Ice thickness

REAL :: saos_bb(n_points,2)
!      Broadband open sea albedo, for methods that need it.
REAL :: spectralsea(n_bands)
!      Spectral albedo scaling factor, for methods that need it.
REAL, PARAMETER :: spec_band_bb(1) = -1.0
!      Array to use as sop and bottom of spectral bands, for Jin11 broadband

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ::              errcode            ! Error code
CHARACTER (LEN = 80) :: errmsg             ! Error message text

IF (lhook) CALL dr_hook('FTSA',zhook_in,zhook_handle)

fhtan = ATAN(ahmax*4.0)

! Open Sea albedo
IF (i_sea_alb_method == 1) THEN
  ! 1 - Briegleb and Ramanathan, 1982, J. Appl. Met.
  ! (doi:10.1175/1520-0450(1982)021<1160:SADVIC>2.0.CO;2)
  saos_bb(:,:) = 0.0
  DO j=1, n_points
    IF (flandg(j) <  1.0) THEN
      saos_bb(j, 1) = 0.026/(cosz(j)**1.7+0.065)                     &
          +0.15*(cosz(j)-0.1)*(cosz(j)-0.5)*(cosz(j)-1.0)
      saos_bb(j,2) = adifc
    END IF
  END DO

  saos = 0.0
  ! this method can only apply spectral scaling if there are 6 bands
  ! as per HadGEM2, as the scaling is hardcoded. This is kept here for 
  ! consistency with the UM radiaiton code.
  IF (l_spec_sea_alb .AND. n_bands == 6) THEN
    ! copy the bb albedo onto bands, with a scaling factor
    spectralsea = (/ 0.000000, 1.444205, 1.799420,                 &
                         0.5470291, 0.000000, 0.000000 /)
    ! copy the broadband albedo onto all bands:
    DO band=1, n_bands
      DO j=1, n_points
        saos(j, 1, band) = saos_bb(j, 1) * spectralsea(band)
        saos(j, 2, band) = saos_bb(j, 2) * spectralsea(band)
      END DO
    END DO
  ELSE
    ! copy the broadband albedo onto all bands:
    DO band=1, n_bands
      DO j=1, n_points
        saos(j, 1, band) = saos_bb(j, 1)
        saos(j, 2, band) = saos_bb(j, 2)
      END DO
    END DO
  END IF

ELSE IF (i_sea_alb_method == 2) THEN
  ! 2 - Modified Barker and Li, 1995, J. Climate,
  ! (doi:10.1175/1520-0442(1995)008<2213:ISOCSS>2.0.CO;2)
  ! modified to agree with Met Office aircraft obs for 
  ! water clarity around the UK
  saos_bb(:,:) = 0.0
  DO j=1, n_points
    IF (flandg(j) <  1.0) THEN
      bx=1.0-cosz(j)
      saos_bb(j, 1) = 0.0315 + 0.0421*bx**2                       &
             + 0.128*bx**3 - 0.04*bx**4                           &
             + (3.12/(5.68) + (0.074*bx/(1.0)))*bx**5
      saos_bb(j,2) = adifc
    END IF
  END DO

  saos = 0.0
  ! this method can only apply spectral scaling if there are 6 bands
  ! as per HadGEM2, as the scaling is hardcoded. This is kept here for 
  ! consistency with the UM radiaiton code.
  IF (l_spec_sea_alb .AND. n_bands == 6) THEN
    ! copy the bb albedo onto bands, with a scaling factor
    spectralsea = (/ 0.000000, 1.444205, 1.799420,                 &
                         0.5470291, 0.000000, 0.000000 /)
    ! copy the broadband albedo onto all bands:
    DO band=1, n_bands
      DO j=1, n_points
        saos(j, 1, band) = saos_bb(j, 1) * spectralsea(band)
        saos(j, 2, band) = saos_bb(j, 2) * spectralsea(band)
      END DO
    END DO
  ELSE
    ! copy the broadband albedo onto all bands:
    DO band=1, n_bands
      DO j=1, n_points
        saos(j, 1, band) = saos_bb(j, 1)
        saos(j, 2, band) = saos_bb(j, 2)
      END DO
    END DO
  END IF

ELSE IF (i_sea_alb_method == 3) THEN
  ! 3 - Jin et al. 2011, Optics Express
  ! (doi:10.1364/OE.19.026429)

  ! the Jin11 albedo parameterisation is a bit more complex/expensive
  ! than the others, so it is worth compressing the input fields onto 
  ! sea points only:
  ! (note sea_pts and sea_index are not all land, and not all sea-ice)
  ALLOCATE(cosz_s(sea_pts))
  ALLOCATE(ws_10m_s(sea_pts))
  ALLOCATE(chloro_s(sea_pts))
  ALLOCATE(saos_s(sea_pts, 2, n_bands))

  DO j = 1, sea_pts
    cosz_s(j) = cosz(sea_index(j))
    ws_10m_s(j) = ws_10m(sea_index(j))
    IF (l_sea_alb_var_chl) THEN
      chloro_s(j) = chloro(sea_index(j))
    END IF
  END DO
  IF (.NOT.l_sea_alb_var_chl) THEN
    ! set the chlorophyll content to the constant value:
    chloro_s = chloro_constant
  END IF

  IF (l_spec_sea_alb) THEN
    ! call the Jin 2011 albedo parameterisation, with spectral bands:
    CALL Jin11_osa(n_bands, spec_band_bot, spec_band_top,                      &
                   sea_pts, cosz_s, ws_10m_s, chloro_s, saos_s)
   
    ! push the output albedo on sea points back onto all points:
    saos = 0.0
    DO band = 1, n_bands
      DO j=1, sea_pts
        saos(sea_index(j),1,band) = saos_s(j,1,band)
        saos(sea_index(j),2,band) = saos_s(j,2,band)
      END DO
    END DO
  ELSE
    ! call the Jin 2011 albedo, for broadband:
    CALL Jin11_osa(1, spec_band_bb, spec_band_bb,                              &
                   sea_pts, cosz_s, ws_10m_s, chloro_s, saos_s)
   
    ! and push the output albedo on sea points back onto all points:
    saos = 0.0
    DO band = 1, n_bands
      DO j = 1, sea_pts
        saos(sea_index(j),1,band) = saos_s(j,1,1)
        saos(sea_index(j),2,band) = saos_s(j,2,1)
      END DO
    END DO
  END IF

  ! and finally deallocate arrays:
  DEALLOCATE(saos_s)
  DEALLOCATE(chloro_s)
  DEALLOCATE(ws_10m_s)
  DEALLOCATE(cosz_s)

ELSE
#if defined(UM_JULES)
  errcode = 1
  errmsg  = 'i_sea_alb_method does not have a valid setting'
  CALL ereport ('ftsa',errcode,errmsg)
#else
  WRITE(6,'(/,a)') 'ERROR: i_sea_alb_method has an invalid setting in ftsa'
  STOP
#endif
END IF


! Sea ice albedo

sa_sice (:,:)   = 0.0
alb_sice(:,:,:) = 0.0

IF (l_cice_alb) THEN

! Parametrisation from Los Alamos sea ice model (CICE vn4.1),
!         their 'default' option
  albice(1)=albicev  ! Direct=diffuse
  albice(2)=albicev
  albice(3)=albicei
  albice(4)=albicei
  albsnow(1)=albsnowv
  albsnow(2)=albsnowv
  albsnow(3)=albsnowi
  albsnow(4)=albsnowi

  DO band = 1, 4
    DO n = 1, nice_use
      DO l = 1, sice_pts_ncat(n)
        ll=sice_index_ncat(l,n)
        j=ssi_index(ll)

        ! Convert effective ice thickness to true ice thickness
        hice(j,n) = di_cat(j,n)                                         &
             - (kappai/kappai_snow) * (s_sea_cat(j,n)/rhosnow)

        ! Bare ice, thickness dependence
        fh = MIN(atan(hice(j,n)*4.0)/fhtan, 1.0)
        albo = albice(band)*fh + adifc*(1.0-fh)

        ! Bare ice, temperature dependence
        fT = MIN(tm-tstar_sice_cat(j,n)-1.0, 0.0)
        alb_sice(ll,n,band) = MAX(albo - dalb_mlt*fT, adifc)

        ! Snow, temperature dependence
        IF (s_sea_cat(j,n) > 0.) THEN
          albs(band) = albsnow(band)
          albs(band) = albs(band) - dalb_mlts(band)*fT
          ! Note rhosnow in next line is required to convert s_sea_cat
          !        from kg/m2 to m
          area_snow = s_sea_cat(j,n)                                    &
               / (s_sea_cat(j,n) + snowpatch*rhosnow)
        ELSE
          area_snow = 0.0
        ENDIF

        ! Combine snow and ice albedo
        alb_sice(ll,n,band) =                                           &
                  alb_sice(ll,n,band) * (1.0 - area_snow)  +            &
                  albs(band)*area_snow

        ! Mean sea ice albedo
        sa_sice(j, band) = sa_sice(j, band) +                           &
          alb_sice(ll,n,band)*sice_frac_ncat(ll,n) / aice(j)
      END DO
    END DO
  END DO

ELSE

! Use old UM schemes
! Derive 3 constants from the basic quantities (supplied in the
! UM namelist RUN_Radiation) for sea-ice albedo calculations:
  tcice = tm - dtice
  ice1 = (alpham-alphac) / dtice
  ice2 = alpham - tm*ice1

  DO n = 1, nice_use
    DO l = 1, sice_pts_ncat(n)
      ll=sice_index_ncat(l,n)
      j=ssi_index(ll)

      IF (l_ssice_albedo) THEN   ! Which albedo scheme?
        ! Effect of snow on sea-ice albedo

        IF (s_sea_cat(j,n) >  0.0) THEN   ! snow on sea ice
                                          ! Cox et al., Clim Dyn,1999
          IF (tstar_sice_cat(j,n) >  tcice) THEN
            snow_albedo=ice2+ice1*tstar_sice_cat(j,n)
          ELSE
            snow_albedo=alphac
          END IF

          alb_sice(ll,n,1)=alphab                                         &
          +(snow_albedo-alphab)*(1.0-EXP(-maskd*s_sea_cat(j,n)))

        ELSE           ! no snow so bare ice only

          IF(l_sice_meltponds) THEN

            ! bare ice, temperature dep. (Ebert &Curry,1993)
            IF (tstar_sice_cat(j,n) > (tm - dt_bare))THEN
              alb_sice(ll,n,1) = alphab+(tstar_sice_cat(j,n)-tm+dt_bare)  &
                                 *dalb_bare_wet
            ELSE      ! surface is dry
              alb_sice(ll,n,1)=alphab
            END IF     ! end melt ponds

          ELSE         ! just use bare ice albedo

            alb_sice(ll,n,1)=alphab

          END IF       ! l_sice_meltponds

          IF(l_sice_scattering) THEN
            ! Semtner modification dry albedo for internal
            ! scattering (Semnter, 1976)
            alb_sice(ll,n,1)=alb_sice(ll,n,1)+beta*(1.0-alb_sice(ll,n,1)) &
                                     *pen_rad_frac
          END IF       ! l_sice_scattering

        END IF         !  any snow on ice

      ELSE             ! default to operational NWP scheme 3.5.1:

        IF ( tstar_sice_cat(j,n)  <   tcice ) THEN
          alb_sice(ll,n,1) = alphac
        ELSE
          alb_sice(ll,n,1) = ice1 * tstar_sice_cat(j,n) + ice2
        END IF

      END IF          ! Which albedo scheme?

      IF (alb_sice(ll,n,1) <  0.0) THEN
        alb_sice(ll,n,1)=alpham
      END IF

!     For this scheme all band albedos are the same:
      DO band = 2, 4
        alb_sice(ll,n,band) = alb_sice(ll,n,1)
      END DO

    END DO  ! l
  END DO    ! n

! Calculate sea ice mean albedo
  DO band = 1, 4
    DO n = 1, nice_use
      DO l = 1, sice_pts_ncat(n)
        ll=sice_index_ncat(l,n)
        j=ssi_index(ll)
        sa_sice(j,band) = sa_sice(j,band) +                             &
               alb_sice(ll,n,band)*sice_frac_ncat(ll,n) / aice(j)
      END DO
    END DO
  END DO

END IF


IF (lhook) CALL dr_hook('FTSA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ftsa
