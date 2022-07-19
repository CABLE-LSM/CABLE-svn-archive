MODULE cable_rad_driv_mod

CONTAINS

SUBROUTINE cable_rad_driver( EffSurfRefl_dif, EffSurfRefl_beam,                &
                       mp, nrb, timestep_len, Clai_thresh, Ccoszen_tols,       &
                       CGauss_w, Cpi, Cpi180, Ctfrz, z0surf_min,               &
                       veg_mask, sunlit_mask, sunlit_veg_mask,                 &
                       jls_standalone, jls_radiation,  SurfaceType,            &
                       LAI_pft_cbl, HGT_pft_cbl, SnowDepth, SnowODepth,        &
                       SnowFlag_3L, SnowDensity, SoilTemp, SnowTemp, SnowAge,  &
                       AlbSoil ,coszen, VegXfang, VegTaul, VegRefl,            &
                       HeightAboveSnow, reducedLAIdue2snow )

!subrs:
USE cbl_albedo_mod,             ONLY: albedo
USE cbl_init_radiation_module,  ONLY: init_radiation

IMPLICIT NONE

! Albedos req'd by JULES - Effective Surface Relectance as seen by atmosphere
REAL, INTENT(OUT) :: EffSurfRefl_dif(mp,nrb)  ! Refl to Diffuse component of rad
                                              ! formerly rad%reffdf
REAL, INTENT(OUT) :: EffSurfRefl_beam(mp,nrb) ! Refl to Beam component of rad
                                              ! formerly rad%reffbm
!model dimensions
INTEGER, INTENT(IN) :: mp         ! total number of "tiles"
INTEGER, INTENT(IN) :: nrb        ! # radiation bands[ 1=VIS,2=NIR,3=LW(legacy, not used)

!---IN: JULES  timestep length in seconds
INTEGER, INTENT(IN) :: timestep_len

!constants
!-------------------------------------------------------------------------------
REAL, INTENT(IN) :: Ccoszen_tols      ! threshold cosine of sun's zenith angle,
                                      ! below which considered SUNLIT
REAL, INTENT(IN) :: Cgauss_w(nrb)     ! Gaussian integration weights
REAL, INTENT(IN) :: Clai_thresh       ! The minimum LAI below which a "cell" is
                                      ! considred  NOT vegetated
REAL, INTENT(IN) :: Cpi               ! PI
REAL, INTENT(IN) :: Cpi180            ! PI in radians
REAL, INTENT(IN) :: Ctfrz             ! freezing temp. of water
REAL, INTENT(IN) :: z0surf_min        ! the minimum roughness of bare soil

LOGICAL, INTENT(IN) :: jls_standalone ! local runtime switch for JULES(/UM) run
LOGICAL, INTENT(IN) :: jls_radiation  ! local runtime switch for radiation path
!-------------------------------------------------------------------------------

!masks
!-------------------------------------------------------------------------------
LOGICAL, INTENT(IN) :: veg_mask(:)        !  vegetated (uses min LAI)
LOGICAL, INTENT(IN) :: sunlit_mask(:)     !  sunlit (uses zenith angle)
LOGICAL, INTENT(IN) :: sunlit_veg_mask(:) !  BOTH sunlit AND  vegetated
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. remapped to "mp"
!-------------------------------------------------------------------------------
REAL, INTENT(IN) :: LAI_pft_cbl(mp)             !LAI -  "limited" and remapped
REAL, INTENT(IN) :: HGT_pft_cbl(mp)             !canopy height -  "limited" and remapped
!-------------------------------------------------------------------------------

REAL, INTENT(IN):: HeightAboveSnow(mp)         ! Height of Canopy above snow (rough%hruff)
                                    ! compute from z0surf_min, HGT_pft_cbl,
                                    ! SnowDepth, SnowDensity
REAL, INTENT(IN) :: reducedLAIdue2snow(mp)      ! Reduced LAI given snow coverage
REAL, INTENT(IN) :: coszen(mp)                  ! cosine zenith angle  (met%coszen)

REAL,INTENT(IN) :: AlbSoil(mp, nrb)              ! soil%AlbSoil

!Prognostics
!-------------------------------------------------------------------------------
REAL,INTENT(IN) :: SnowDepth(mp)               ! Total Snow depth - water eqivalent -
                                    ! packed from snow_surft (ssnow%snowd)
REAL :: SnowODepth(mp)              ! Total Snow depth before any update
!this timestep (ssnow%Osnowd)
REAL,INTENT(IN) :: SnowDensity(mp)             ! Total Snow density (assumes 1 layer
!describes snow cover) (ssnow%ssdnn)
REAL,INTENT(IN) :: SoilTemp(mp)                ! Soil Temperature of top layer (soil%tgg)
REAL,INTENT(IN) :: SnowTemp(mp)                ! Snow Temperature of top layer (soil%tgg)
REAL,INTENT(IN) :: SnowAge(mp)                 ! Snow age (assumes 1 layer describes snow
!cover) (ssnow%snage)
INTEGER,INTENT(IN):: SnowFlag_3L(mp)           ! Flag to treat snow as 3 layer - if enough
! snow present. Updated depending on total depth (ssnow%isflag)
!-------------------------------------------------------------------------------

!Vegetation parameters
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: SurfaceType(mp)
REAL :: VegXfang(mp)                !leaf angle PARAMETER (veg%xfang)
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
!-------------------------------------------------------------------------------

!local:
REAL, ALLOCATABLE :: ExtCoeff_beam(:)           ! rad%extkb,
REAL, ALLOCATABLE :: ExtCoeff_dif(:)            ! rad%extkd
REAL, ALLOCATABLE :: EffExtCoeff_beam(:, :)     ! rad%extkbm
REAL, ALLOCATABLE :: EffExtCoeff_dif(:, :)      ! rad%extkdm,
REAL, ALLOCATABLE :: CanopyTransmit_dif(:, :)   ! rad%cexpkdm
REAL, ALLOCATABLE :: CanopyTransmit_beam(:, :)  ! rad%cexpkbm
REAL, ALLOCATABLE :: CanopyRefl_dif(:, :)       ! rad%rhocdf
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)       ! rad%rhocbm
REAL, ALLOCATABLE :: RadFbeam(:, :)             ! rad%fbeam
REAL, ALLOCATABLE :: RadAlbedo(:, :)            ! rad%albedo
REAL, ALLOCATABLE :: AlbSnow(:, :)              ! ssnow%AlbSoilsn
REAL, ALLOCATABLE :: c1(:, :)
REAL, ALLOCATABLE :: rhoch(:, :)
REAL, ALLOCATABLE :: xk(:, :)
!NOT used on rad/albedo path. Need to fulfill arg list with dummy
INTEGER, ALLOCATABLE :: metDoY(:)               ! can pass DoY from current_time
REAL, ALLOCATABLE  :: SW_down(:,:)           !dummy

CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_driver"

LOGICAL :: cbl_standalone = .FALSE.

! alloc/zero each timestep
! metDoY, SW_down, RadFbeaam, RadAlbedo NOT used on rad/albedo path.
! Nevertheless,  need to fulfill later arg list(s) with dumies
CALL alloc_albedo_vars( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,         &
                        EffExtCoeff_dif, CanopyTransmit_dif,                   &
                        CanopyTransmit_beam, CanopyRefl_dif,CanopyRefl_beam,   &
                        c1, rhoch, xk, AlbSnow,                                &
                        RadFbeam, RadAlbedo, metDoY, SW_down, mp, nrb )

!Defines Extinction Coefficients to use in calculation of Canopy
!Reflectance/Transmitance.
CALL init_radiation( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,            &
                     EffExtCoeff_dif, RadFbeam, c1, rhoch, xk,                 &
                     mp,nrb, Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180, &
                     cbl_standalone, jls_standalone, jls_radiation, subr_name, &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
                     VegXfang, VegTaul, VegRefl, coszen, metDoY, SW_down,      &
                     reducedLAIdue2snow )

!Finally call albedo to get what we really need to fill contract with JULES
!Defines 4-"band" albedos [VIS/NIR bands. direct beam/diffuse components] from
!considering albedo of Ground (snow?) and Canopy Reflectance/Transmitance.

CALL Albedo( AlbSnow, AlbSoil,                                                 &
             mp, nrb, timestep_len,                                            &
             jls_radiation,                                                    &
             veg_mask, sunlit_mask, sunlit_veg_mask,                           &
             Ccoszen_tols, cgauss_w, Ctfrz,                                    &
             SurfaceType, VegRefl, VegTaul,                                    &
             coszen, reducedLAIdue2snow,                                       &
             SnowDepth, SnowODepth, SnowFlag_3L,                               &
             SnowDensity, SoilTemp, SnowTemp, SnowAge,                         &
             xk, c1, rhoch,                                                    &
             RadFbeam, RadAlbedo,                                              &
             ExtCoeff_dif, ExtCoeff_beam,                                      &
             EffExtCoeff_dif, EffExtCoeff_beam,                                &
             CanopyRefl_dif,CanopyRefl_beam,                                   &
             CanopyTransmit_dif, CanopyTransmit_beam,                          &
             EffSurfRefl_dif, EffSurfRefl_beam )

CALL flush_albedo( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,              &
                         EffExtCoeff_dif, CanopyTransmit_dif,                  &
                         CanopyTransmit_beam, CanopyRefl_dif, CanopyRefl_beam, &
                         RadFbeam, RadAlbedo, AlbSnow, c1, rhoch, xk )


RETURN

END SUBROUTINE cable_rad_driver

SUBROUTINE alloc_albedo_vars( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,   &
                        EffExtCoeff_dif, CanopyTransmit_dif,                   &
                        CanopyTransmit_beam, CanopyRefl_dif,CanopyRefl_beam,   &
                        c1, rhoch, xk, AlbSnow,                                &
                        RadFbeam, RadAlbedo, metDoY, SW_down, mp, nrb )
IMPLICIT NONE

INTEGER :: mp, nrb
!local to CABLE and can be flushed every timestep
REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, ALLOCATABLE :: AlbSnow(:,:)
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)
! these are dummies in JULES rad call but req'd to load arg lists
INTEGER, ALLOCATABLE :: metDoY(:)        ! can pass DoY from current_time
REAL, ALLOCATABLE :: SW_down(:,:)        ! dummy
REAL, ALLOCATABLE :: RadFbeam(:,:)
REAL, ALLOCATABLE :: RadAlbedo(:,:)

IF (.NOT. ALLOCATED(ExtCoeff_beam) )      ALLOCATE(ExtCoeff_beam(mp) )
IF (.NOT. ALLOCATED(ExtCoeff_dif) )       ALLOCATE(ExtCoeff_dif(mp) )
IF (.NOT. ALLOCATED(EffExtCoeff_beam) )   ALLOCATE(EffExtCoeff_beam(mp, nrb) )
IF (.NOT. ALLOCATED(EffExtCoeff_dif) )    ALLOCATE(EffExtCoeff_dif(mp, nrb) )
IF (.NOT. ALLOCATED(CanopyTransmit_dif))  ALLOCATE(CanopyTransmit_dif(mp, nrb))
IF (.NOT. ALLOCATED(CanopyTransmit_beam)) ALLOCATE(CanopyTransmit_beam(mp, nrb))
IF (.NOT. ALLOCATED(CanopyRefl_dif) )     ALLOCATE(CanopyRefl_dif(mp, nrb) )
IF (.NOT. ALLOCATED(CanopyRefl_beam) )    ALLOCATE(CanopyRefl_beam(mp, nrb) )
IF (.NOT. ALLOCATED(AlbSnow) )            ALLOCATE(AlbSnow(mp, nrb) )
IF (.NOT. ALLOCATED(c1) )                 ALLOCATE(c1(mp, nrb) )
IF (.NOT. ALLOCATED(rhoch) )              ALLOCATE(rhoch(mp, nrb) )
IF (.NOT. ALLOCATED(xk) )                 ALLOCATE(xk(mp, nrb) )
IF (.NOT. ALLOCATED(metDoY) )             ALLOCATE(metDoY(mp) )
IF (.NOT. ALLOCATED(SW_down) )            ALLOCATE(SW_down(mp,nrb) )
IF (.NOT. ALLOCATED(RadFbeam) )           ALLOCATE(RadFbeam(mp, nrb) )
IF (.NOT. ALLOCATED(RadAlbedo) )          ALLOCATE(RadAlbedo(mp, nrb) )

ExtCoeff_beam(:) = 0.0; ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0; EffExtCoeff_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 0.0; CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0; CanopyRefl_beam(:,:) = 0.0
AlbSnow(:,:) = 0.0
rhoch(:,:) = 0.0; xk(:,:) = 0.0; c1(:,:) = 0.0
RadFbeam(:,:) = 0.0; RadAlbedo(:,:) = 0.0
SW_down(:,:) = 0.0; metDoY(:) = 0.0  !can pass DoY from current_time%

RETURN

END SUBROUTINE alloc_albedo_vars


SUBROUTINE flush_albedo( ExtCoeff_beam, ExtCoeff_dif, EffExtCoeff_beam,        &
                         EffExtCoeff_dif, CanopyTransmit_dif,                  &
                         CanopyTransmit_beam, CanopyRefl_dif, CanopyRefl_beam, &
                         RadFbeam, RadAlbedo, AlbSnow, c1, rhoch, xk )
IMPLICIT NONE

!these local to CABLE and can be flushed every timestep
REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, ALLOCATABLE :: RadFbeam(:,:)
REAL, ALLOCATABLE :: RadAlbedo(:,:)
REAL, ALLOCATABLE :: AlbSnow(:,:)
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)

IF ( ALLOCATED (ExtCoeff_beam)       ) DEALLOCATE (ExtCoeff_beam )
IF ( ALLOCATED (ExtCoeff_dif)        ) DEALLOCATE (ExtCoeff_dif )
IF ( ALLOCATED (EffExtCoeff_beam)    ) DEALLOCATE (EffExtCoeff_beam )
IF ( ALLOCATED (EffExtCoeff_dif)     ) DEALLOCATE (EffExtCoeff_dif )
IF ( ALLOCATED (CanopyTransmit_dif)  ) DEALLOCATE (CanopyTransmit_dif )
IF ( ALLOCATED (CanopyTransmit_beam) ) DEALLOCATE (CanopyTransmit_beam )
IF ( ALLOCATED (CanopyRefl_dif)      ) DEALLOCATE (CanopyRefl_dif )
IF ( ALLOCATED (CanopyRefl_beam)     ) DEALLOCATE (CanopyRefl_beam )
IF ( ALLOCATED (RadFbeam)            ) DEALLOCATE (RadFbeam )
IF ( ALLOCATED (RadAlbedo)           ) DEALLOCATE (RadAlbedo )
IF ( ALLOCATED (AlbSnow)             ) DEALLOCATE (AlbSnow )
IF ( ALLOCATED (c1)                  ) DEALLOCATE (c1 )
IF ( ALLOCATED (rhoch)               ) DEALLOCATE (rhoch )
IF ( ALLOCATED (xk)                  ) DEALLOCATE (xk )

END SUBROUTINE flush_albedo

!==============================================================================


END MODULE cable_rad_driv_mod

