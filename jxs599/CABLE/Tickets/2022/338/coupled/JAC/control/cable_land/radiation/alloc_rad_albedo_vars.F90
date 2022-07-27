MODULE alloc_rad_albedo_vars_mod

IMPLICIT NONE

CONTAINS

! Allocate vars in mp format
SUBROUTINE alloc_local_vars( EffSurfRefl_dif, EffSurfRefl_beam,mp, nrb,        &
                             reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,     &
                             HeightAboveSnow, coszen, ExtCoeff_beam,           &
                             ExtCoeff_dif, EffExtCoeff_beam, EffExtCoeff_dif,  &
                             CanopyTransmit_dif, CanopyTransmit_beam,          &
                             CanopyRefl_dif,CanopyRefl_beam, c1, rhoch, xk,    &
                             AlbSnow, RadFbeam, RadAlbedo, metDoY, SW_down )
IMPLICIT NONE

INTEGER :: mp, nrb
REAL, INTENT(OUT), ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: LAI_pft_cbl(:)      ! Prescribed LAI 
REAL, INTENT(OUT), ALLOCATABLE :: HGT_pft_cbl(:)      ! Prescribed canopy height 
REAL, INTENT(OUT), ALLOCATABLE :: reducedLAIdue2snow(:) ! Eff. LAI given snow
REAL, INTENT(OUT), ALLOCATABLE :: HeightAboveSnow(:)  ! Canopy hgt above snow
REAL, INTENT(OUT), ALLOCATABLE :: coszen(:)
REAL, INTENT(OUT), ALLOCATABLE :: ExtCoeff_beam(:)
REAL, INTENT(OUT), ALLOCATABLE :: ExtCoeff_dif(:)
REAL, INTENT(OUT), ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: AlbSnow(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: c1(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: rhoch(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: xk(:,:)
! these are dummies in JULES rad call but req'd to load arg lists
REAL, INTENT(OUT), ALLOCATABLE :: SW_down(:,:)        ! dummy
REAL, INTENT(OUT), ALLOCATABLE :: RadFbeam(:,:)
REAL, INTENT(OUT), ALLOCATABLE :: RadAlbedo(:,:)
INTEGER, ALLOCATABLE :: metDoY(:)        ! can pass DoY from current_time

IF ( .NOT. ALLOCATED(reducedLAIdue2snow))  ALLOCATE(reducedLAIdue2snow(mp) )
IF ( .NOT. ALLOCATED(HeightAboveSnow) )    ALLOCATE(HeightAboveSnow(mp) )
IF ( .NOT. ALLOCATED(LAI_pft_cbl) )        ALLOCATE(LAI_pft_cbl(mp) )
IF ( .NOT. ALLOCATED(HGT_pft_cbl) )        ALLOCATE(HGT_pft_cbl(mp) )
IF ( .NOT. ALLOCATED(coszen) )             ALLOCATE(coszen(mp) )
IF ( .NOT. ALLOCATED(EffSurfRefl_dif) )    ALLOCATE(EffSurfRefl_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(EffSurfRefl_beam) )   ALLOCATE(EffSurfRefl_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(ExtCoeff_beam) )      ALLOCATE(ExtCoeff_beam(mp) )
IF ( .NOT. ALLOCATED(ExtCoeff_dif) )       ALLOCATE(ExtCoeff_dif(mp) )
IF ( .NOT. ALLOCATED(EffExtCoeff_beam) )   ALLOCATE(EffExtCoeff_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(EffExtCoeff_dif) )    ALLOCATE(EffExtCoeff_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(CanopyTransmit_dif))  ALLOCATE(CanopyTransmit_dif(mp, nrb))
IF ( .NOT. ALLOCATED(CanopyTransmit_beam)) ALLOCATE(CanopyTransmit_beam(mp,nrb))
IF ( .NOT. ALLOCATED(CanopyRefl_dif) )     ALLOCATE(CanopyRefl_dif(mp, nrb) )
IF ( .NOT. ALLOCATED(CanopyRefl_beam) )    ALLOCATE(CanopyRefl_beam(mp, nrb) )
IF ( .NOT. ALLOCATED(AlbSnow) )            ALLOCATE(AlbSnow(mp, nrb) )
IF ( .NOT. ALLOCATED(c1) )                 ALLOCATE(c1(mp, nrb) )
IF ( .NOT. ALLOCATED(rhoch) )              ALLOCATE(rhoch(mp, nrb) )
IF ( .NOT. ALLOCATED(xk) )                 ALLOCATE(xk(mp, nrb) )
IF ( .NOT. ALLOCATED(metDoY) )             ALLOCATE(metDoY(mp) )
IF ( .NOT. ALLOCATED(SW_down) )            ALLOCATE(SW_down(mp,nrb) )
IF ( .NOT. ALLOCATED(RadFbeam) )           ALLOCATE(RadFbeam(mp, nrb) )
IF ( .NOT. ALLOCATED(RadAlbedo) )          ALLOCATE(RadAlbedo(mp, nrb) )

EffSurfRefl_dif(:,:) = 0.0; EffSurfRefl_beam(:,:) = 0.0
LAI_pft_cbl(:) = 0.0; HGT_pft_cbl(:) = 0.0; coszen(:) = 0.0
reducedLAIdue2snow(:) = 0.0; HeightAboveSnow(:) = 0.0  
ExtCoeff_beam(:) = 0.0; ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0; EffExtCoeff_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 0.0; CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0; CanopyRefl_beam(:,:) = 0.0
AlbSnow(:,:) = 0.0
rhoch(:,:) = 0.0; xk(:,:) = 0.0; c1(:,:) = 0.0
RadFbeam(:,:) = 0.0; RadAlbedo(:,:) = 0.0
SW_down(:,:) = 0.0; metDoY(:) = 0  !can pass DoY from current_time%

RETURN

END SUBROUTINE alloc_local_vars

!flush memory
SUBROUTINE flush_local_vars( EffSurfRefl_dif, EffSurfRefl_beam, SnowDepth,     &
                             SnowDensity, SoilTemp, SnowAge, AlbSoil,          &
                             reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,     &
                             HeightAboveSnow, coszen, ExtCoeff_beam,           & 
                             ExtCoeff_dif, EffExtCoeff_beam, EffExtCoeff_dif,  &
                             CanopyTransmit_dif, CanopyTransmit_beam,          &
                             CanopyRefl_dif, CanopyRefl_beam, RadFbeam,        &
                             RadAlbedo, AlbSnow, c1, rhoch, xk, metDoY, SW_down) )

IMPLICIT NONE

REAL, INTENT(INOUT), ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: SnowDepth(:)
REAL, INTENT(INOUT), ALLOCATABLE :: SnowDensity(:)
REAL, INTENT(INOUT), ALLOCATABLE :: SoilTemp(:)
REAL, INTENT(INOUT), ALLOCATABLE :: SnowAge( :)
REAL, INTENT(INOUT), ALLOCATABLE :: AlbSoil(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: reducedLAIdue2snow(:)
REAL, INTENT(INOUT), ALLOCATABLE :: HeightAboveSnow(:)
REAL, INTENT(INOUT), ALLOCATABLE :: LAI_pft_cbl(:)
REAL, INTENT(INOUT), ALLOCATABLE :: HGT_pft_cbl(:)
REAL, INTENT(INOUT), ALLOCATABLE :: coszen(:)
!these local to CABLE and can be flushed every timestep
REAL, INTENT(INOUT), ALLOCATABLE :: ExtCoeff_beam(:)
REAL, INTENT(INOUT), ALLOCATABLE :: ExtCoeff_dif(:)
REAL, INTENT(INOUT), ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: EffExtCoeff_dif(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: CanopyTransmit_dif(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: CanopyTransmit_beam(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: CanopyRefl_dif(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: CanopyRefl_beam(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: RadFbeam(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: RadAlbedo(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: AlbSnow(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: c1(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: rhoch(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: xk(:,:)
REAL, INTENT(INOUT), ALLOCATABLE :: SW_down(:,:)        ! dummy
INTEGER, INTENT(INOUT),, ALLOCATABLE :: metDoY(:)  ! pass DoY from current_time

IF ( ALLOCATED(EffSurfRefl_dif)     ) DEALLOCATE ( EffSurfRefl_dif )
IF ( ALLOCATED(EffSurfRefl_beam)    ) DEALLOCATE ( EffSurfRefl_beam )
IF ( ALLOCATED(SnowDepth)           ) DEALLOCATE( SnowDepth )
IF ( ALLOCATED(SnowDensity)         ) DEALLOCATE( SnowDensity )
IF ( ALLOCATED(SoilTemp)            ) DEALLOCATE( SoilTemp )
IF ( ALLOCATED(SnowAge)             ) DEALLOCATE( SnowAge )
IF ( ALLOCATED(AlbSoil)             ) DEALLOCATE( AlbSoil )
IF ( ALLOCATED(reducedLAIdue2snow)  ) DEALLOCATE( reducedLAIdue2snow )
IF ( ALLOCATED(HeightAboveSnow)     ) DEALLOCATE( HeightAboveSnow )
IF ( ALLOCATED(LAI_pft_cbl)         ) DEALLOCATE( LAI_pft_cbl )
IF ( ALLOCATED(HGT_pft_cbl)         ) DEALLOCATE( HGT_pft_cbl )
IF ( ALLOCATED(coszen)              ) DEALLOCATE( coszen )
IF ( ALLOCATED (ExtCoeff_beam)      ) DEALLOCATE (ExtCoeff_beam )
IF ( ALLOCATED (ExtCoeff_dif)       ) DEALLOCATE (ExtCoeff_dif )
IF ( ALLOCATED (EffExtCoeff_beam)   ) DEALLOCATE (EffExtCoeff_beam )
IF ( ALLOCATED (EffExtCoeff_dif)    ) DEALLOCATE (EffExtCoeff_dif )
IF ( ALLOCATED (CanopyTransmit_dif) ) DEALLOCATE (CanopyTransmit_dif )
IF ( ALLOCATED (CanopyTransmit_beam)) DEALLOCATE (CanopyTransmit_beam )
IF ( ALLOCATED (CanopyRefl_dif)     ) DEALLOCATE (CanopyRefl_dif )
IF ( ALLOCATED (CanopyRefl_beam)    ) DEALLOCATE (CanopyRefl_beam )
IF ( ALLOCATED (RadFbeam)           ) DEALLOCATE (RadFbeam )
IF ( ALLOCATED (RadAlbedo)          ) DEALLOCATE (RadAlbedo )
IF ( ALLOCATED (AlbSnow)            ) DEALLOCATE (AlbSnow )
IF ( ALLOCATED (c1)                 ) DEALLOCATE (c1 )
IF ( ALLOCATED (rhoch)              ) DEALLOCATE (rhoch )
IF ( ALLOCATED (xk)                 ) DEALLOCATE (xk )
IF ( ALLOCATED (SW_down)            ) DEALLOCATE (SW_down )
IF ( ALLOCATED (MetDoY)             ) DEALLOCATE (MetDoY )

RETURN
END SUBROUTINE flush_local_vars

END MODULE alloc_rad_albedo_vars_mod
