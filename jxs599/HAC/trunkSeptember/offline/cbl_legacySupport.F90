module legacy_support_mod

contains

subroutine legacy_support( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
RadFbeam, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam, &
canopy, rad, ssnow )

USE cable_def_types_mod,        ONLY : canopy_type, &
                                    radiation_type, &
                                    soil_snow_type

implicit none
integer :: mp
integer :: nrb 
REAL :: ExtCoeff_beam(mp)
REAL :: ExtCoeff_dif(mp)
REAL :: EffExtCoeff_beam(mp,nrb)
REAL :: EffExtCoeff_dif(mp,nrb)

REAL :: CanopyTransmit_dif(mp,nrb)  
REAL :: CanopyTransmit_beam(mp,nrb)  
REAL :: CanopyRefl_dif(mp,nrb)  
REAL :: CanopyRefl_beam(mp,nrb)  

REAL :: EffSurfRefl_dif(mp,nrb)  
REAL :: EffSurfRefl_beam(mp,nrb)  

REAL :: reducedLAIdue2snow(mp)
REAL :: RadFbeam(mp,nrb)  
REAL :: RadAlbedo(mp,nrb)  
REAL :: AlbSnow(mp,nrb)  
REAL :: AlbSoil(mp,nrb)  
!REAL, ALLOCATABLE :: coszen(mp:)
!REAL, ALLOCATABLE :: metDoY(mp:)
!REAL, ALLOCATABLE :: SW_down(mp,:)  

! CABLE model variables
TYPE (canopy_type),    INTENT(INOUT) :: canopy
TYPE (radiation_type), INTENT(INOUT) :: rad
TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
!TYPE (roughness_type), INTENT(INOUT) :: rough

canopy%vlaiw  = reducedLAIdue2snow 
rad%extkb     = ExtCoeff_beam 
rad%extkd     = ExtCoeff_dif 
rad%extkbm    = EffExtCoeff_beam 
rad%extkdm    = EffExtCoeff_dif 
rad%cexpkdm   = CanopyRefl_dif 
rad%cexpkbm   = CanopyRefl_beam 
rad%rhocdf     = CanopyTransmit_dif 
rad%rhocbm     = CanopyTransmit_beam 
rad%reffdf    = EffSurfRefl_dif 
rad%reffbm    = EffSurfRefl_beam 
rad%fbeam     = RadFbeam 
rad%Albedo    = RadAlbedo 
ssnow%Albsoilsn  = AlbSnow 
!soil%AlbSoil    = AlbSoil 
!met%coszen    = coszen 
!met%doy       = metDoY 

End subroutine legacy_support

End module legacy_support_mod
