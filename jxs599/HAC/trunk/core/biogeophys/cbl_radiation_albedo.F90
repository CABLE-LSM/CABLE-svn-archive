module radiation_albedo_mod

public

REAL, ALLOCATABLE, SAVE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE, SAVE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE, SAVE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE, SAVE :: EffExtCoeff_dif(:,:)

REAL, ALLOCATABLE, SAVE :: CanopyTransmit_dif(:,:)  
REAL, ALLOCATABLE, SAVE :: CanopyTransmit_beam(:,:)  
REAL, ALLOCATABLE, SAVE :: CanopyRefl_dif(:,:)  
REAL, ALLOCATABLE, SAVE :: CanopyRefl_beam(:,:)  

REAL, ALLOCATABLE, SAVE :: EffSurfRefl_dif(:,:)  
REAL, ALLOCATABLE, SAVE :: EffSurfRefl_beam(:,:)  

REAL, ALLOCATABLE, SAVE :: reducedLAIdue2snow(:)
REAL, ALLOCATABLE, SAVE :: coszen(:)
REAL, ALLOCATABLE, SAVE :: VegXfang(:)
REAL, ALLOCATABLE, SAVE :: VegTaul(:,:)
REAL, ALLOCATABLE, SAVE :: VegRefl(:,:)
REAL, ALLOCATABLE, SAVE :: metDoY(:)
REAL, ALLOCATABLE, SAVE :: SW_down(:,:)  
REAL, ALLOCATABLE, SAVE :: RadFbeam(:,:)  
REAL, ALLOCATABLE, SAVE :: RadAlbedo(:,:)  
REAL, ALLOCATABLE, SAVE :: AlbSnow(:,:)  
REAL, ALLOCATABLE, SAVE :: AlbSoil(:,:)  

!co-efficients used throughout init_radiation ` called from _albedo as well
REAL, ALLOCATABLE, SAVE :: c1(:,:)
REAL, ALLOCATABLE, SAVE :: rhoch(:,:)
REAL, ALLOCATABLE, SAVE :: xk(:,:)     ! extinct. coef.for beam rad. and black leaves

contains

subroutine allocate_rad_albedo( mp, nrb )
implicit none
integer :: mp
integer :: nrb 

!IF(.NOT. ALLOCATED( ) ) ALLOCATE( () )
IF(.NOT. ALLOCATED( reducedLAIdue2snow ) ) ALLOCATE( reducedLAIdue2snow(mp) )
IF(.NOT. ALLOCATED( coszen ) )      ALLOCATE( coszen(mp) )
IF(.NOT. ALLOCATED( VegXfang ) )    ALLOCATE( VegXfang(mp) )
IF(.NOT. ALLOCATED( metDoY) )     ALLOCATE( metDoY(mp) )
IF(.NOT. ALLOCATED( ExtCoeff_beam) )   ALLOCATE( ExtCoeff_beam(mp) )
IF(.NOT. ALLOCATED( ExtCoeff_dif) )    ALLOCATE( ExtCoeff_dif(mp) )


!IF(.NOT. ALLOCATED( ) ) ALLOCATE( (mp,nrb) )
IF(.NOT. ALLOCATED( VegTaul ) )     ALLOCATE( VegTaul(mp,nrb) )
IF(.NOT. ALLOCATED( VegRefl ) )     ALLOCATE( VegRefl(mp,nrb) )
IF(.NOT. ALLOCATED( EffExtCoeff_beam) )   ALLOCATE( EffExtCoeff_beam(mp,nrb) )
IF(.NOT. ALLOCATED( EffExtCoeff_dif) )    ALLOCATE( EffExtCoeff_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyRefl_dif) )  ALLOCATE( CanopyRefl_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyRefl_beam) ) ALLOCATE( CanopyRefl_beam(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyTransmit_dif) )  ALLOCATE( CanopyTransmit_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyTransmit_beam) ) ALLOCATE( CanopyTransmit_beam(mp,nrb) )
IF(.NOT. ALLOCATED( EffSurfRefl_dif) )  ALLOCATE( EffSurfRefl_dif(mp,nrb) )
IF(.NOT. ALLOCATED( EffSurfRefl_beam) ) ALLOCATE( EffSurfRefl_beam(mp,nrb) )
IF(.NOT. ALLOCATED( c1 ) )             ALLOCATE( c1(mp,nrb) )
IF(.NOT. ALLOCATED( xk ) )             ALLOCATE( xk(mp,nrb) )
IF(.NOT. ALLOCATED( rhoch ) )          ALLOCATE( rhoch(mp,nrb) )
IF(.NOT. ALLOCATED( SW_down) )          ALLOCATE( SW_down(mp,nrb) )
IF(.NOT. ALLOCATED( RadFbeam) )          ALLOCATE( RadFbeam(mp,nrb) )
IF(.NOT. ALLOCATED( RadAlbedo) )          ALLOCATE( RadAlbedo(mp,nrb) )
IF(.NOT. ALLOCATED( AlbSnow) )          ALLOCATE( AlbSnow(mp,nrb) )
IF(.NOT. ALLOCATED( AlbSoil) )          ALLOCATE( AlbSoil(mp,nrb) )

End subroutine allocate_rad_albedo

End module radiation_albedo_mod
