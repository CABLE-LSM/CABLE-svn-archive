module radiation_albedo_mod

contains

subroutine allocate_rad_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
coszen,        &
!VegXfang, VegTaul, VegRefl, 
c1, rhoch, &
RadFbeam, xk, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam &
)
implicit none
integer :: mp
integer :: nrb 

REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)

REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)  
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)  

REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)  
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)  

REAL, ALLOCATABLE :: reducedLAIdue2snow(:)
REAL, ALLOCATABLE :: coszen(:)
!REAL, ALLOCATABLE :: VegXfang(:)
!REAL, ALLOCATABLE :: VegTaul(:,:)
!REAL, ALLOCATABLE :: VegRefl(:,:)
REAL, ALLOCATABLE :: metDoY(:)
REAL, ALLOCATABLE :: SW_down(:,:)  
REAL, ALLOCATABLE :: RadFbeam(:,:)  
REAL, ALLOCATABLE :: RadAlbedo(:,:)  
REAL, ALLOCATABLE :: AlbSnow(:,:)  
REAL, ALLOCATABLE :: AlbSoil(:,:)  

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)     ! extinct. coef.for beam rad. and black leaves


!IF(.NOT. ALLOCATED( ) ) ALLOCATE( () )
IF(.NOT. ALLOCATED( reducedLAIdue2snow ) ) ALLOCATE( reducedLAIdue2snow(mp) )
IF(.NOT. ALLOCATED( coszen ) )             ALLOCATE( coszen(mp) )
!IF(.NOT. ALLOCATED( VegXfang ) )           ALLOCATE( VegXfang(mp) )
IF(.NOT. ALLOCATED( metDoY) )              ALLOCATE( metDoY(mp) )
IF(.NOT. ALLOCATED( ExtCoeff_beam) )       ALLOCATE( ExtCoeff_beam(mp) )
IF(.NOT. ALLOCATED( ExtCoeff_dif) )        ALLOCATE( ExtCoeff_dif(mp) )

!IF(.NOT. ALLOCATED( ) ) ALLOCATE( (mp,nrb) )
!IF(.NOT. ALLOCATED( VegTaul ) )            ALLOCATE( VegTaul(mp,nrb) )
!IF(.NOT. ALLOCATED( VegRefl ) )            ALLOCATE( VegRefl(mp,nrb) )
IF(.NOT. ALLOCATED( EffExtCoeff_beam) )    ALLOCATE( EffExtCoeff_beam(mp,nrb) )
IF(.NOT. ALLOCATED( EffExtCoeff_dif) )     ALLOCATE( EffExtCoeff_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyRefl_dif) )      ALLOCATE( CanopyRefl_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyRefl_beam) )     ALLOCATE( CanopyRefl_beam(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyTransmit_dif) )  ALLOCATE( CanopyTransmit_dif(mp,nrb) )
IF(.NOT. ALLOCATED( CanopyTransmit_beam) ) ALLOCATE( CanopyTransmit_beam(mp,nrb) )
IF(.NOT. ALLOCATED( EffSurfRefl_dif) )     ALLOCATE( EffSurfRefl_dif(mp,nrb) )
IF(.NOT. ALLOCATED( EffSurfRefl_beam) )    ALLOCATE( EffSurfRefl_beam(mp,nrb) )
IF(.NOT. ALLOCATED( c1 ) )                 ALLOCATE( c1(mp,nrb) )
IF(.NOT. ALLOCATED( xk ) )                 ALLOCATE( xk(mp,nrb) )
IF(.NOT. ALLOCATED( rhoch ) )              ALLOCATE( rhoch(mp,nrb) )
IF(.NOT. ALLOCATED( SW_down) )             ALLOCATE( SW_down(mp,nrb) )
IF(.NOT. ALLOCATED( RadFbeam) )            ALLOCATE( RadFbeam(mp,nrb) )
IF(.NOT. ALLOCATED( RadAlbedo) )           ALLOCATE( RadAlbedo(mp,nrb) )
IF(.NOT. ALLOCATED( AlbSnow) )             ALLOCATE( AlbSnow(mp,nrb) )
IF(.NOT. ALLOCATED( AlbSoil) )             ALLOCATE( AlbSoil(mp,nrb) )

End subroutine allocate_rad_albedo

subroutine deallocate_rad_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
coszen,        &
!VegXfang, VegTaul, VegRefl, 
c1, rhoch, &
RadFbeam, xk, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam &
)
implicit none
integer :: mp
integer :: nrb 
REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)

REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)  
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)  

REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)  
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)  

REAL, ALLOCATABLE :: reducedLAIdue2snow(:)
REAL, ALLOCATABLE :: coszen(:)
!REAL, ALLOCATABLE :: VegXfang(:)
!REAL, ALLOCATABLE :: VegTaul(:,:)
!REAL, ALLOCATABLE :: VegRefl(:,:)
REAL, ALLOCATABLE :: metDoY(:)
REAL, ALLOCATABLE :: SW_down(:,:)  
REAL, ALLOCATABLE :: RadFbeam(:,:)  
REAL, ALLOCATABLE :: RadAlbedo(:,:)  
REAL, ALLOCATABLE :: AlbSnow(:,:)  
REAL, ALLOCATABLE :: AlbSoil(:,:)  

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)     ! extinct. coef.for beam rad. and black leaves



IF( ALLOCATED( reducedLAIdue2snow ) ) DEALLOCATE( reducedLAIdue2snow )
IF( ALLOCATED( coszen ) )             DEALLOCATE( coszen )
!IF( ALLOCATED( VegXfang ) )           DEALLOCATE( VegXfang )
IF( ALLOCATED( metDoY) )              DEALLOCATE( metDoY )
IF( ALLOCATED( ExtCoeff_beam) )       DEALLOCATE( ExtCoeff_beam )
IF( ALLOCATED( ExtCoeff_dif) )        DEALLOCATE( ExtCoeff_dif )
!IF( ALLOCATED( VegTaul ) )            DEALLOCATE( VegTaul )
!IF( ALLOCATED( VegRefl ) )            DEALLOCATE( VegRefl )
IF( ALLOCATED( EffExtCoeff_beam) )    DEALLOCATE( EffExtCoeff_beam )
IF( ALLOCATED( EffExtCoeff_dif) )     DEALLOCATE( EffExtCoeff_dif )
IF( ALLOCATED( CanopyRefl_dif) )      DEALLOCATE( CanopyRefl_dif )
IF( ALLOCATED( CanopyRefl_beam) )     DEALLOCATE( CanopyRefl_beam )
IF( ALLOCATED( CanopyTransmit_dif) )  DEALLOCATE( CanopyTransmit_dif )
IF( ALLOCATED( CanopyTransmit_beam) ) DEALLOCATE( CanopyTransmit_beam )
IF( ALLOCATED( EffSurfRefl_dif) )     DEALLOCATE( EffSurfRefl_dif )
IF( ALLOCATED( EffSurfRefl_beam) )    DEALLOCATE( EffSurfRefl_beam )
IF( ALLOCATED( c1 ) )                 DEALLOCATE( c1 )
IF( ALLOCATED( xk ) )                 DEALLOCATE( xk )
IF( ALLOCATED( rhoch ) )              DEALLOCATE( rhoch )
IF( ALLOCATED( SW_down) )             DEALLOCATE( SW_down )
IF( ALLOCATED( RadFbeam) )            DEALLOCATE( RadFbeam )
IF( ALLOCATED( RadAlbedo) )           DEALLOCATE( RadAlbedo )
IF( ALLOCATED( AlbSnow) )             DEALLOCATE( AlbSnow )
IF( ALLOCATED( AlbSoil) )             DEALLOCATE( AlbSoil )

End subroutine deallocate_rad_albedo

End module radiation_albedo_mod
