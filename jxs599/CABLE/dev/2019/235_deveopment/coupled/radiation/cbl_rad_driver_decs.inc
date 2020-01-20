
!model dimensions
!-------------------------------------------------------------------------------
!JaC:todo:ultimatelty get this from JaC~
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
!-------------------------------------------------------------------------------

!constants
!-------------------------------------------------------------------------------
real :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)               !Gaussian integration weights
real :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
real :: Cpi                         !PI - describing the ratio of circumference to diameter
real :: Cpi180                      !PI in radians
real :: z0surf_min                  !the minimum roughness of bare soil
LOGICAL :: jls_standalone           !runtime switch defined in cable_*main routines signifying this is a JULES(/UM) run 
LOGICAL :: jls_radiation            !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
!-------------------------------------------------------------------------------

!masks
!-------------------------------------------------------------------------------
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_mask(mp)          ! this "mp" is sunlit (uses zenith angle)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. remapped to "mp"
!-------------------------------------------------------------------------------
integer:: surface_type(mp)          ! Integer index of Surface type (veg%iveg)
real :: LAI_pft_cbl(mp)             !LAI -  "limited" and remapped
real :: HGT_pft_cbl(mp)             !canopy height -  "limited" and remapped
!-------------------------------------------------------------------------------

real :: HeightAboveSnow(mp)         !Height of Canopy above snow (rough%hruff)
                                    !compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity
real :: reducedLAIdue2snow(mp)      ! Reduced LAI given snow coverage

!Forcing
!-------------------------------------------------------------------------------
REAL :: MetTk(mp)                   !Air Temperture at surface - atmospheric forcing (met%tk)
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
integer :: metDoY(mp)                  !Day of the Year - not always available (met%doy)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)
!-------------------------------------------------------------------------------

!Prognostics
!-------------------------------------------------------------------------------
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (ssnow%snowd)
REAL :: SnowODepth(mp)              !Total Snow depth before any update this timestep (ssnow%Osnowd)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (ssnow%ssdnn)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer (soil%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (ssnow%snage)
integer:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer  - if enough present. Updated depending on total depth (ssnow%isflag)
!-------------------------------------------------------------------------------
                                                                                           
! Albedos
!-------------------------------------------------------------------------------
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------
REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings [computed in albedo() ]
!-------------------------------------------------------------------------------
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
!-------------------------------------------------------------------------------

!Variables shared primarily between radiation and albedo and possibly elsewhere
!-------------------------------------------------------------------------------
!Extinction co-efficients compued in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-efficient for Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-efficient for Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%rhodf   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%rhobm)    
!-------------------------------------------------------------------------------

!Vegetation parameters
!-------------------------------------------------------------------------------
REAL :: VegXfang(mp)                !leaf angle PARAMETER (veg%xfang)
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
!-------------------------------------------------------------------------------