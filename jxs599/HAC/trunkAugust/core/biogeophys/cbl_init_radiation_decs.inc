integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
real :: Clai_thresh             !threshold LAI below which considered UN-vegetated
real :: Ccoszen_tols            !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
real :: Cpi                     !PI - from cable_math_constants originally
real :: Cpi180                  !PI in radians - from cable_math_constants originally
LOGICAL :: jls_standalone       !runtime switch defined in cable_*main routines signifying this is jules_standalone
LOGICAL :: jls_radiation        !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
!masks
logical :: veg_mask(mp)         !vegetated mask [formed by comparrisson of LAI CLAI_thresh ]
logical :: sunlit_mask(mp)      !sunlit mask [formed by comparrisson of coszen to coszen_tols i.e. is the sun up]
logical :: sunlit_veg_mask(mp)  !combined mask - BOTH sunlit and vegetated

REAL :: reducedLAIdue2snow(mp)         !Effective LAI given (potential sno coverage)
REAL :: coszen(mp)              ! cosine zenith angle of sun

REAL :: ExtCoeff_beam(mp)       !"raw" Extinction co-efficient for Direct Beam component of SW radiation
REAL :: ExtCoeff_dif(mp)        !"raw"Extinction co-efficient for Diffuse component of SW radiation
REAL :: EffExtCoeff_beam(mp,nrb)!Effective Extinction co-efficient for Direct Beam component of SW radiation
REAL :: EffExtCoeff_dif(mp,nrb) !Effective Extinction co-efficient for Diffuse component of SW radiation

integer :: metDoY(mp)           !Day of the Year [formerly met%doy]
REAL :: SW_down(mp,nrb)         !Downward SW radiation [formerly met%fsd]
REAL :: RadFbeam(mp,nrb)        !Beam Fraction of Downward SW radiation [formerly rad%fbeam]

!vaegetation parameters input via namelist
REAL :: VegXfang(mp)
REAL :: VegTaul(mp,nrb)
REAL :: VegRefl(mp,nrb)

!co-efficients used throughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)              ! extinct. coef.for beam rad. and black leaves

character(len=*) :: subr_name !where am i called from
