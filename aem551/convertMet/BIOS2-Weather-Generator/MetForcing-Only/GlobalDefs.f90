MODULE GlobalDefs
USE TypeDef
! jtk561, to define some globally appplicalble stuff
real(sp) :: delT=60.0*60.0 ! subdiurnal timestep (in s)
integer(i4b) :: np=841*681     ! number of spatial elements
integer(i4b) :: ntime=24    ! number of subdiurnal steps in 24 hr
character(LEN=100) :: base_dir='/g/data1/fj2/AWAP/'
real :: cellsize=0.05
integer :: ncols=841
integer :: nrows=681
real :: yllcorner=-44.025,xllcorner=111.975 
END MODULE GlobalDefs
