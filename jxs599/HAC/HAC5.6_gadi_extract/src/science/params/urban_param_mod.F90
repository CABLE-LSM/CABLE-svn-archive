! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
! This file belongs in section: Land

MODULE urban_param

IMPLICIT NONE

! Parameters for the MacDonald et al. (2008) formulation of displacement height
! and effective roughness length for momentum
REAL, PARAMETER :: a       = 4.43
REAL, PARAMETER :: cdz     = 1.2     ! Drag coefficient
REAL, PARAMETER :: kappa2  = 0.16    ! Von Karman constant**2
REAL, PARAMETER :: z0m_mat = 0.05    ! Material roughness length for momentum

! Note: z0m_mat has a value of 0.005 in original UM version, but 0.05 in
! original JULES version. z0m_mat = 0.05 was used in inter-comparison for VL92

REAL, PARAMETER ::                                                            &
   emiss       = 1.0,           & ! Emissivity sky
   omega_day   = 7.272e-5         ! Angular velocity of the Earth
                                  ! wrt sun( s-1 )

! At the moment set urban materials here. Could be re-written to be read from
! a look up table for different fabrics. If this was done these would need to
! be made into arrays(land_points) to store the values here and the code
! re-written to enable this.

REAL, PARAMETER ::                                                            &
   ! Road material = asphalt
   diffus_road  = 0.38e-6,     & ! Road: Thermal diffusivity (m2 s-1)
   cap_road     = 1.94e6,      & ! Road: Volumetric heat capacity (J K-1 m-3)

   ! Wall material = brick
   diffus_wall = 0.61e-6,      & ! Wall: Thermal diffusivity (m2 s-1)
   cap_wall    = 1.37e6,       & ! Wall: Volumetric heat capacity (J K-1 m-3)

   ! Roof material = clay
   diffus_roof = 0.47e-6,      & ! Roof: Thermal diffusivity (m2 s-1)
   cap_roof    = 1.77e6,       & ! Roof: Volumetric heat capacity (J K-1 m-3)
   dz_roof_p   = 0.02            ! Physical depth of roof as opposed to

REAL, ALLOCATABLE, SAVE ::                                                    &
   hgt_gb(:),        &  ! Building height
   hwr_gb(:),        &  ! Height to width ratio
   wrr_gb(:),        &  ! Width ratio
   disp_gb(:),       &  ! Displacemnet height
   ztm_gb(:),        &  ! Roughness length
   albwl_gb(:),      &  ! Wall albedo
   albrd_gb(:),      &  ! Road albedo
   emisw_gb(:),      &  ! Wall emissivity
   emisr_gb(:)          ! Road emissivity

REAL ::                                                                       &
   anthrop_heat_scale = 1.0 ! Scales anthropogenic heat source of roof &
                            ! canyon from being equally distributed (= 1.0)
                            ! to all being released in the canyon (= 0.0).
                            ! Takes a value between 0.0 - 1.0

!-----------------------------------------------------------------------
! Set up a namelist to allow URBAN-2T parameters to be set
!
! In standalone, all parameters except anthrop_heat_scale are populated
! from the non-veg parameters in init_nonveg
!-----------------------------------------------------------------------
  NAMELIST  / jules_urban2t_param/                                            &
     anthrop_heat_scale

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='URBAN_PARAM'

CONTAINS
SUBROUTINE print_nlist_jules_urban2t_param()
USE jules_print_mgr, ONLY: jules_print
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('urban_param',                                               &
    'Contents of namelist jules_urban2t_param')

WRITE(lineBuffer,*)' anthrop_heat_scale = ',anthrop_heat_scale
CALL jules_print('urban_param',lineBuffer)

CALL jules_print('urban_param',                                               &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_urban2t_param

#if defined(UM_JULES) && !defined(LFRIC)
SUBROUTINE read_nml_jules_urban2t_param (unitnumber)

! Description:
!  Read the JULES_URBAN2T_PARAM namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_URBAN2T_PARAM'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_real = 1

TYPE my_namelist
  SEQUENCE
  REAL :: anthrop_heat_scale
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in = n_real)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_urban2t_param, IOSTAT = errorstatus,   &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_urban2t_param", iomessage)

  my_nml % anthrop_heat_scale = anthrop_heat_scale
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  anthrop_heat_scale = my_nml % anthrop_heat_scale
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_urban2t_param
#endif

END MODULE urban_param
