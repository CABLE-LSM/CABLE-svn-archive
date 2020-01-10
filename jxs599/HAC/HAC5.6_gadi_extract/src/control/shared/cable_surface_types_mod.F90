! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cable_surface_types_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains CABLE surface type information and a namelist for setting them
!   The CABLE equivalent of jules_surface_types_mod
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE 
!-----------------------------------------------------------------------------

USE max_dimensions,   ONLY:                                                   &
    elev_tile_max,                                                            &
    ntype_max

USE missing_data_mod, ONLY: imdi

USE jules_surface_types_mod, ONLY: nnpft, ncpft, npft, nnvg, ntype, urban,    &
                                   ice, soil

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
!  nnpft,                                                                      &
                ! Number of natural pfts
                !   Derived from namelist inputs as npft - ncpft
!  ncpft = 0,                                                                  &
                ! Number of crop pfts
!  npft  = imdi,                                                               &
                ! Number of plant functional types
!  nnvg  = imdi,                                                               &
                ! Number of non-vegetation surface types
!  ntype,                                                                      &
                ! Total number of surface types
!  urban = imdi,                                                               &
                ! Index of the urban surface type
  lakes  = imdi,                                                              &
                ! Index of the lake surface type
!  ice   = imdi 
                ! Index of the ice surface type
  barren= imdi
                ! Index of the barren surface type

!-----------------------------------------------------------------------------
! UM only: The following are only required for UM flexible tiles.
!-----------------------------------------------------------------------------
! The UM uses the following to identify the original PFTs in the UM output.

integer ::                                                                    &
   evergreen_needleleaf = imdi,                                               &
                  ! Index of surface type 'evergreen needleleaf'
   evergreen_broadleaf  = imdi,                                               &
                  ! Index of surface type 'evergreen broadleaf'
   deciduous_needleleaf = imdi,&
                  ! Index of surface type 'deciduous needleleaf'
   deciduous_broadleaf  = imdi,&
                  ! Index of surface type 'deciduous broadleaf'
   shrub            = imdi,                                                  &
                  ! Index of surface type 'shrub'
   c3_grassland=imdi,&
                  ! Index of surface type 'C3 grassland'
   c4_grassland=imdi,&
                  ! Index of surface type 'C4 grassland'
   tundra=imdi,&
                  ! Index of surface type 'tundra'
   c3_cropland=imdi,&
                  ! Index of surface type 'C3 cropland'
   c4_cropland=imdi,&
                  ! Index of surface type 'C4 cropland'
   wetland=imdi,&
                  ! Index of surface type 'wetland'
   empty1=imdi,&
                  ! Index of surface type 'empty 1'
   empty2=imdi,&
                  ! Index of surface type 'empty 2'

! Allows new tiles to be added to the UM without a code change.
    usr_type(ntype_max)  = imdi,                                              &
                  ! Index of user defined surface types
    tile_map_ids(ntype_max)    = imdi,                                        &
                  ! User specified tile map
    tile_map_pslevs(ntype_max) = imdi


!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / cable_surface_types/                                              &
    npft, ncpft, nnvg,                                                        &
    evergreen_needleleaf, evergreen_broadleaf, deciduous_needleleaf,          &
    deciduous_broadleaf, shrub, c3_grassland, c4_grassland,                   &
    tundra, c3_cropland, c4_cropland, wetland, empty1, empty2,                &
    barren, urban, lakes, ice

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CABLE_SURFACE_TYPES_MOD'

CONTAINS

SUBROUTINE check_cable_surface_types()

USE max_dimensions, ONLY: npft_max, ncpft_max, nnvg_max

USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks CABLE_SURFACE_TYPES namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE 
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: i ! Loop counter

INTEGER :: errorstatus

!-----------------------------------------------------------------------------
! Check that the given values are less than the fixed values for IO
!-----------------------------------------------------------------------------
IF ( npft > npft_max ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "npft > npft_max - increase npft_max and recompile")
END IF
IF ( ncpft > ncpft_max ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "ncpft > ncpft_max - increase ncpft_max and recompile")
END IF
IF ( nnvg > nnvg_max ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "nnvg > nnvg_max - increase nnvg_max and recompile")
END IF
IF ( ncpft > npft ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "ncpft > npft - total number of PFTs must be >= " //           &
               "number of crop PFTs")
END IF
! If these are true, then we automatically have:
!   ntype <= ntype_max (since ntype(_max) = npft(_max) + nnvg(_max))
!   nnpft <= nnpft_max (since nnpft(_max) = npft(_max) - ncpft(_max))

!-----------------------------------------------------------------------------
! Check values for the specific surface types are sensible
!-----------------------------------------------------------------------------
! PFT surface types must come before non-veg types, so if urban, lakes, ice or
! urban are given (i.e. > 0) then they must be > npft
! A soil type is required
IF ( urban > 0 .AND. ( urban <= npft .OR. urban > ntype ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "urban tile is given but is out of range")
END IF

IF ( lakes > 0 .AND. ( lakes <= npft .OR. lakes > ntype ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "lakes tile is given but is out of range")
END IF

IF ( ice > 0 .AND. ( ice <= npft .OR. ice > ntype ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "ice tile is given but is out of range")
END IF

END SUBROUTINE check_cable_surface_types


SUBROUTINE print_nlist_cable_surface_types()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

INTEGER :: i ! Loop counter

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------

CALL jules_print('cable_surface_types',                                       &
                 'Contents of namelist cable_surface_types')

WRITE(lineBuffer, *) '  npft = ', npft
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  nnvg = ', nnvg
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  evergreen_needleleaf = ', evergreen_needleleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  evergreen_broadleaf = ', evergreen_broadleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  deciduous_needleleaf = ', deciduous_broadleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  shrub = ', shrub
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c3_grassland = ', c3_grassland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c4_grassland = ', c4_grassland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  tundra = ', tundra
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c3_cropland = ', c3_cropland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c4_cropland = ', c4_cropland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  wetland = ', wetland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  empty1 = ', empty1
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  empty2 = ', empty2
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  barren = ', barren
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  urban = ', urban
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  lakes = ', lakes
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  ice = ', ice
CALL jules_print('cable_surface_types', lineBuffer)

END SUBROUTINE print_nlist_cable_surface_types

#if defined(UM_JULES) && !defined(LFRIC)
SUBROUTINE read_nml_cable_surface_types (unitnumber)

! Description:
!  Read the CABLE_SURFACE_TYPES namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_CABLE_SURFACE_TYPES'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 25 + (2 * ntype_max) + (2 * elev_tile_max)
    !CABLE_LSM:showed as a CABLE modification in 10.6. ONLY makes sense IF n_int was indeed increased
    !to 20 in earlier 8.5 and change was carried anyway
    !INTEGER, PARAMETER :: n_int = 20

!CABLE_LSM:
TYPE my_namelist
  SEQUENCE
  INTEGER :: npft
  INTEGER :: ncpft
  INTEGER :: nnvg
  INTEGER :: evergreen_needleleaf
  INTEGER :: evergreen_broadleaf
  INTEGER :: deciduous_needleleaf
  INTEGER :: deciduous_broadleaf
  INTEGER :: shrub
  INTEGER :: c3_grassland
  INTEGER :: c4_grassland
  INTEGER :: tundra
  INTEGER :: c3_cropland
  INTEGER :: c4_cropland
  INTEGER :: wetland
  INTEGER :: empty1
  INTEGER :: empty2
  INTEGER :: barren
  INTEGER :: urban
  INTEGER :: lakes
  INTEGER :: ice
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = cable_surface_types, IOSTAT = errorstatus,   &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist cable_surface_types",iomessage)

      !CABLE_LSM:
  my_nml % npft                 = npft
  my_nml % ncpft                = ncpft
  my_nml % nnvg                 = nnvg
  my_nml % evergreen_needleleaf = evergreen_needleleaf
  my_nml % evergreen_broadleaf  =  evergreen_broadleaf
  my_nml % deciduous_needleleaf = deciduous_needleleaf
  my_nml % deciduous_broadleaf  = deciduous_broadleaf
  my_nml % shrub                = shrub
  my_nml % c3_grassland         = c3_grassland
  my_nml % c4_grassland         = c4_grassland
  my_nml % tundra               = tundra
  my_nml % c3_cropland          = c3_cropland
  my_nml % c4_cropland          = c4_cropland
  my_nml % wetland              = wetland
  my_nml % empty1               = empty1
  my_nml % empty2               = empty2
  my_nml % barren               = barren
  my_nml % urban                = urban
  my_nml % lakes                = lakes
  my_nml % ice                  = ice
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  npft                 = my_nml % npft
  ncpft                = my_nml % ncpft
  nnvg                 = my_nml % nnvg
  evergreen_needleleaf = my_nml % evergreen_needleleaf
  evergreen_broadleaf  = my_nml % evergreen_broadleaf
  deciduous_needleleaf = my_nml % deciduous_needleleaf
  deciduous_broadleaf  = my_nml % deciduous_broadleaf
  shrub                = my_nml % shrub
  c3_grassland         = my_nml % c3_grassland
  c4_grassland         = my_nml % c4_grassland
  tundra               = my_nml % tundra
  c3_cropland          = my_nml % c3_cropland
  c4_cropland          = my_nml % c4_cropland
  wetland              = my_nml % wetland
  empty1               = my_nml % empty1
  empty2               = my_nml % empty2
  barren               = my_nml % barren
  urban                = my_nml % urban
  lakes                = my_nml % lakes
  ice                  = my_nml % ice

END IF

! JULES expects the soil field, equivalent to CABLE's barren, to be set
soil = barren

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_cable_surface_types
#endif

SUBROUTINE set_derived_variables_cable_surface_types()

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Derive ntype and nnpft from the namelist values
!-----------------------------------------------------------------------------
ntype = npft + nnvg
nnpft = npft - ncpft

RETURN

END SUBROUTINE set_derived_variables_cable_surface_types

END MODULE cable_surface_types_mod
