MODULE soilin_pars_mod_cbl

USE grid_constants_mod_cbl, ONLY : nsoil_max   ! # of soil types [9]

IMPLICIT NONE

PUBLIC :: soilin_type
PUBLIC :: soilin

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Soil parameters I/O:
TYPE soilin_type
  REAL ::                                                                     &
       silt(nsoil_max),                                                        &
       clay(nsoil_max),                                                        &
       sand(nsoil_max),                                                        &
       swilt(nsoil_max),                                                       &
       sfc(nsoil_max),                                                         &
       ssat(nsoil_max),                                                        &
       bch(nsoil_max),                                                         &
       hyds(nsoil_max),                                                        &
       sucs(nsoil_max),                                                        &
       rhosoil(nsoil_max),                                                     &
       css(nsoil_max)
END TYPE soilin_type

!Instantiate types
TYPE(soilin_type) :: soilin !read from namelist

END MODULE soilin_pars_mod_cbl

