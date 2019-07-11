! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for picking up BL, and other, options

! Description:
!   Permissible settings for BL options.
!   This module replaces blopt8a.h for settings required by land surface.
!   It also includes variables from bl_options_mod, bl_diags_mod 
!   and swapable_field_mod

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE blopt8a

    IMPLICIT NONE

#include "blopt8a.h"

!---------------------------------------------------------------
! Duplicating later BL options from bl_options_mod
! for use in the standalone JULES code.
!
! NOTE THAT DEFAULT VALUES SHOULD BE SUPPLIED HERE
!

  INTEGER :: ISrfExCnvGust = 0
!                      ! Switch to include the effect of convective
!                      ! downdraughts on surface exchange
!                      ! OFF (=0) => not used: only boundary-layer
!                      !   gustiness is considered (oiginal version)
!                      ! IP_SrfExWithCnv (=1) the impact of gustiness
!                      !   due to boundary layer eddies is reduced
!                      !   relative to the above, but eddies driven
!                      !   by convective downdraughts are included
  INTEGER :: Fric_heating = 0
!                      ! Switch to apply heating source from turbulence 
!                      ! dissipation
!                      ! OFF (=0) => not used
!                      ! ON  (=1) => used
!

  REAL :: Max_Stress_Grad = 0.05
!                      ! Maximum implied stress gradient across the 
!                      ! boundary layer, used to limit the explicit 
!                      ! stress applied in non-local scalings (m/s2)
!---------------------------------------------------------------
! From bl_diags_mod: only need dTfric which appears in imp_solver
      TYPE Strnewbldiag
        Logical :: L_dTfric
        Real, Pointer :: dTfric(:, :, :)
!                    188      Heating increment from turbulence dissipation
      END TYPE Strnewbldiag
!----------------------------------------------------
! From swapable_field_mod
  TYPE swapable_field_pointer_type
    INTEGER :: field_type
    INTEGER :: levels
    INTEGER :: rows
    LOGICAL :: vector
    REAL, POINTER :: field(:, :, :)
    REAL, POINTER :: field_2d(:, :)
  END TYPE swapable_field_pointer_type


END MODULE blopt8a
