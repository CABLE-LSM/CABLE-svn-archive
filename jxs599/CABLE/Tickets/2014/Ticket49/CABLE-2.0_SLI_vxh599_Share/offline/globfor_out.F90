SUBROUTINE GLOBFOR_OUT(mp, pop, casapool,  veg, rad, cleaf_max, npp_ann, gpp_ann, stemnpp_ann, &
     leafnpp_ann)

  USE cable_common_module,  ONLY: cable_user
  USE POP_Types,            Only: POP_TYPE
  USE casavariable,         ONLY: casa_pool
  USE cable_def_types_mod,  ONLY: veg_parameter_type, radiation_type


  IMPLICIT NONE

  INTEGER,                  INTENT(IN) :: mp
  TYPE(casa_pool),          INTENT(IN) :: casapool
  TYPE(POP_TYPE),           INTENT(IN) :: POP
  TYPE(veg_parameter_type), INTENT(IN) :: veg  ! vegetation parameters
  TYPE(radiation_type),     INTENT(IN) :: rad
  REAL, DIMENSION(1:mp),    INTENT(IN) :: cleaf_max, npp_ann, gpp_ann, stemnpp_ann, leafnpp_ann

  INTEGER :: i,x, k_output
  REAL    :: summe
  LOGICAL, SAVE :: first_GF = .true.

  if (first_GF) then
     do x = 1,mp
        write(2000,"(2f8.2)") rad%latitude(x), rad%longitude(x)
     enddo
     first_GF = .false.
  endif

  k_output = 42
  DO x = 1, SIZE(pop%pop_grid(:))
     if (pop%pop_grid(x)%patch(k_output)%Layer(1)%density.gt.1e-9) then
        write(10000+x,"(i6, 13e16.6, i6)")  pop%pop_grid(x)%patch(k_output)%age(1), &
             log10(pop%pop_grid(x)%patch(k_output)%Layer(1)%density*10000), &
             log10(pop%pop_grid(x)%patch(k_output)%Layer(1)%biomass/0.49/0.7/pop%pop_grid(x)%patch(k_output)%Layer(1)%density), & ! log(total biomass)
             log10(pop%pop_grid(x)%patch(k_output)%Layer(1)%biomass/0.49/pop%pop_grid(x)%patch(k_output)%Layer(1)%density), & ! log(stem biomass)
             log10(cleaf_max(x)/1000/0.49/pop%pop_grid(x)%patch(k_output)%Layer(1)%density), &
             cleaf_max(x), casapool%cplant(x,2), casapool%cplant(x,3), npp_ann(x), gpp_ann(x), stemnpp_ann(x), &
             stemnpp_ann(x)/npp_ann(x), leafnpp_ann(x), leafnpp_ann(x)/npp_ann(x), veg%iveg(x)
     endif
  END DO



END SUBROUTINE GLOBFOR_OUT

