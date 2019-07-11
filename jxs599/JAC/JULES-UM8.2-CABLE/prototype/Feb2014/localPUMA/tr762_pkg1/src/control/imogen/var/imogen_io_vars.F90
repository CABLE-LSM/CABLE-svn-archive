#if !defined(UM_JULES)
MODULE imogen_io_vars

  IMPLICIT NONE

  INTEGER ::                                                      &
    YR_EMISS(300)
             ! Years in which CO2 emissions are prescribed
             ! (only first NYR_EMISS array components are used)
             ! For this code, years must be sequential

  REAL ::                                                         &
    C_EMISS(300)
             ! Values of CO2 emissions read in (up to NYR_EMISS)

END MODULE imogen_io_vars
#endif
