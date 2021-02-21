MODULE cable_fprint_type_mod

IMPLICIT NONE

TYPE fprint_type

  REAL, DIMENSION(:), POINTER ::                                              &
       fprint1,  & 
       fprint2,  & 
       fprint3,  & 
       fprint4,  & 
       fprint5,  & 
       fprint6,  & 
       fprint7,  & 
       fprint8,  & 
       fprint9

END TYPE fprint_type

!Instantiation:
TYPE ( fprint_type ) :: fprintx 

END MODULE cable_fprint_type_mod
