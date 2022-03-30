MODULE cable_soil_params_mod

USE grid_constants_mod_cbl, ONLY: nsoil_max   ! # of soil types [9]
IMPLICIT NONE 

TYPE csoilin_type

  CHARACTER(LEN=70) ::  desc(nsoil_max) 
   REAL, DIMENSION(nsoil_max) ::                                        &
      silt,    & !
      clay,    & !
      sand,    & !
      swilt,   & !
      sfc,     & !
      ssat,    & !
      bch,     & !
      hyds,    & !
      sucs,    & !
      rhosoil, & !
      css,     & !
      c3         !

END TYPE csoilin_type

TYPE(csoilin_type), SAVE  :: csoilin

TYPE soilin_type

  CHARACTER(LEN=70) ::  desc(nsoil_max)

   REAL, DIMENSION(nsoil_max) ::                                        &
      silt,    & !
      clay,    & !
      sand,    & !
      swilt,   & !
      sfc,     & !
      ssat,    & !
      bch,     & !
      hyds,    & !
      sucs,    & !
      rhosoil, & !
      css,     & !
      c3         !

END TYPE soilin_type

CHARACTER(LEN=70), DIMENSION(nsoil_max) ::                                 &
   soil_desc     ! decriptns of soil type 

TYPE(soilin_type), SAVE  :: soilin

CONTAINS


subroutine cable_soil_params()

! Gets parameter values for each vegetation type and soil type.
USE cable_def_types_mod, ONLY : mstype

implicit none

integer :: ERROR
integer :: namelist_unit
integer :: j
CHARACTER(LEN=*), parameter :: iomessage='something wrong with your soil params file' 
CHARACTER(LEN=*), parameter :: nml_dir='./' 
CHARACTER(LEN=*), PARAMETER :: routinename='SOILN_CABLE'

NAMELIST / cable_soilparm / soilin
!hard-wired #  of soil types, promote to nml
mstype = nsoil_max 
 
!SOIL parameters: description and corresponding variable name in code. 
!SOIL parameters are assigned as TYPE soilin% but later used as soil%

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
write (6,*) "Reading CABLE_SOILPARM namelist..."

OPEN( namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'cable_soilparm.nml'),          &
      STATUS='old', POSITION='rewind', ACTION='read', IOSTAT  = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error opening  CABLE_SOILPARM namelist..."

READ(namelist_unit, NML = cable_soilparm, IOSTAT = ERROR )
IF ( ERROR /= 0 ) write (6,*) "Error reading  CABLE_SOILPARM namelist..."
                 
CLOSE(namelist_unit, IOSTAT = ERROR)

soil_desc = soilin%desc

!SOIL: Coarse sand/Loamy sand                                                
! =========================================================
    csoilin%silt( 1) =        0.080000
    csoilin%clay( 1) =        0.090000
    csoilin%sand( 1) =        0.830000
   csoilin%swilt( 1) =        0.072000
     csoilin%sfc( 1) =        0.143000
    csoilin%ssat( 1) =        0.398000
     csoilin%bch( 1) =        4.200000
    csoilin%hyds( 1) =        0.000166
    csoilin%sucs( 1) =       -0.106000
 csoilin%rhosoil( 1) =     1600.000000
     csoilin%css( 1) =      850.000000
 
 !SOIL: Medium clay loam/silty clay loam/silt loam                            
 !=========================================================
    csoilin%silt( 2) =        0.330000
    csoilin%clay( 2) =        0.300000
    csoilin%sand( 2) =        0.370000
   csoilin%swilt( 2) =        0.216000
     csoilin%sfc( 2) =        0.301000
    csoilin%ssat( 2) =        0.479000
     csoilin%bch( 2) =        7.100000
    csoilin%hyds( 2) =        0.000004
    csoilin%sucs( 2) =       -0.591000
 csoilin%rhosoil( 2) =     1600.000000
     csoilin%css( 2) =      850.000000
 
 !SOIL: Fine clay                                                             
 !=========================================================
    csoilin%silt( 3) =        0.170000
    csoilin%clay( 3) =        0.670000
    csoilin%sand( 3) =        0.160000
   csoilin%swilt( 3) =        0.286000
     csoilin%sfc( 3) =        0.367000
    csoilin%ssat( 3) =        0.482000
     csoilin%bch( 3) =       11.400000
    csoilin%hyds( 3) =        0.000001
    csoilin%sucs( 3) =       -0.405000
 csoilin%rhosoil( 3) =     1381.000000
     csoilin%css( 3) =      850.000000
 
 !SOIL: Coarse-medium sandy loam/loam                                         
 !=========================================================
    csoilin%silt( 4) =        0.200000
    csoilin%clay( 4) =        0.200000
    csoilin%sand( 4) =        0.600000
   csoilin%swilt( 4) =        0.135000
     csoilin%sfc( 4) =        0.218000
    csoilin%ssat( 4) =        0.443000
     csoilin%bch( 4) =        5.150000
    csoilin%hyds( 4) =        0.000021
    csoilin%sucs( 4) =       -0.348000
 csoilin%rhosoil( 4) =     1373.000000
     csoilin%css( 4) =      850.000000
 
 !SOIL: Coarse-fine sandy clay                                                
 !=========================================================
    csoilin%silt( 5) =        0.060000
    csoilin%clay( 5) =        0.420000
    csoilin%sand( 5) =        0.520000
   csoilin%swilt( 5) =        0.219000
     csoilin%sfc( 5) =        0.310000
    csoilin%ssat( 5) =        0.426000
     csoilin%bch( 5) =       10.400000
    csoilin%hyds( 5) =        0.000002
    csoilin%sucs( 5) =       -0.153000
 csoilin%rhosoil( 5) =     1476.000000
     csoilin%css( 5) =      850.000000
 
 !SOIL: Medium-fine silty clay                                                
 !=========================================================
    csoilin%silt( 6) =        0.250000
    csoilin%clay( 6) =        0.480000
    csoilin%sand( 6) =        0.270000
   csoilin%swilt( 6) =        0.283000
     csoilin%sfc( 6) =        0.370000
    csoilin%ssat( 6) =        0.482000
     csoilin%bch( 6) =       10.400000
    csoilin%hyds( 6) =        0.000001
    csoilin%sucs( 6) =       -0.490000
 csoilin%rhosoil( 6) =     1521.000000
     csoilin%css( 6) =      850.000000
 
 !SOIL: Coarse-medium-fine sandy clay loam                                    
 !=========================================================
    csoilin%silt( 7) =        0.150000
    csoilin%clay( 7) =        0.270000
    csoilin%sand( 7) =        0.580000
   csoilin%swilt( 7) =        0.175000
     csoilin%sfc( 7) =        0.255000
    csoilin%ssat( 7) =        0.420000
     csoilin%bch( 7) =        7.120000
    csoilin%hyds( 7) =        0.000006
    csoilin%sucs( 7) =       -0.299000
 csoilin%rhosoil( 7) =     1373.000000
     csoilin%css( 7) =      850.000000
 
 !SOIL: Organic peat                                                          
 !=========================================================
    csoilin%silt( 8) =        0.700000
    csoilin%clay( 8) =        0.170000
    csoilin%sand( 8) =        0.130000
   csoilin%swilt( 8) =        0.395000
     csoilin%sfc( 8) =        0.450000
    csoilin%ssat( 8) =        0.451000
     csoilin%bch( 8) =        5.830000
    csoilin%hyds( 8) =        0.000800
    csoilin%sucs( 8) =       -0.356000
 csoilin%rhosoil( 8) =     1537.000000
     csoilin%css( 8) =     1920.000000
 
 !SOIL: Permanent ice                                                         
 !=========================================================
    csoilin%silt( 9) =        0.330000
    csoilin%clay( 9) =        0.300000
    csoilin%sand( 9) =        0.370000
   csoilin%swilt( 9) =        0.216000
     csoilin%sfc( 9) =        0.301000
    csoilin%ssat( 9) =        0.479000
     csoilin%bch( 9) =        7.100000
    csoilin%hyds( 9) =        0.000001
    csoilin%sucs( 9) =       -0.153000
  csoilin%rhosoil( 9) =      910.000000
     csoilin%css( 9) =     2100.000000    




do j=1, 9

!IF( csoilin%desc(j)    /=  soilin%desc(j)     ) write(6,*) "jh:j, csoilin%desc(j)   , soilin%desc(j)    )  ",j, csoilin%desc(j)   , soilin%desc(j)    
IF( csoilin%silt(j)    /=  soilin%silt(j)     ) write(6,*) "jh:j, csoilin%silt(j)   , soilin%silt(j)    )  ",j, csoilin%silt(j)   , soilin%silt(j)    
IF( csoilin%clay(j)    /=  soilin%clay(j)     ) write(6,*) "jh:j, csoilin%clay(j)   , soilin%clay(j)    )  ",j, csoilin%clay(j)   , soilin%clay(j)    
IF( csoilin%sand(j)    /=  soilin%sand(j)     ) write(6,*) "jh:j, csoilin%sand(j)   , soilin%sand(j)    )  ",j, csoilin%sand(j)   , soilin%sand(j)    
IF( csoilin%swilt(j)   /=  soilin%swilt(j)    ) write(6,*) "jh:j, csoilin%swilt(j)  , soilin%swilt(j)   )  ",j, csoilin%swilt(j)  , soilin%swilt(j)    
IF( csoilin%sfc(j)     /=  soilin%sfc(j)      ) write(6,*) "jh:j, csoilin%sfc(j)    , soilin%sfc(j)     )  ",j, csoilin%sfc(j)    , soilin%sfc(j)     
IF( csoilin%ssat(j)    /=  soilin%ssat(j)     ) write(6,*) "jh:j, csoilin%ssat(j)   , soilin%ssat(j)    )  ",j, csoilin%ssat(j)   , soilin%ssat(j)      
IF( csoilin%bch(j)     /=  soilin%bch(j)      ) write(6,*) "jh:j, csoilin%bch(j)    , soilin%bch(j)     )  ",j, csoilin%bch(j)    , soilin%bch(j)        
IF( csoilin%hyds(j)    /=  soilin%hyds(j)     ) write(6,*) "jh:j, csoilin%hyds(j)   , soilin%hyds(j)    )  ",j, csoilin%hyds(j)   , soilin%hyds(j)      
IF( csoilin%sucs(j)    /=  soilin%sucs(j)     ) write(6,*) "jh:j, csoilin%sucs(j)   , soilin%sucs(j)    )  ",j, csoilin%sucs(j)   , soilin%sucs(j)      
IF( csoilin%rhosoil(j) /=  soilin%rhosoil(j)  ) write(6,*) "jh:j, csoilin%rhosoil(j), soilin%rhosoil(j) )  ",j, csoilin%rhosoil(j), soilin%rhosoil(j) 
IF( csoilin%css(j)     /=  soilin%css(j)      ) write(6,*) "jh:j, csoilin%css(j)    , soilin%css(j)     )  ",j, csoilin%css(j)    , soilin%css(j)        

end do


End subroutine cable_soil_params

END MODULE cable_soil_params_mod

