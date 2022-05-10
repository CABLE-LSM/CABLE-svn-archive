MODULE snow_aging_mod

  !This routine evaluates the effective age of the snow pack
  !It is called from cbl_soilsnow_main and the output SnowAge
  !is used in cbl_snow_albedo

PUBLIC snow_aging

CONTAINS

  SUBROUTINE snow_aging(SnowAge,mp,dels,SnowDepth,SnowODepth,SnowTemp, &
       SoilTemp,SnowFlag_3L, surface_type,soil_type)

    !may need to eliminate this USE data statement  
    USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ

    IMPLICIT NONE

    REAL, INTENT(IN OUT)   :: SnowAge(mp)
    INTEGER, INTENT(IN) :: mp
    REAL, INTENT(IN)    :: dels
    REAL, INTENT(IN)    :: SnowDepth(mp)
    REAL, INTENT(IN)    :: SnowODepth(mp)
    REAL, INTENT(IN)    :: SnowTemp(mp)
    REAL, INTENT(IN)    :: SoilTemp(mp)
    INTEGER, INTENT(IN) :: SnowFlag_3L(mp)
    INTEGER, INTENT(IN) :: surface_type(mp) 
    INTEGER, INTENT(IN) :: soil_type(mp) 

    !hard wired index to be eliminated
    INTEGER, PARAMETER :: perm_ice = 9
    !model parameter shared across subroutines -> cable_phys_constants
    REAL, PARAMETER :: snow_depth_thresh = 1.0

    !working variables - could be converted to scalars
    REAL, DIMENSION(mp) ::                             &
         ar1,     &  ! factor for crystal growth  (-ve)
         ar2,     &  ! factor for freezing of melt water
         ar3,     &  ! factor for accumulation of dirt
         dnsnow,  &  ! depth of new snow albedo
         dtau        !change in effective snow age

    INTEGER :: i        !looping variable
    REAL    :: tmp      !local effective surface temperauture

    !initialise at values if there is no snow
    ar1(:)=0.0
    ar2(:)=0.0
    ar3(:)=0.0
    DO i=1,mp
       IF (SnowDepth(i)>snow_depth_thresh) THEN

          !depth of new snow (in cm H20)
          dnsnow(i) = MIN (1.0, 0.1 * MAX( 0.0, SnowDepth(i) - SnowODepth(i) ) )

          ! Snow age depends on snow crystal growth, freezing of melt water,
          ! accumulation of dirt and amount of new snow.
          tmp = SnowFlag_3L(i) * SnowTemp(i) + ( 1 - SnowFlag_3L(i) ) * SoilTemp(i)
          tmp = MIN( tmp, CTFRZ )
          ar1(i) = 5000.0 * (1.0 / (CTFRZ-0.01) - 1.0 / tmp) ! crystal growth  (-ve)
          ar2(i) = 10.0 * ar1(i)                             ! freezing of melt water

          IF (soil_type(i) == perm_ice) THEN
             ! permanent ice case
             ar3(i) = 0.0000001
             
             !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
             !dnsnow = max (dnsnow(i), 0.5) !increase refreshing of snow in Antarctic
             dnsnow(i) = 1.0

          ELSE
             ! accumulation of dirt
             ar3(i) = 0.1

          END IF

          !update snow age
          dtau(i) = 1.0e-6 * (EXP( ar1(i) ) + EXP( ar2(i) ) + ar3(i) ) * dels
          SnowAge(i) = MAX(0.0,(SnowAge(i)+dtau(i))*(1.0-dnsnow(i)))

       END IF

    END DO

    RETURN

  END SUBROUTINE snow_aging

END MODULE snow_aging_mod
