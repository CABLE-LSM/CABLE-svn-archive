MODULE read_dump_casa_mod

IMPLICIT NONE

#define  UM_BUILD YES

#ifndef UM_BUILD

CONTAINS

SUBROUTINE read_casa_dump(  ncfile, casamet, casaflux,phen, climate, ncall, kend, allATonce )
      USE netcdf
      USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
      USE casadimension,         ONLY : mplant,mdyear
      USE casavariable,          ONLY : casa_met, casa_flux
      USE phenvariable
      USE cable_common_module,  ONLY:  CABLE_USER
      USE casa_ncdf_module,     ONLY : get_var_ncr2, &
                                        get_var_ncr3, stderr_nc
      IMPLICIT NONE

      TYPE (casa_flux), INTENT(INOUT) :: casaflux
      TYPE (casa_met), INTENT(inout)  :: casamet
      TYPE (phen_variable),         INTENT(INOUT) :: phen
      TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
      INTEGER, INTENT(in)             :: kend, ncall
      CHARACTER(len=*), INTENT(in)    :: ncfile
      LOGICAL, INTENT(in)             :: allATonce

      !netcdf IDs/ names
      INTEGER            :: num_vars
      INTEGER, PARAMETER :: num_dims=3
      INTEGER, SAVE                        :: ncrid  ! netcdf file ID
      !INTEGER , DIMENSION(num_vars)        :: varrID ! (1) tvair, (2) pmb


      !vars
      CHARACTER, DIMENSION(:), POINTER :: var_name*15

      REAL     , DIMENSION(mp)        :: lat, lon
      REAL(r_2), DIMENSION(mp)        :: tairk,  cgpp, mtemp, Ndep, Pdep
      REAL(r_2), DIMENSION(mp,ms)     :: tsoil, moist
      REAL(r_2), DIMENSION(mp,mplant) :: crmplant
      REAL(r_2), DIMENSION(mp)        :: phenphase, phendoyphase1, &
           phendoyphase2,  phendoyphase3,  phendoyphase4
      INTEGER :: ncok,  idoy

      num_vars=14

      !Add extra mtemp variable when running with climate
      IF(cable_user%CALL_climate) THEN
          num_vars=num_vars+1
      ENDIF

      allocate(var_name(num_vars))

      !Variable names
      var_name =  (/"lat          ", &
                    "lon          ", &
                    "casamet_tairk", &
                    "tsoil        ", &
                    "moist        ", &
                    "cgpp         ", &
                    "crmplant     ", &
                    "phenphase    ", &
                    "phendoyphase1", &
                    "phendoyphase2", &
                    "phendoyphase3", &
                    "phendoyphase4", &
                    "Ndep         ", &
                    "Pdep         "/)

      !Add extra mtemp variable when running with climate
      IF (cable_user%CALL_climate) THEN
        var_name(num_vars)="mtemp"
      ENDIF

 IF ( allATonce .OR. ncall .EQ. 1 ) THEN
         ncok = NF90_OPEN(TRIM(ncfile), nf90_nowrite, ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'re-opening ', ncfile)
      ENDIF
      IF ( allATonce ) THEN
         DO idoy=1,mdyear

            CALL get_var_ncr2(ncrid, var_name(3), tairk   , idoy )
            CALL get_var_ncr3(ncrid, var_name(4), tsoil   , idoy ,ms)
            CALL get_var_ncr3(ncrid, var_name(5), moist   , idoy ,ms)
            CALL get_var_ncr2(ncrid, var_name(6), cgpp    , idoy )
            CALL get_var_ncr3(ncrid, var_name(7), crmplant, idoy ,mplant)
            CALL get_var_ncr2(ncrid, var_name(8), phenphase, idoy)
            CALL get_var_ncr2(ncrid, var_name(9), phendoyphase1, idoy)
            CALL get_var_ncr2(ncrid, var_name(10), phendoyphase2, idoy)
            CALL get_var_ncr2(ncrid, var_name(11), phendoyphase3, idoy)
            CALL get_var_ncr2(ncrid, var_name(12), phendoyphase4, idoy)
            CALL get_var_ncr2(ncrid, var_name(13), Ndep   , idoy )
            CALL get_var_ncr2(ncrid, var_name(14), Pdep   , idoy )
            !amu561 this need to be in if-block
            IF (cable_user%CALL_climate ) THEN
               CALL get_var_ncr2(ncrid, var_name(15), mtemp   , idoy )
            ENDIF

            casamet%Tairkspin(:,idoy) = tairk
            casamet%cgppspin (:,idoy) = cgpp
            casamet%crmplantspin_1(:,idoy) = crmplant(:,1)
            casamet%crmplantspin_2(:,idoy) = crmplant(:,2)
            casamet%crmplantspin_3(:,idoy) = crmplant(:,3)
            casamet%Tsoilspin_1(:,idoy)    = tsoil(:,1)
            casamet%Tsoilspin_2(:,idoy)    = tsoil(:,2)
            casamet%Tsoilspin_3(:,idoy)    = tsoil(:,3)
            casamet%Tsoilspin_4(:,idoy)    = tsoil(:,4)
            casamet%Tsoilspin_5(:,idoy)    = tsoil(:,5)
            casamet%Tsoilspin_6(:,idoy)    = tsoil(:,6)
            casamet%moistspin_1(:,idoy)    = moist(:,1)
            casamet%moistspin_2(:,idoy)    = moist(:,2)
            casamet%moistspin_3(:,idoy)    = moist(:,3)
            casamet%moistspin_4(:,idoy)    = moist(:,4)
            casamet%moistspin_5(:,idoy)    = moist(:,5)
            casamet%moistspin_6(:,idoy)    = moist(:,6)
            phen%phasespin(:,idoy) = int(phenphase)
            phen%doyphasespin_1(:,idoy) = int(phendoyphase1)
            phen%doyphasespin_2(:,idoy) = int(phendoyphase2)
            phen%doyphasespin_3(:,idoy) = int(phendoyphase3)
            phen%doyphasespin_4(:,idoy) = int(phendoyphase4)
            casaflux%Nmindep = Ndep
            casaflux%Pdep = Pdep
            !amu561 this need to be in if-block
            IF(cable_user%CALL_climate) THEN
               casamet%mtempspin(:,idoy) = mtemp
            ENDIF
         END DO
      ELSE

         CALL get_var_ncr2(ncrid, var_name(3), tairk   ,ncall )
         CALL get_var_ncr3(ncrid, var_name(4), tsoil   ,ncall , ms)
         CALL get_var_ncr3(ncrid, var_name(5), moist   ,ncall , ms)
         CALL get_var_ncr2(ncrid, var_name(6), cgpp    ,ncall )
         CALL get_var_ncr3(ncrid, var_name(7), crmplant,ncall , mplant)
         CALL get_var_ncr2(ncrid, var_name(8), phenphase    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(9), phendoyphase1    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(10), phendoyphase2    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(11), phendoyphase3    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(12), phendoyphase4    ,ncall )
         CALL get_var_ncr2(ncrid, var_name(13), Ndep   , ncall )
         CALL get_var_ncr2(ncrid, var_name(14), Pdep   , ncall )
         IF(cable_user%CALL_climate) THEN
            CALL get_var_ncr2(ncrid, var_name(15), mtemp   , ncall )
         ENDIF

         casamet%tairk     = tairk
         casamet%tsoil     = tsoil
         casamet%moist     = moist
         casaflux%cgpp     = cgpp
         casaflux%crmplant = crmplant
         phen%phase = int(phenphase)
         phen%doyphase(:,1) = int(phendoyphase1)
         phen%doyphase(:,2) = int(phendoyphase2)
         phen%doyphase(:,3) = int(phendoyphase3)
         phen%doyphase(:,4) = int(phendoyphase4)
         casaflux%Nmindep = Ndep
         casaflux%Pdep = Pdep
         IF (cable_user%CALL_climate) THEN
            climate%mtemp_max = mtemp
         ENDIF

      ENDIF

      IF ( allATonce .OR. ncall .EQ. kend ) THEN
         ncok = NF90_CLOSE(ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'closing ', ncfile)
      ENDIF

      deallocate(var_name)

   END SUBROUTINE read_casa_dump
#endif

END MODULE read_dump_casa_mod
