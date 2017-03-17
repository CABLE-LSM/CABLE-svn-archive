SUBROUTINE read_casa_dump(  ncfile, casamet, casaflux,phen, climate, ncall, kend, allATonce )
      USE netcdf
      USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
      USE casadimension,         ONLY : mplant,mdyear
      USE casavariable,          ONLY : casa_met, casa_flux
      USE phenvariable
#     ifndef UM_BUILD
      USE cable_diag_module,     ONLY : get_var_ncr2, &
                                        get_var_ncr3, stderr_nc
#     endif
      IMPLICIT NONE

      TYPE (casa_flux), INTENT(INOUT) :: casaflux
      TYPE (casa_met), INTENT(inout)  :: casamet
      TYPE (phen_variable),         INTENT(INOUT) :: phen
      TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
      INTEGER, INTENT(in)             :: kend, ncall
      CHARACTER(len=*), INTENT(in)    :: ncfile
      LOGICAL, INTENT(in)             :: allATonce

      !netcdf IDs/ names
      INTEGER, PARAMETER :: num_vars=14
      INTEGER, PARAMETER :: num_dims=3
      INTEGER, SAVE                        :: ncrid  ! netcdf file ID
      INTEGER , DIMENSION(num_vars)        :: varrID ! (1) tvair, (2) pmb

      !vars
      CHARACTER(len=*), DIMENSION(num_vars), PARAMETER :: &
            var_name =  (/  "lat          ", &
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
                            "mtemp        ", &
                            "Ndep         " /)

      REAL     , DIMENSION(mp)        :: lat, lon
      REAL(r_2), DIMENSION(mp)        :: tairk,  cgpp, mtemp, Ndep
      REAL(r_2), DIMENSION(mp,ms)     :: tsoil, moist
      REAL(r_2), DIMENSION(mp,mplant) :: crmplant
      REAL(r_2), DIMENSION(mp)        :: phenphase, phendoyphase1, &
           phendoyphase2,  phendoyphase3,  phendoyphase4
      INTEGER :: ncok,  idoy

#     ifndef UM_BUILD

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
            CALL get_var_ncr2(ncrid, var_name(13), mtemp   , idoy )
            CALL get_var_ncr2(ncrid, var_name(14), Ndep   , idoy )

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
            casamet%mtempspin(:,idoy) = mtemp
            casaflux%Nmindep = Ndep
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
         CALL get_var_ncr2(ncrid, var_name(13), mtemp   , ncall )
         CALL get_var_ncr2(ncrid, var_name(14), Ndep   , ncall )

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
         climate%mtemp_max = mtemp
         casaflux%Nmindep = Ndep

      ENDIF

      IF ( allATonce .OR. ncall .EQ. kend ) THEN
         ncok = NF90_CLOSE(ncrid)
         IF (ncok /= nf90_noerr ) CALL stderr_nc(ncok,'closing ', ncfile)
      ENDIF
#     endif

   END SUBROUTINE read_casa_dump


