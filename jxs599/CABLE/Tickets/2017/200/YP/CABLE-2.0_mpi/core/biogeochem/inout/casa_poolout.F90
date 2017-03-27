SUBROUTINE casa_poolout(ktau,veg,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_ncdf_module,    ONLY: def_dims
  USE cable_common_module, ONLY: filename
  USE cable_write_module
  USE cable_checks_module, ONLY: ranges
  USE cable_abort_module,  ONLY: nc_abort
  USE netcdf
  IMPLICIT NONE
  INTEGER,               INTENT(IN)    :: ktau
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

  ! variables for netcdf
  INTEGER, PARAMETER            :: num_dims = 3
  INTEGER                       :: ncid, ncok, mp_restart, mpID
  INTEGER                       :: soID, areaID, laiID, slaID, phaseID
  INTEGER                       :: ClabID, CplantID, ClitterID, CsoilID
  INTEGER                       :: NplantID, NlitterID, NsoilID, NsminID
  INTEGER                       :: PplantID, PlitterID, PsoilID
  INTEGER                       :: PslabID, PssorbID, PsoccID
  INTEGER                       :: CbalID, NbalID, PbalID
  INTEGER,  DIMENSION(num_dims) :: dimID, dim_len
  CHARACTER(LEN=12),DIMENSION(num_dims) :: dim_name 
  CHARACTER(LEN=99)             :: ncfile
  REAL,     DIMENSION(mp)       :: dummy

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol	1	61.3
  ! Andisol	2	103.9
  ! Aridisol	3	92.8
  ! Entisol	4	136.9
  ! Gellisol	5	98.2
  ! Histosol	6	107.6
  ! Inceptisol	7	84.1
  ! Mollisol	8	110.1
  ! Oxisol	9	35.4	
  ! Spodosol	10	41.0	
  ! Ultisol	11	51.5	
  ! Vertisol	12	190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/
   !
   ! estimated based on Yang, Post and Jain (2013)
!   Soiltype     soilnumber soil P(g P/m2  top 50 cm)
!   Alfisol     1       400
!   Andisol     2       426   
!   Aridisol    3       352 
!   Entisol     4       490 
!   Gellisol    5       403 
!   Histosol    6       441 
!   Inceptisol  7       501 
!   Mollisol    8       358   
!   Oxisol      9       96   
!   Spodosol    10      364     
!   Ultisol     11      272   
!   Vertisol    12      430   
!  DATA psorder/400.0,426.0,352.0,490.0,403.0,441.0,501.0,358.0,96.0,364.0,272.0,430.0/
!  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
!  DATA fracpLab/0.07,0.04,0.08,0.10,0.08,0.10,0.12,0.05,0.05,0.06,0.06,0.05/
!  DATA fracPsorb/0.30,0.44,0.69,0.53,0.37,0.14,0.24,0.32,0.15,0.21,0.17,0.35/
!  DATA fracPocc/0.38,0.22,0.18,0.22,0.38,0.42,0.23,0.44,0.60,0.30,0.51,0.48/
!  DATA fracPorg/0.25,0.30,0.05,0.15,0.17,0.34,0.41,0.19,0.20,0.43,0.26,0.12/
!  DATA xpsoil50/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/

  PRINT *, 'Within casa_poolout, mp,ktau = ', mp,ktau
!  nout=103
!  OPEN(nout,file=casafile%cnpepool)
!  PRINT *, 'Opened file ', casafile%cnpepool

  casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
  casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
  casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)
   if(casamet%iveg2(npt)>0 ) then
    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0
      if(casamet%iveg2(npt)==grass) then
         casapool%nplant(npt,wood) = 0.0
         casapool%nlitter(npt,cwd) = 0.0
      endif
    ENDIF 

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%Pplant(npt,:) = casapool%Nplant(npt,:)  &
                             / casapool%ratioNPplant(npt,:)
      casapool%Plitter(npt,:)= casapool%Nlitter(npt,:) &
                             / casapool%ratioNPlitter(npt,:)
      casapool%Psoil(npt,:)  = casapool%Nsoil(npt,:)   &
                             / casapool%ratioNPsoil(npt,:)
      casapool%Psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%Psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%Psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
      if(casamet%iveg2(npt)==grass) then
         casapool%Pplant(npt,wood) = 0.0
         casapool%Plitter(npt,cwd) = 0.0
      endif
    ENDIF
   else
     casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0
     casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
     casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0
     casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
     casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0
     casapool%psoil(npt,:) = 0.0
     casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0
     casapool%psoilocc(npt) = 0.0
     casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
   endif 

    !WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
    !    casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
    !    casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
    !    casabiome%sla(veg%iveg(npt)), phen%phase(npt), casapool%clabile(npt), &
    !    casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
    !    casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
    !    casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
    !    casapool%plitter(npt,:), casapool%psoil(npt,:),         &
    !    casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt),&
    !    casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
  ENDDO

!  CLOSE(nout)

!92    format(5(i6,',',2x),5(f15.6,',',2x),i6,',',2x,100(f15.6,',',2x))

  ! open CABLE restart file (netcdf format)
  ncfile = filename%restart_out
  ncok = NF90_OPEN(ncfile, NF90_WRITE, ncid)
  IF (ncok /= NF90_NOERR) CALL nc_abort(ncok,'Error opening '//TRIM(ncfile))

  ! getting info for the existing dimensions
  ncok = NF90_INQ_DIMID(ncid,'mp',mpID)
  IF(ncok /= NF90_NOERR) THEN
    ncok = NF90_INQ_DIMID(ncid,'mp_patch',mpID)
    IF(ncok /= NF90_NOERR)  CALL nc_abort &
       (ncok,'Error finding mp or mp_patch dimension in ' //TRIM(ncfile))
  END IF
  ncok = NF90_INQUIRE_DIMENSION(ncid,mpID,len=mp_restart)
  IF(ncok /= NF90_NOERR) CALL nc_abort &
       (ncok,'Error finding total number of patches in ' //TRIM(ncfile))
  ! Check that mp_restart = mp from default/met values
  IF(mp_restart /= mp) CALL abort('Number of patches in '// TRIM(ncfile)// &
       ' does not equal to number in default/met file settings.')

  ! get into define mode
  ncok = NF90_REDEF(ncid)
  IF(ncok /= NF90_NOERR) CALL nc_abort &
       (ncok,'Error starting define mode in '//TRIM(ncfile))

  ! define new dimensions
  dim_len(1) = mplant
  dim_len(2) = mlitter
  dim_len(3) = msoil
  dim_name   = (/ "pools_plant", &
                  "pools_litter", &
                  "pools_soil" /)
  CALL def_dims(num_dims, ncid, dimID, dim_len, dim_name )

  ! define new variables
!  ncok = NF90_DEF_VAR(ncid, 'soilOrder', NF90_INT, (/mpID/), soID)
!  IF(ncok /= NF90_NOERR) CALL nc_abort &
!       (ncok,'Error defining soil order in '//TRIM(ncfile))
!  ncok = NF90_PUT_ATT(ncid, soID, "longname", "soil order")
!  IF(ncok /= NF90_NOERR) CALL nc_abort &
!       (ncok,'Error defining attribute of soil order in '//TRIM(ncfile))
!! **** or use **** !
  CALL define_ovar(ncid, soID, 'soilOrder', '-', 'soil order', &
                   .TRUE., 'integer', 0, 0, 0, mpID, 0, .TRUE.)
!  CALL define_ovar(ncid, areaID, 'areacell', '1.0E-9 m2', 'area of tile', &
!                   .TRUE., 'real', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, areaID, 'areacell', 'm2', 'area of tile', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, laiID, 'LAI', '-', 'Leaf Area Index', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, slaID, 'SLA', 'm2', 'Specific Leaf Area', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, phaseID, 'phase', '-', 'phenological phase', &
                   .TRUE., 'integer', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, ClabID, 'Clabile', 'gC/m2', 'labile C pool', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, CplantID, 'CASA_Cplant', 'gC/m2', 'plant C pools', &
                   .TRUE., dimID(1), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, ClitterID, 'Clitter', 'gC/m2', 'litter C pools', &
                   .TRUE., dimID(2), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, CsoilID, 'CASA_Csoil', 'gC/m2', 'soil C pools', &
                   .TRUE., dimID(3), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, NplantID, 'Nplant', 'gN/m2', 'plant N pools', &
                   .TRUE., dimID(1), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, NlitterID, 'Nlitter', 'gN/m2', 'litter N pools', &
                   .TRUE., dimID(2), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, NsoilID, 'Nsoil', 'gN/m2', 'soil N pools', &
                   .TRUE., dimID(3), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, NsminID, 'Nsoilmin', 'gN/m2', 'mineral N in soil', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PplantID, 'Pplant', 'gP/m2', 'plant P pools', &
                   .TRUE., dimID(1), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PlitterID, 'Plitter', 'gP/m2', 'litter P pools', &
                   .TRUE., dimID(2), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PsoilID, 'Psoil', 'gP/m2', 'soil P pools', &
                   .TRUE., dimID(3), 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PslabID, 'Psoillab', 'gP/m2', 'labile P in soil', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PssorbID, 'Psoilsorb', 'gP/m2', 'adsorbed P in soil', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PsoccID, 'Psoilocc', 'gP/m2', 'occluded P in soil', &
                   .TRUE., 'r2', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, CbalID, 'sumCbal', 'gC/m2', 'Accumulated C balance', &
                   .TRUE., 'real', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, NbalID, 'sumNbal', 'gN/m2', 'Accumulated N balance', &
                   .TRUE., 'real', 0, 0, 0, mpID, 0, .TRUE.)
  CALL define_ovar(ncid, PbalID, 'sumPbal', 'gP/m2', 'Accumulated P balance', &
                   .TRUE., 'real', 0, 0, 0, mpID, 0, .TRUE.)

  ! End netcdf define mode:
  ncok = NF90_ENDDEF(ncid)
  IF(ncok /= NF90_NOERR) CALL nc_abort(ncok, 'Error redefining restart file '  &
                 //TRIM(filename%restart_out)// '(SUBROUTINE casa_poolout)')

  CALL write_ovar(ncid, soID, 'soilOrder', REAL(casamet%isorder,4),  &
                   ranges%SoilOrder, .TRUE., 'integer', .TRUE.)
  CALL write_ovar(ncid, areaID, 'areacell', casamet%areacell, ranges%area, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, laiID, 'LAI', casamet%glai, ranges%lai, &
                  .TRUE., 'cnp', .TRUE.)
  dummy(:) = casabiome%sla(veg%iveg(:))
  CALL write_ovar(ncid, slaID, 'SLA', dummy, ranges%sla, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, phaseID, 'phase', REAL(phen%phase,4), ranges%phase, &
                  .TRUE., 'integer', .TRUE.)
  CALL write_ovar(ncid, ClabID, 'Clabile', casapool%clabile, ranges%Clab, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, CplantID, 'CASA_Cplant', casapool%cplant,ranges%Cplant,&
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, ClitterID, 'Clitter', casapool%clitter, ranges%Clitter,&
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, CsoilID, 'CASA_Csoil', casapool%csoil, ranges%Csoil, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, NplantID, 'Nplant', casapool%nplant, ranges%Nplant, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, NlitterID, 'Nlitter', casapool%nlitter, ranges%Nlitter,&
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, NsoilID, 'Nsoil', casapool%nsoil, ranges%Nsoil, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, NsminID, 'Nsoilmin', casapool%nsoilmin, ranges%Nsmin, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PplantID, 'Pplant', casapool%pplant, ranges%Pplant, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PlitterID, 'Plitter', casapool%plitter, ranges%Plitter,&
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PsoilID, 'Psoil', casapool%psoil, ranges%Psoil, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PslabID, 'Psoillab', casapool%psoillab, ranges%Pslab, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PssorbID, 'Psoilsorb', casapool%psoilsorb,  &
                  ranges%Pssorb, .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, PsoccID, 'Psoilocc', casapool%psoilocc, ranges%Psocc, &
                  .TRUE., 'cnp', .TRUE.)
  CALL write_ovar(ncid, CbalID, 'sumCbal', casabal%sumcbal, ranges%Cbal, &
                  .TRUE., 'real', .TRUE.)
  CALL write_ovar(ncid, NbalID, 'sumNbal', casabal%sumnbal, ranges%Nbal, &
                  .TRUE., 'real', .TRUE.)
  CALL write_ovar(ncid, PbalID, 'sumPbal', casabal%sumpbal, ranges%Pbal, &
                  .TRUE., 'real', .TRUE.)

  ! Close restart file
  ncok = NF90_CLOSE(ncid)

END SUBROUTINE casa_poolout


