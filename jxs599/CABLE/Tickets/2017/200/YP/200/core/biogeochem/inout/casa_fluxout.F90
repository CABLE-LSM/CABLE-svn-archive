! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)
!SUBROUTINE casa_fluxout(myear,clitterinput,csoilinput)
  USE cable_def_types_mod
!  USE cableDeclare, ONLY: veg, soil
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
!  USE casaDeclare
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  INTEGER,               INTENT(IN)    :: myear
!  REAL(r_2),    INTENT(IN) :: clitterinput(mp,3),csoilinput(mp,3)

  ! local variables
  INTEGER  npt,nout
  REAL(r_2) xyear, totGPP, totNPP

  totGPP =0.0
  totNPP =0.0
  nout=104
  xyear=1.0/FLOAT(myear)
  casabal%FCgppyear=casabal%FCgppyear * xyear
  casabal%FCnppyear=casabal%FCnppyear * xyear
  casabal%FCrmleafyear=casabal%FCrmleafyear * xyear
  casabal%FCrmwoodyear=casabal%FCrmwoodyear * xyear
  casabal%FCrmrootyear=casabal%FCrmrootyear * xyear
  casabal%FCrgrowyear=casabal%FCrgrowyear * xyear
  casabal%FCrsyear=casabal%FCrsyear * xyear
  casabal%FCneeyear=casabal%FCneeyear * xyear
  casabal%FNdepyear=casabal%FNdepyear * xyear
  casabal%FNfixyear=casabal%FNfixyear * xyear
  casabal%FNsnetyear=casabal%FNsnetyear * xyear
  casabal%FNupyear=casabal%FNupyear * xyear
  casabal%FNleachyear=casabal%FNleachyear * xyear
  casabal%FNlossyear=casabal%FNlossyear * xyear
  casabal%FPweayear=casabal%FPweayear * xyear
  casabal%FPdustyear=casabal%FPdustyear * xyear
  casabal%FPsnetyear=casabal%FPsnetyear * xyear
  casabal%FPupyear=casabal%FPupyear * xyear
  casabal%FPleachyear=casabal%FPleachyear * xyear
  casabal%FPlossyear=casabal%FPlossyear * xyear
!  clitterinput = clitterinput * xyear
!  csoilinput   = csoilinput   * xyear

  print *, 'writing CNP fluxes out to file ', casafile%cnpflux
  OPEN(nout,file=casafile%cnpflux)
    DO npt =1,mp
      SELECT CASE(icycle)
      CASE(1)

        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%Fcnppyear(npt),  &
            casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
            casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &
!            clitterinput(npt,:),csoilinput(npt,:)

      CASE(2)
        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%FCnppyear(npt),                                 &
            casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
            casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
            casabal%FCrsyear(npt), casabal%FCneeyear(npt),          &
!        clitterinput(npt,:),csoilinput(npt,:), &
        casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

      CASE(3)
        WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
        casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
        casabal%FCnppyear(npt),                                  &
        casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
        casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
        casabal%FCrsyear(npt),   casabal%FCneeyear(npt),         &
!        clitterinput(npt,:),csoilinput(npt,:), &
       casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt),&
       casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt),&
       casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt),&
       casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)

      END SELECT
      totGPP = totGPP+casabal%Fcgppyear(npt)* casamet%areacell(npt)


      totNPP = totNPP+casabal%Fcnppyear(npt)* casamet%areacell(npt)
    ENDDO

    print *, 'totGPP global = ', totGPP*(1.0e-15)
    print *, 'totNPP global = ', totNPP*(1.0e-15)
  CLOSE(nout)
92    format(5(i6,',',2x),100(f15.6,',',2x))
END SUBROUTINE casa_fluxout


