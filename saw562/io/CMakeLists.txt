cmake_minimum_required(VERSION 3.4)
enable_language(Fortran)

project(cable)

set(CMAKE_C_COMPILER       mpicc)
set(CMAKE_Fortran_COMPILER mpif90)

if (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -Wextra -fimplicit-none -fbacktrace")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check --coverage")

elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -warn all -fp-model precise -diag-error=6717 -diag-disable=7712")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -xHost")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check all -ftz -fpe0 -check noarg_temp_created")

endif(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)

add_library(offline
    core/biogeochem/CASAONLY_LUC.F90
    core/biogeochem/POP.F90
    core/biogeochem/POPLUC.F90
    core/biogeochem/cable_phenology.F90
    core/biogeochem/casa_cable.F90
    core/biogeochem/casa_cnp.F90
    core/biogeochem/casa_inout.F90
    core/biogeochem/casa_variable.F90
    core/biogeochem/pop_io.F90
    core/biogeochem/pop_mpi.F90
    core/biogeochem/spincasacnp.F90
    core/biogeophys/cable_air.F90
    core/biogeophys/cable_albedo.F90
    core/biogeophys/cable_canopy.F90
    core/biogeophys/cable_carbon.F90
    core/biogeophys/cable_cbm.F90
    core/biogeophys/cable_climate.F90
    core/biogeophys/cable_common.F90
    core/biogeophys/cable_constants.F90
    core/biogeophys/cable_data.F90
    core/biogeophys/cable_define_types.F90
    core/biogeophys/cable_diag.F90
    core/biogeophys/cable_radiation.F90
    core/biogeophys/cable_roughness.F90
    core/biogeophys/cable_sli_main.F90
    core/biogeophys/cable_sli_numbers.F90
    core/biogeophys/cable_sli_roots.F90
    core/biogeophys/cable_sli_solve.F90
    core/biogeophys/cable_sli_utils.F90
    core/biogeophys/cable_soilsnow.F90
    offline/cable_LUC_EXPT.F90
    offline/cable_abort.F90
    offline/cable_checks.F90
    offline/cable_cru_TRENDY.F90
    offline/cable_initialise.F90
    offline/cable_input.F90
    offline/cable_iovars.F90
    offline/cable_mpicommon.F90
    offline/cable_mpimaster.F90
    offline/cable_mpiworker.F90
    offline/cable_output.F90
    offline/cable_parameters.F90
    offline/cable_plume_mip.F90
    offline/cable_read.F90
    offline/cable_weathergenerator.F90
    offline/cable_write.F90
    comms.f90
    save.f90
    cable_comms.f90
    )

add_executable(cable offline/cable_driver.F90)
target_link_libraries(cable offline -lnetcdff -lnetcdf)

add_executable(cable_mpi offline/cable_mpidrv.F90)
target_link_libraries(cable_mpi offline -lnetcdff -lnetcdf)

add_executable(test_comms test_comms.f90)
target_link_libraries(test_comms offline)

add_executable(diff diff.f90)
target_link_libraries(diff offline)
