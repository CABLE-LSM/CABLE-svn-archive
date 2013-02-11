#!/bin/csh
ar r libcable.a cable_explicit_driver.o cable_implicit_driver.o cable_rad_driver.o cable_hyd_driver.o cable_common_module.o cable_diag.o cable_define_dimensions.o cable_math_constants.o cable_other_constants.o cable_photosynthetic_constants.o cable_physical_constants.o cable_define_types.o cable_variables.o cable_iovars.o cable_soilsnow.o cable_air.o cable_albedo.o cable_radiation.o cable_roughness.o cable_carbon.o cable_canopy.o cable_cbm.o cable_um_tech.o cable_um_init_subrs.o cable_um_init.o 

/bin/cp libcable.a ~/LIBCABLE

