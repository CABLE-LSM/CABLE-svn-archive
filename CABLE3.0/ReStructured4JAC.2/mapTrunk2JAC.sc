###------------------------------------------------------------------------
###r7308 | jxs599 | 2020-07-03 19:11:07 +1000 (Fri, 03 Jul 2020) | 1 line
###Changed paths:
svn mkdir params
svn mkdir science
svn mkdir science/albedo
svn mkdir science/canopy
svn mkdir science/casa-cnp
svn mkdir science/coupled
svn mkdir science/gw_hydro
svn mkdir science/misc
svn mkdir science/pop
svn mkdir science/radiation
svn mkdir science/roughness
svn mkdir science/sli
svn mkdir science/soilsnow
svn mkdir util
##
###------------------------------------------------------------------------
###r7309 | jxs599 | 2020-07-03 19:26:15 +1000 (Fri, 03 Jul 2020) | 1 line
###Changed paths:
###mv to renamed directories
svn mv UM coupled
##
###------------------------------------------------------------------------
###r7310 | jxs599 | 2020-07-03 19:50:38 +1000 (Fri, 03 Jul 2020) | 1 line
###Changed paths:
###re-distribute caontents of core/utils/
svn mv core/utils/cable_pft_params.F90 offline/
svn mv core/utils/cable_soil_params.F90 offline/
svn mv core/biogeochem/CASAONLY_LUC.F90 science/casa-cnp/
svn mv core/biogeochem/CSIRO_BSD_MIT_License_v2.0_CABLE.txt science/casa-cnp/
svn mv core/biogeochem/cable_iovars_CMIP6.F90 science/casa-cnp/
svn mv core/biogeochem/cable_phenology.F90 science/casa-cnp/
svn mv core/biogeochem/casa_cable.F90 science/casa-cnp/
svn mv core/biogeochem/casa_cnp.F90 science/casa-cnp/
svn mv core/biogeochem/casa_dimension.F90 science/casa-cnp/
svn mv core/biogeochem/casa_inout.F90 science/casa-cnp/
svn mv core/biogeochem/casa_param.F90 science/casa-cnp/
svn mv core/biogeochem/casa_variable.F90 science/casa-cnp/
svn mv core/biogeochem/spincasacnp.F90 science/casa-cnp/
svn mv core/biogeochem/POP.F90 science/pop/
svn mv core/biogeochem/POPLUC.F90 science/pop/
svn mv core/biogeochem/pop_io.F90 science/pop/
svn mv core/biogeochem/pop_mpi.F90 science/pop/
svn mv core/utils/diag util/diag 

svn rm core/biogeochem
svn rm core/utils

#------------------------------------------------------------------------
#r7312 | jxs599 | 2020-07-05 10:44:12 +1000 (Sun, 05 Jul 2020) | 1 line
#Changed paths:
#mods to build script (serial ONLY) to gather files
#------------------------------------------------------------------------
#r7311 | jxs599 | 2020-07-03 20:17:03 +1000 (Fri, 03 Jul 2020) | 1 line
#Changed paths:
#re-distribute contents of core/biogeophys/
svn mv core/biogeophys/cable_cbm.F90 offline/cable_cbm.F90
svn mv core/biogeophys/cable_constants.F90 params/cable_constants.F90 
svn mv core/biogeophys/cable_data.F90 params/cable_data.F90 
svn mv core/biogeophys/CSIRO_BSD_MIT_License_v2.0_CABLE.txt science/CSIRO_BSD_MIT_License_v2.0_CABLE.txt 
svn mv core/biogeophys/cable_albedo.F90 science/albedo/cable_albedo.F90 
svn mv core/biogeophys/cable_canopy.F90 science/canopy/cable_canopy.F90 
svn mv core/biogeophys/cable_gw_hydro.F90 science/gw_hydro/cable_gw_hydro.F90 
svn mv core/biogeophys/cable_psm.F90 science/gw_hydro/cable_psm.F90 
svn mv core/biogeophys/cable_air.F90 science/misc/cable_air.F90 
svn mv core/biogeophys/cable_carbon.F90 science/misc/cable_carbon.F90 
svn mv core/biogeophys/cable_climate.F90 science/misc/
svn mv core/biogeophys/cable_radiation.F90 science/radiation/cable_radiation.F90 
svn mv core/biogeophys/cable_roughness.F90 science/roughness/cable_roughness.F90 
svn mv core/biogeophys/cable_sli_main.F90 science/sli/cable_sli_main.F90 
svn mv core/biogeophys/cable_sli_numbers.F90 science/sli/cable_sli_numbers.F90 
svn mv core/biogeophys/cable_sli_roots.F90 science/sli/cable_sli_roots.F90 
svn mv core/biogeophys/cable_sli_solve.F90 science/sli/cable_sli_solve.F90 
svn mv core/biogeophys/cable_sli_utils.F90 science/sli/cable_sli_utils.F90 
svn mv core/biogeophys/cbl_soilsnow_init_special.F90 science/soilsnow/cbl_soilsnow_init_special.F90 
svn mv core/biogeophys/cbl_soilsnow_main.F90 science/soilsnow/cbl_soilsnow_main.F90 
svn mv core/biogeophys/cbl_soilsnow_subrs.F90 science/soilsnow/cbl_soilsnow_subrs.F90 
svn mv core/biogeophys/cable_common.F90 util/cable_common.F90 
svn mv core/biogeophys/cable_define_types.F90 util/cable_define_types.F90 

svn rm core/biogeophys/
svn rm core


