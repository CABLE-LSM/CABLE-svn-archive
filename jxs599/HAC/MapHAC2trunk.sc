bincp science_cable/ /short/p66/jxs599/CABLE_validation/src/226_test_albedo/science
bincp control/cable/ /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/
bincp params/cable/  /short/p66/jxs599/CABLE_validation/src/226_test_albedo/params/

bincp util/cable/diag /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util/diag

cp util/cable/cable_data.F90          /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util/
cp util/cable/cable_pft_params.F90    /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util/
cp util/cable/cable_soil_params.F90   /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util/
cp util/cable/cable_common.F90        /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util/
cp util/cable/cable_constants.F90     /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util
cp util/cable/cbl_masks.F90           /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util
cp util/cable/cable_define_types.F90  /short/p66/jxs599/CABLE_validation/src/226_test_albedo/util

cp util/cable/cbl_radiation_albedo.F90      /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cbl_allocate_types.F90        /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cable_pack_mod.F90            /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cable_um_tech.F90             /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cable_jules_links_mod.F90     /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cbl_cable_wide.F90            /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/
cp util/cable/cable_gather_UM_data_decs.F90 /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/util/

bincp control/cable/initialization/ /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/initialization/
bincp initialisation/cable /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/initialization/

cp initialisation/standalone/ancillaries/init_cable.inc     /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/initialization/
cp initialisation/standalone/params/init_soilparm_cable.inc /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/initialization/
cp initialisation/standalone/init_cable_mod.F90             /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/initialization/

cp control/shared/cable_surface_types_mod.F90 /short/p66/jxs599/CABLE_validation/src/226_test_albedo/coupled/






