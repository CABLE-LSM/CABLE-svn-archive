# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(parallel_mod.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(parallel_mod.o): $(ARC)(logging_mod.o) $(ARC)(io_constants.o) $(ARC)(string_utils_mod.o) $(ARC)(grid_utils_mod.o) $(ARC)(model_grid_mod.o) scatter_land_field.inc $(ARC)(theta_field_sizes.o) $(ARC)(ancil_info.o) is_master_task.inc decompose_domain.inc gather_land_field.inc


# End of file.
