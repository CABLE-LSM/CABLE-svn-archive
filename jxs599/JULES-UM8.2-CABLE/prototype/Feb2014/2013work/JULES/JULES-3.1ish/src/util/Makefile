# This is a generated file - do not modify
        
include $(JULESDIR)/Makefile.common.mk

OBJS = $(ARC)(dictionary_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(datetime_mod.o) $(ARC)(logging_mod.o) $(ARC)(templating_mod.o)

all : $(OBJS)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

##############################
## Dependencies             ##
##############################
$(ARC)(dictionary_mod.o): $(ARC)(logging_mod.o)
$(ARC)(string_utils_mod.o): $(ARC)(io_constants.o)
$(ARC)(datetime_mod.o): $(ARC)(logging_mod.o) $(ARC)(string_utils_mod.o)
$(ARC)(logging_mod.o): $(ARC)(io_constants.o)
$(ARC)(templating_mod.o): $(ARC)(datetime_mod.o) $(ARC)(string_utils_mod.o) $(ARC)(io_constants.o)


# End of file.
