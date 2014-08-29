all: cable libcable.a

FC          = mpif90
FFLAGS     += -fpp -i8 -r8

coresrc    := $(shell find core    -name '*.F90')
offlinesrc := $(shell find offline -name '*.F90')
umsrc      := $(shell find UM      -name '*.F90')

libcable.a: $(patsubst %.F90,%.o,$(coresrc) $(umsrc))
	echo $^

cable: $(patsubst %.F90,%.o,$(coresrc) $(offlinesrc))

%.o:%.F90
	$(FC) $(FFLAGS) -c -o $@ $<

# Generate module dependencies (compiles the module as a side-effect)
# If the build fails remove the dep file
%.dep:%.F90
	$(FC) $(FFLAGS) -c -gen-dep -o $*.o $< > $@ || $(RM) $@

-include $(patsubst %.F90,%.dep,$(coresrc) $(umsrc) $(offlinesrc))
