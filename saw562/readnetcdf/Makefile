all: cable libcable.a

GENDEPS += -gen-dep -syntax-only

coresrc    := $(shell find core    -name '*.F90')
offlinesrc := $(shell find offline -name '*.F90')
umsrc      := $(shell find UM      -name '*.F90')

libcable.a: $(patsubst %.F90,%.o,$(coresrc) $(umsrc))

cable: $(patsubst %.F90,%.o,$(coresrc) $(offlinesrc))

%.o:%.F90
	$(FC) $(FFLAGS) -c -o $@ $<

%.d:%.F90
	$(FC) $(FFLAGS) -syntax-only -gen-dep $<

-include $(patsubst %.F90,%.d,$(coresrc) $(umsrc) $(offlinesrc))
