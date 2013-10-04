FC = ifort
#FC = g95
XFLAGS = 
OPT = 
NCDIR = /usr/local/intel/lib/     # netcdf library address on shine-cl

.SUFFIXES : .o .f90

.f90.o:
	$(FC)  -c $(OPT) $(XFLAGS) $<

PROG = convertMetTB3

OBJS = convertMetTB3.o

#convertrun: $(PROG)
#	$(PROG)

$(PROG): $(OBJS)
	$(FC) -o $(PROG) $(OBJS) $(XFLAGS) -L$(NCDIR) -lnetcdf

convertMetTB3.o: convertMetTB3.f90
	$(FC) -c $(OPT) $(XFLAGS) convertMetTB3.f90

clean:
	rm -f $(PROG) convertMetTB3.o 
	rm -f fort.*
	ln -s /usr/local/intel/include/netcdf.mod ./netcdf.mod   # for shine-cl
