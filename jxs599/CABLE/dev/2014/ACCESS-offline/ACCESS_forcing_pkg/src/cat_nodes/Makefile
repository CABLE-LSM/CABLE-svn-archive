#program specific structure
AUX= subr/

#compiler independent
FLINK = -c 

#compiler specific

FC = ifort
FFLAGS = -O0 -fpp -r8 -i8

#FC = gfortran 
#FFLAGS = -O0 -x f95-cpp-input  -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -ffree-form -ffree-line-length-0    

OBJ= cat_common.o cat_read.o cat_write.o cat_nodes.o

cat_Nnodes: cat_Nnodes.f90
	$(FC) $(FFLAGS) -c $(AUX)cat_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)cat_read.f90	
	$(FC) $(FFLAGS) -c $(AUX)cat_write.f90	
	$(FC) $(FFLAGS) -c $(AUX)cat_nodes.f90	
	$(FC) $(FFLAGS) -c cat_Nnodes.f90 
	$(FC) -o cat_Nnodes cat_Nnodes.o $(OBJ) 
clean:
	rm -f *.o *~ *.mod





##FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng
#cherax?
#FPGPLOTLIBS = -L/cs/home/csdar/srb001/Public/pgplot -lpgplot -L/usr/lib -lpng -L/usr/lib -lz -L/usr/X11R6/lib/ -lX11  
# xe
#FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  


