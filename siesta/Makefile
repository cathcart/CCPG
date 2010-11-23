FC=ifort
FCFLAGS=-O2
CC=gcc
CCFLAGS=

SUB_DIR_OBJECTS=gsl_interface_c.o gsl_interface_f.o global.o mesh.o splines.o siesta_pass.o
LIBRARIES= -lgsl -lgslcblas -lm 
CINCLUDES= -I.

all: $(SUB_DIR_OBJECTS)

include make.depend


%:%.o
	$(FC) $(FCFLAGS) -o $@ $(LIBRARIES) $^ 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
%.o: %.F90
	$(FC) $(FCFLAGS) -c $<
%.o: %.c
	$(CC) $(CCFLAGS) -c $(CINCLUDES) $<

clean:
	rm -f *.o *.mod *.MOD $(SUB_DIR_OBJECTS)
