FC=ifort
FCFLAGS=-O2 -static
CC=gcc
CCFLAGS=

SUB_DIR_OBJECTS= gsl_interface_c.o gsl_interface_f.o global.o mesh.o splines.o siesta_pass.o
LIBRARIES= -lgsl -lgslcblas -lm 
CINCLUDES= -I.

#all: $(SUB_DIR_OBJECTS)
all: siesta_pass.a

BUILDDIR=.

include make.depend

siesta_pass.a:$(SUB_DIR_OBJECTS)
	ar rc siesta_pass.a $(SUB_DIR_OBJECTS);cp siesta_pass.a ../.;cp *.o ../.;cp *.mod ../.

$(BUILDDIR)/%:%.o
	$(FC) $(FCFLAGS) -o $@ $(LIBRARIES) $^ 

$(BUILDDIR)/%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
$(BUILDDIR)/%.o: %.F90
	$(FC) $(FCFLAGS) -c $<
$(BUILDDIR)/%.o: %.c
	$(CC) $(CCFLAGS) -c $(CINCLUDES) $<

clean:
	rm -f *.o *.mod *.MOD $(SUB_DIR_OBJECTS)
