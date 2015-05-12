# CCPG - Cross Code Pseudopotential Generator

CCPG is a modification to the [QuantumEspresso](http://www.quantum-espresso.org/) pseudopotential generator to allow it to also produce pseudopotentials in the [Siesta](http://departments.icmab.es/leem/siesta/) format.

### Building and compiling info
Ok so compile without gfortran, or open mpi.
compile and run with gsl

must set FC to

export FC=ifort

build with:
```cd siesta;make;cd ../;make thing to compile though```

## NB! pay attention to the pseudotype number
