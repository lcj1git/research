#
#  makefile to compile phfmol.x and related programs
#
#  Usage:  make all      - compile all programs
#          make phfmol.x - compile only main program
#          make clean    - clean object files
#
#  CAJH, 01.2013


#  Fortran compiler and flags

SHELL   = /bin/bash
FC      = mpif90
FCFLAGS = -g -fbacktrace -Wall -O3 -fopenmp
LDFLAGS = -fopenmp -llapack -lblas


#  rule to compile source code

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FCFLAGS) -c $< -o $@


#  define objects to build

all: phfmol.x


#  rule to clean object and module files

clean:
	rm -f *.o *.mod


#  ==  phfmol.x  ==

MAIN_phf = phfmol.o
OBJT_phf = constants.o  i2sint.o  grid.o  spin.o  util.o \
	   linesrch.o  lbfgs.o  linalg.o  trdnmt.o  purcart.o \
	   erictr.o  spatpr.o  iguess.o  hamilt.o  diagh11.o \
	   phfscf.o

phfmol.x: $(OBJT_phf) $(MAIN_phf)
	$(FC) $(LDFLAGS) $(MAIN_phf) $(OBJT_phf) -o phfmol.x


