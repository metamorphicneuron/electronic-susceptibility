# lindhard_chi: Lindhard susceptibility (Fortran/MPI) and its Brillouin-zone overlay
#
#   make                      builds bin/lindhard_chi.x and bin/chi_bz.x with mpif90
#   make FC=mpiifx FFLAGS=-O3 Intel compilers (or mpiifort)
#   make debug                build with run-time checks
#   make test                 compare against numpy (needs python3 + numpy)
#   make test MPIRUN="srun -n 2"
#   make clean                remove bin/ and build/

FC      = mpif90
FFLAGS  = -O3
MPIRUN  = mpirun -np 2

all: bin/lindhard_chi.x bin/chi_bz.x

# compiled inside build/, so that the module files (.mod) stay there
bin/lindhard_chi.x: src/lindhard_chi.f90
	mkdir -p bin build
	cd build && $(FC) $(FFLAGS) -o ../bin/lindhard_chi.x ../src/lindhard_chi.f90

# serial; built with the same compiler
bin/chi_bz.x: src/chi_bz.f90
	mkdir -p bin build
	cd build && $(FC) $(FFLAGS) -o ../bin/chi_bz.x ../src/chi_bz.f90

debug: FFLAGS = -O0 -g -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
debug: clean all

test: all
	python3 tests/run_test.py bin/lindhard_chi.x "$(MPIRUN)"

clean:
	rm -rf bin build

.PHONY: all debug test clean

