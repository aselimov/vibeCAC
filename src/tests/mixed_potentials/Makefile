FC=mpif90

#gfortran flags
ARCH := $(shell uname -m)
FFLAGS=-O2 -g -fbacktrace -fcheck=all
ifeq ($(ARCH),x86_64)
FFLAGS += -mcmodel=large
endif
#FFLAGS=-O3 -g 
TEST=mixed_potentials
OBJECTS=$(wildcard *.o)
MAINS=at_main.o read_dump.o
AT_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))
CG_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))

.SUFFIXES: .c .f .f90 .F90 .o

.PHONY: run_$(TEST) 
run_$(TEST): $(TEST)_at_test 
	./$(TEST)_at_test
	mpirun -np 4 ./$(TEST)_at_test
	mpirun -np 8 ./$(TEST)_at_test

$(TEST)_at_test: at_main.o $(AT_OBJECTS)
	$(FC) $(FFLAGS) $(AT_OBJECTS) read_dump.o at_main.o -o $@

%.o : %.mod

.f90.o: 
	$(FC) $(FFLAGS) -c $< 	

at_main.o :read_dump.o


