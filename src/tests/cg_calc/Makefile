FC=mpif90

#gfortran flags
ARCH := $(shell uname -m)
FFLAGS=-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal
ifeq ($(ARCH),x86_64)
FFLAGS += -mcmodel=large
endif
TEST=cohesive_energy
OBJECTS=$(wildcard *.o)
MAINS=cg_main.o
AT_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))
CG_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))

.SUFFIXES: .c .f .f90 .F90 .o

.PHONY: run_$(TEST) 
run_$(TEST): $(TEST)_cg_test
	./$(TEST)_cg_test

$(TEST)_cg_test: cg_main.o $(CG_OBJECTS)
	$(FC) $(FFLAGS) $(CG_OBJECTS) cg_main.o -o $@

%.o : %.mod

.f90.o: 
	$(FC) $(FFLAGS) -c $< 	





