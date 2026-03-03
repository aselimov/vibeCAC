FC=mpif90

#gfortran flags
ARCH := $(shell uname -m)
FFLAGS=-O0 -g -fbacktrace -fcheck=all
ifeq ($(ARCH),x86_64)
FFLAGS += -mcmodel=large
endif
TEST=pair_force
OBJECTS=$(wildcard *.o)
MAINS=at_main.o cg_main.o
AT_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))
CG_OBJECTS=$(filter-out $(MAINS),$(OBJECTS))

.SUFFIXES: .c .f .f90 .F90 .o

.PHONY: run_$(TEST) 
run_$(TEST): $(TEST)_at_test  
	./$(TEST)_at_test

$(TEST)_at_test: at_main.o $(AT_OBJECTS)
	$(FC) $(FFLAGS) $(AT_OBJECTS) at_main.o -o $@

%.o : %.mod

.f90.o: 
	$(FC) $(FFLAGS) -c $< 	





