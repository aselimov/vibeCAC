.DEFAULT_GOAL := all

FC=mpif90
OBJDIR=obj
SRCS := $(wildcard *.f90)
OBJECTS := $(addprefix $(OBJDIR)/,$(SRCS:%.f90=%.o))
ARCH := $(shell uname -m)
FFLAGS=-Wall -O2 -g -fbacktrace -fcheck=all -J$(OBJDIR)
ifeq ($(ARCH),x86_64)
FFLAGS += -mcmodel=large
endif
#
#FFLAGS=-O3 -g -J$(OBJDIR) 

#----------------- DEPENDENCIES -----------------#
# GENERATED USING https://github.com/ZedThree/fort_depend.py **requires python3**
# > pip install fortdepend
# > fortdepend -o Makefile.dep -i mpi -b obj/
include Makefile.dep

#----------------- DEFAULTS -----------------#
	
all: CAC

.PHONY: deps
CAC: $(OBJECTS) $(OBJDIR)/main.o
	$(FC) $(FFLAGS) $(OBJECTS) -o $@

$(OBJDIR)/%.o: %.f90
	@mkdir -p $(@D) 
	$(FC) $(FFLAGS) -c -o $@ $< 

.f90.o: 
	$(FC) $(FFLAGS) -c $< 


#deps:
#	@fortdepend -o Makefile.dep -i mpi -b obj -w

deps:
	@makedepf90 -b obj *.f90 > Makefile.dep
#----------------- TEST RULES -----------------#

COMMS_TEST_OBJECTS := $(addprefix $(OBJDIR)/,$(COMMS_TEST))

MULTICOMM_TEST:=$(COMMS_TEST)
MULTICOMM_TEST_OBJECTS:=$(addprefix $(OBJDIR)/,$(MULTICOMM_TEST))

ELEMENTS_TEST=parameters.o elements.o box.o errors.o math.o logger.o
ELEMENTS_TEST_OBJECTS := $(addprefix $(OBJDIR)/,$(ELEMENTS_TEST))


DUMP_TEST=$(POT_TEST) dump.o
DUMP_TEST_OBJECTS := $(addprefix $(OBJDIR)/,$(DUMP_TEST))



.PHONY:	tests
tests: test_cohesive_energy test_pair_force test_multi_pair_force test_fs_pair_force test_cg_calc 

test_cohesive_energy: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/cohesive_energy
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/cohesive_energy
	$(MAKE) -C ./tests/cohesive_energy

test_gsf: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/gsf
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/gsf
	$(MAKE) -C ./tests/gsf

test_pair_force: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/pair_force
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/pair_force
	$(MAKE) -C ./tests/pair_force

test_multi_pair_force: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/multi_pair_force
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/multi_pair_force
	$(MAKE) -C ./tests/multi_pair_force

test_fs_pair_force: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/fs_pair_force
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/fs_pair_force
	$(MAKE) -C ./tests/fs_pair_force

test_morse_force: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/morse_force
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/morse_force
	$(MAKE) -C ./tests/morse_force

test_cg_calc: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/cg_calc
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/cg_calc
	$(MAKE) -C ./tests/cg_calc


test_mixed_potential: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	@cp $(TESTOBJECTS) ./tests/mixed_potentials
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./tests/mixed_potentials
	$(MAKE) -C ./tests/mixed_potentials


.PHONY: unit-test
unit-test: $(OBJECTS)
	$(eval TESTOBJECTS := $(filter-out $(OBJDIR)/main.o,$(OBJECTS)))
	$(MAKE) -C ./unit_tests clean
	@cp $(TESTOBJECTS) ./unit_tests
	@cp $(subst .o,.mod,$(TESTOBJECTS)) ./unit_tests
	$(MAKE) -C ./unit_tests run UNIT_TEST=$(if $(UNIT_TEST),$(UNIT_TEST),math) FFLAGS="$(if $(UNIT_TEST_FFLAGS),$(UNIT_TEST_FFLAGS),-Wall -O0 -g -fcheck=all -fbacktrace)"

.PHONY: unit-test-coverage
unit-test-coverage:
	$(MAKE) clean
	$(MAKE) unit-test UNIT_TEST=$(if $(UNIT_TEST),$(UNIT_TEST),math) FFLAGS="-Wall -O0 -g -fcheck=all -fbacktrace --coverage -J$(OBJDIR)" UNIT_TEST_FFLAGS="-Wall -O0 -g -fcheck=all -fbacktrace --coverage"
	$(MAKE) -C ./unit_tests coverage-report UNIT_TEST=$(if $(UNIT_TEST),$(UNIT_TEST),math)

#----------------- CLEAN UP -----------------#

.PHONY: clean
clean:
	$(RM) *.mod *.o
	$(RM) $(OBJDIR)/*.mod $(OBJDIR)/*.o CAC
	@$(RM) -rf obj/

.PHONY: cleantest
cleantest:
	$(RM) tests/*/*.o tests/*/*.mod tests/*/*_test

.PHONY: clean-all
clean-all: clean cleantest

# DEBUGGING VARIABLE PRINT
print-%  : ; @echo $* = $($*)
