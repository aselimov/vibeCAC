.SUFFIXES:

FC = mpif90
GCOV ?= $(shell command -v gcov-15 2>/dev/null || command -v gcov)
OBJDIR ?= obj
FFLAGS ?= -Wall -O0 -g -fcheck=all -fbacktrace -J$(OBJDIR)

UNIT_TEST ?= math

CORE_OBJECTS = parameters.o
ifeq ($(UNIT_TEST),math)
CORE_OBJECTS += math.o
endif
ifeq ($(UNIT_TEST),str)
CORE_OBJECTS += str.o
endif

ifeq ($(UNIT_TEST),time)
CORE_OBJECTS += logger.o errors.o time.o
endif

ifeq ($(UNIT_TEST),box)
CORE_OBJECTS += logger.o errors.o box.o
endif

ifeq ($(UNIT_TEST),dump_logic)
CORE_OBJECTS += dump_logic.o
endif

ifeq ($(UNIT_TEST),input_validation)
CORE_OBJECTS += input_validation.o
endif

ifeq ($(UNIT_TEST),atom_types)
CORE_OBJECTS += str.o logger.o atom_types.o
endif

TEST_OBJECTS = fixtures.o assertions.o test_$(UNIT_TEST).o
EXECUTABLE = ut_$(UNIT_TEST)

.PHONY: all run coverage-report clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(CORE_OBJECTS) $(TEST_OBJECTS)
	$(FC) $(FFLAGS) $^ -o $@

run: $(EXECUTABLE)
	./$(EXECUTABLE)

COVERAGE_THRESHOLD ?= 15
BRANCH_COVERAGE_THRESHOLD ?= 5

coverage-report: $(EXECUTABLE)
	./$(EXECUTABLE)
	@if [ "$(UNIT_TEST)" = "math" ]; then \
		OUT=$$(cd .. && $(GCOV) -b -c math.f90 -o obj); \
		LINE_PCT=$$(printf "%s\n" "$$OUT" | awk '/Lines executed:/{gsub("executed:","",$$2); gsub("%","",$$2); print $$2; exit}'); \
		BRANCH_PCT=$$(printf "%s\n" "$$OUT" | awk '/Branches executed:/{gsub("executed:","",$$2); gsub("%","",$$2); print $$2; exit}'); \
		test -n "$$LINE_PCT" || ( echo "FAIL: could not parse gcov line coverage output. Ensure GCOV matches mpif90 toolchain."; exit 1 ); \
		test -n "$$BRANCH_PCT" || ( echo "FAIL: could not parse gcov branch coverage output. Ensure GCOV matches mpif90 toolchain."; exit 1 ); \
		echo "math.f90 line coverage: $${LINE_PCT}%"; \
		echo "math.f90 branch coverage: $${BRANCH_PCT}%"; \
		awk -v pct="$$LINE_PCT" -v th="$(COVERAGE_THRESHOLD)" 'BEGIN{exit !(pct+0 >= th+0)}' || \
		( echo "FAIL: math.f90 line coverage $${LINE_PCT}% is below threshold $(COVERAGE_THRESHOLD)%"; exit 1 ); \
		awk -v pct="$$BRANCH_PCT" -v th="$(BRANCH_COVERAGE_THRESHOLD)" 'BEGIN{exit !(pct+0 >= th+0)}' || \
		( echo "FAIL: math.f90 branch coverage $${BRANCH_PCT}% is below threshold $(BRANCH_COVERAGE_THRESHOLD)%"; exit 1 ); \
	fi

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) *.o *.mod ut_*
