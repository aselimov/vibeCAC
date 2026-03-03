.SUFFIXES:

FC = mpif90
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

ifeq ($(UNIT_TEST),dump_logic)
CORE_OBJECTS += dump_logic.o
endif

ifeq ($(UNIT_TEST),input_validation)
CORE_OBJECTS += input_validation.o
endif

TEST_OBJECTS = assertions.o test_$(UNIT_TEST).o
EXECUTABLE = ut_$(UNIT_TEST)

.PHONY: all run coverage-report clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(CORE_OBJECTS) $(TEST_OBJECTS)
	$(FC) $(FFLAGS) $^ -o $@

run: $(EXECUTABLE)
	./$(EXECUTABLE)

COVERAGE_THRESHOLD ?= 5

coverage-report: $(EXECUTABLE)
	./$(EXECUTABLE)
	@if [ "$(UNIT_TEST)" = "math" ]; then \
		OUT=$$(cd .. && gcov -b -c math.f90 -o obj); \
		PCT=$$(printf "%s\n" "$$OUT" | awk '/Lines executed:/{gsub("executed:","",$$2); gsub("%","",$$2); print $$2; exit}'); \
		echo "math.f90 line coverage: $${PCT}%"; \
		awk -v pct="$$PCT" -v th="$(COVERAGE_THRESHOLD)" 'BEGIN{exit !(pct+0 >= th+0)}' || \
		( echo "FAIL: math.f90 coverage $${PCT}% is below threshold $(COVERAGE_THRESHOLD)%"; exit 1 ); \
	fi

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) *.o *.mod ut_*
