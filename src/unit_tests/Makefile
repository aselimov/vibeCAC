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

TEST_OBJECTS = assertions.o test_$(UNIT_TEST).o
EXECUTABLE = ut_$(UNIT_TEST)

.PHONY: all run clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(CORE_OBJECTS) $(TEST_OBJECTS)
	$(FC) $(FFLAGS) $^ -o $@

run: $(EXECUTABLE)
	./$(EXECUTABLE)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) *.o *.mod ut_*
