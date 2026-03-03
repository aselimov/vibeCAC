.PHONY: unit-test unit-test-coverage tests clean

unit-test:
	$(MAKE) -C src unit-test UNIT_TEST=$(UNIT_TEST) UNIT_TEST_FFLAGS="$(UNIT_TEST_FFLAGS)"

unit-test-coverage:
	$(MAKE) -C src unit-test-coverage UNIT_TEST=$(UNIT_TEST)

tests:
	$(MAKE) -C src tests

clean:
	$(MAKE) -C src clean
	$(MAKE) -C src/unit_tests clean
