# Contributing Unit Tests

This guide describes how to write, run, and maintain fast, deterministic unit tests under `src/unit_tests/`.

## Local workflow

From `src/`:

1. Run a specific suite:
   - `make unit-test UNIT_TEST=math`
   - `make unit-test UNIT_TEST=str`
2. Run coverage gate for the baseline deterministic target:
   - `make unit-test-coverage UNIT_TEST=math`

## Test file template

Use one file per module/function target:

- file: `src/unit_tests/test_<target>.f90`
- executable: `ut_<target>`

Skeleton:

```fortran
program test_example
    use parameters
    use assertions
    implicit none

    call test_case_name()
    write(*,*) "All example unit tests passed."

contains

    subroutine test_case_name()
        ! Arrange
        ! ... prepare inputs

        ! Act
        ! ... call target procedure

        ! Assert
        call assert_true(.true., "FAIL: expected behavior")
    end subroutine test_case_name

end program test_example
```

## Debugging failures

- Keep debug flags enabled (`-O0 -g -fcheck=all -fbacktrace`).
- Re-run a single target repeatedly:
  - `make unit-test UNIT_TEST=<target>`
- For coverage-related failures, run:
  - `make unit-test-coverage UNIT_TEST=math COVERAGE_THRESHOLD=<value>`

## Flaky-test policy

Unit tests must be deterministic and should not rely on wall-clock timing, random seeds without control, file-system race conditions, or MPI rank ordering side effects.

If a flaky test is found:

1. Mark it as quarantined by commenting out its invocation in the test driver and add a `TODO(flaky)` note with issue/owner.
2. Open a tracking issue immediately.
3. Fix SLA:
   - critical paths: within 2 business days
   - non-critical paths: within 5 business days
4. Re-enable only after 10 consecutive local passes.

## Ownership and rotation

- Primary owner: maintainers touching `src/unit_tests/` and `src/Makefile` in the same PR.
- Rotation: review responsibility rotates each release cycle among active maintainers.
- Any change to deterministic helper logic in `src/*.f90` should include unit tests or an explicit PR exemption rationale.
