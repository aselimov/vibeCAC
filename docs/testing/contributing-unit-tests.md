# Contributing Unit Tests

This guide describes how to write, run, and maintain fast, deterministic unit tests under `src/unit_tests/`.

## Local workflow

From `src/`:

1. Run a specific suite:
   - `make unit-test UNIT_TEST=math`
   - `make unit-test UNIT_TEST=str`
2. Run coverage gate for the baseline deterministic target:
   - `make unit-test-coverage UNIT_TEST=math`

## Quick local commands

From `src/`:

1. Run a single module/unit target:
   - `make unit-test UNIT_TEST=<target>`
2. Re-run failed tests only (where supported):
   - `ctest --rerun-failed --output-on-failure`
   - Note: the default Make-based unit harness does not track prior failures; use `UNIT_TEST=<target>` to re-run a specific failing target.
3. Collect local coverage output:
   - `make unit-test-coverage UNIT_TEST=math`
   - `gcov-15 -b -c math.f90 -o obj` (or `gcov -b -c math.f90 -o obj` if `gcov-15` is unavailable)

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
  - `make unit-test-coverage UNIT_TEST=math COVERAGE_THRESHOLD=<value> BRANCH_COVERAGE_THRESHOLD=<value>`

## Flaky-test policy

Unit tests must be deterministic and should not rely on wall-clock timing, random seeds without control, file-system race conditions, or MPI rank ordering side effects.

If a flaky test is found:

1. Mark it as quarantined by commenting out its invocation in the test driver and add a `TODO(flaky)` note with issue/owner.
2. Open a tracking issue immediately.
3. Fix SLA:
   - critical paths: within 2 business days
   - non-critical paths: within 5 business days
4. Re-enable only after 10 consecutive local passes.

## Meaningful coverage review rules

- Assertion-free unit tests are disallowed. Every `src/unit_tests/test_*.f90` file must include at least one `call assert_*` check.
- Pure smoke tests are disallowed for logic-heavy routines (`compute`, `potential`, `integration`, `box`, `temp`, `neighbors`, `thermo`, `berendsen`, `math`). Tests for these targets must include at least one non-tautological behavioral assertion.
- Validator/parser targets must include at least one explicit negative-case assertion (invalid/rejected/missing/error path).
- CI enforces this with `.github/scripts/enforce_unit_assertions.sh`.

## Regression-test intake workflow

- Every production bug fix in deterministic logic must include a reproducer unit test in `src/unit_tests/test_*.f90`.
- Each regression test must include a traceability comment directly above the test body using:
  - `! Regression: BUGREF(<issue-or-bug-id>) <short summary>`
- If no external issue tracker ID exists, use a stable local reference:
  - `BUGREF(local:YYYYMMDD-<slug>)`
- Review cadence: once per release cycle, run:
  - `rg -n "Regression: BUGREF\(" src/unit_tests/test_*.f90`
- During review, merge duplicate reproducers and remove obsolete cases only when equivalent coverage remains.

## Ownership and rotation

- Primary owner: maintainers touching `src/unit_tests/` and `src/Makefile` in the same PR.
- Rotation: review responsibility rotates each release cycle among active maintainers.
- Any change to deterministic helper logic in `src/*.f90` should include unit tests or an explicit PR exemption rationale.
