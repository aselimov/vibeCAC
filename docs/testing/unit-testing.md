# Unit Testing Design Note

This document defines a minimal, repeatable unit-testing workflow for deterministic Fortran procedures in `src/`.

## 1) Initial inventory and classification

The first wave targets deterministic helpers with low side effects.

| File/module | Candidate procedure(s) | Classification | Priority |
| --- | --- | --- | --- |
| `src/math.f90` | `identity_mat`, `cross_product`, `triple_product` | pure math/helper | High |
| `src/math.f90` | `matrix_normal`, `sort_array` | stateful logic (deterministic) | Medium |
| `src/input_parser.f90` | parser/token helpers | parser/input validation | Medium |
| `src/comms.f90` | MPI communication routines | MPI-dependent | Low (unit tests deferred) |
| `src/read_data.f90` | file parsing and import | file-I/O-dependent | Low (requires test seams) |

## 2) Unit-test architecture and conventions

- Unit tests live in `src/unit_tests/` (separate from integration tests in `src/tests/`).
- File naming:
  - test source: `test_<module_or_function>.f90`
  - executable: `ut_<module_or_function>`
- Test case structure: one clearly named procedure per test case with explicit `Arrange / Act / Assert` comments.
- Assertion helpers are provided in `src/unit_tests/assertions.f90`:
  - `assert_true`
  - `assert_false`
  - `assert_equal_int`
  - `assert_close_real`
- Failure messages follow `FAIL: <test intent>` and include actual/expected where applicable.

### Floating-point tolerance policy

- Default tolerance for deterministic math helpers: `1.0e-12_wp`.
- If a unit needs looser bounds, document the reason inline with the assertion.

## 3) Framework evaluation and selection

### Candidates

- **pFUnit**
  - Pros: mature xUnit feature set, parameterized tests, MPI support.
  - Cons: heavier setup and toolchain integration overhead for this repository.
- **test-drive (lightweight Fortran harness style)**
  - Pros: minimal ceremony, simple compiler integration, CI-friendly for fast deterministic tests.
  - Cons: fewer advanced testing features than pFUnit.

### Selection

For Phase 1, we adopt a **test-drive style lightweight harness** using local assertion utilities and simple test executables in `src/unit_tests/`. This keeps setup small and compatible with the current `mpif90`-based workflow while still giving deterministic function-level coverage.

## 4) Build and execution

From `src/`:

- Run the default unit test target:
  - `make unit-test`
- Run a specific unit test file/executable:
  - `make unit-test UNIT_TEST=math`

Both paths use debug-friendly flags:

- `-O0 -g -fcheck=all -fbacktrace`

## 5) Bootstrap for local developers

- Ensure an MPI Fortran compiler wrapper is available (`mpif90`).
- Build/run from the `src/` directory:
  1. `make unit-test`
  2. Optionally target a specific test using `UNIT_TEST=<name>`.

No external framework installation is required for the current Phase 1 harness.

## 6) Coverage reporting and quality gates

From `src/`, run:

- `make unit-test-coverage`

This compiles production and unit-test objects with `--coverage`, executes the selected unit test (default `UNIT_TEST=math`), then emits `gcov` reports.

### Initial coverage target

- Initial line-coverage threshold for the first deterministic target (`math.f90`) is **15%**.
- The threshold is enforced by `src/unit_tests/Makefile` and can be adjusted via:
  - `make unit-test-coverage COVERAGE_THRESHOLD=<percent>`
- Coverage is generated with a GCC-compatible `gcov` binary (for example `gcov-15`) to match the `mpif90`/GNU Fortran toolchain.

### Gate policy for touched code

- New or modified deterministic helper logic should include/extend a corresponding unit test under `src/unit_tests/`.
- If unit tests are not practical (for MPI-heavy or I/O-heavy changes), add an explicit exemption note in the PR description explaining why and what broader test validates the behavior.

## 7) Contributor guide and maintenance policy

See `docs/testing/contributing-unit-tests.md` for:

- how to write and run unit tests locally,
- a template for adding a new unit-test file,
- flaky-test quarantine and fix SLA policy,
- ownership and review rotation guidance for test infrastructure.

## 8) Phase 1 retrospective (2026-03-03)

Phase 1 goals (framework + build plumbing + representative tests) were met with the current lightweight harness and `make unit-test` integration.

Retrospective outcomes and adjustments before scale-up:

- Keep the local assertion harness for deterministic units; defer pFUnit adoption unless MPI-heavy unit coverage becomes a priority.
- Keep `UNIT_TEST=<name>` selective execution as the default debug path and require new unit targets to support it.
- Preserve the initial `math.f90` coverage gate baseline and ratchet only when adding deterministic tests, not during unrelated infrastructure changes.
- Continue extracting pure decision helpers from stateful modules before adding tests, to avoid brittle tests tied to global state.
