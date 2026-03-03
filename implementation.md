# Function-Level Unit Testing Implementation Plan

## Objective
Build a repeatable, automated unit-testing workflow for individual Fortran procedures/modules in `src/`, while preserving existing scenario/integration tests under `src/tests/`.

## Scope and Constraints
- Focus on **function-level/module-level** verification (small, deterministic tests).
- Keep compatibility with existing compiler/toolchain assumptions (`mpif90`, GNU Fortran flags in current makefiles).
- Minimize disruption to production source layout.

## Step-by-Step Requirements

- [x] **Inventory and classify testable units**
  - [x] Enumerate candidate procedures in `src/*.f90` by module/file.
  - [x] Tag each as: pure math/helper, stateful logic, MPI-dependent, file-I/O-dependent.
  - [x] Prioritize units that are deterministic and side-effect light (e.g., math, parser/token logic, utility transforms).

- [x] **Define unit-test architecture and conventions**
  - [x] Create a dedicated unit-test tree (e.g., `src/unit_tests/`) separate from current integration-style tests in `src/tests/`
  - [x] Establish naming rules:
    - [x] test files: `test_<module_or_function>.f90`
    - [x] test executables: `ut_<module_or_function>`
  - [x] Standardize assertion style (`assert_equal`, `assert_close`, `assert_true`, `assert_false`) and failure messaging format.
  - [x] Require one test case = one clearly named procedure/subroutine with explicit Arrange/Act/Assert sections.

- [x] **Choose and integrate a Fortran unit testing framework**
  - [x] Evaluate `pFUnit` vs. `test-drive` (or equivalent lightweight harness) against project constraints:
    - [x] compiler support,
    - [x] MPI compatibility needs,
    - [x] setup complexity,
    - [x] CI friendliness.
  - [x] Select one framework and document rationale in a short `docs/testing/unit-testing.md` design note.
  - [x] Add bootstrap/install instructions for local developer environments.

- [x] **Create build system support for unit tests**
  - [x] Add a top-level build target for unit tests (e.g., `make unit-test`).
  - [x] Add compile/link rules for unit-test binaries that link required production modules.
  - [x] Ensure unit-test builds use debug-friendly flags (`-O0 -g -fcheck=all -fbacktrace` or equivalent).
  - [x] Add selective execution support (single test file, single test case where framework permits).

- [ ] **Improve testability seams in production code**
  - [ ] Refactor tightly coupled routines to allow injection/mocking of external dependencies where practical (I/O, global state, MPI wrappers).
  - [ ] Isolate pure computation from orchestration code into smaller procedures.
  - [ ] Add explicit interfaces and reduce hidden dependencies to simplify direct procedure invocation in tests.
  - [ ] Preserve runtime behavior by validating refactors against existing `src/tests` workloads.

- [x] **Implement first wave of baseline unit tests (high-value targets)**
  - [x] Add initial tests for deterministic numeric helpers (e.g., in `math.f90`, related utility modules).
  - [x] Add tests for parser/input validation edge cases where deterministic.
  - [x] Add boundary-condition tests:
    - [x] zero/empty inputs,
    - [x] min/max valid ranges,
    - [x] invalid argument handling.
  - [x] Add floating-point tolerance policy and verify tests use consistent epsilons.

- [x] **Define coverage and quality gates**
  - [x] Add line/function coverage reporting workflow (e.g., `gcov`/`lcov` compatible with Fortran toolchain).
  - [x] Set initial minimum coverage targets for unit-test scope (starting threshold, then ratchet policy).

- [x] **Integrate unit tests into CI and developer workflow**
  - [x] Add CI job(s) to compile and run unit tests on every PR.
  - [x] Add separate CI stage for existing heavier integration tests to keep signal clear.
  - [x] Enforce fast-feedback runtime target for unit suite (goal: a few minutes).
  - [x] Publish test and coverage artifacts in CI for debugging failures.

- [x] **Document developer usage and maintenance policy**
  - [x] Provide templates/examples for adding a new unit test file.
  - [x] Define policy for flaky tests (identification, quarantine, fix SLA).
  - [x] Define ownership/rotation for maintaining test infrastructure.

- [ ] **Rollout plan**
  - [x] Phase 1: framework + build plumbing + 5–10 representative unit tests.
  - [ ] Phase 2: cover core deterministic modules.
  - [ ] Phase 3: extend seams/mocks for currently hard-to-test logic and raise coverage gates.
  - [ ] Conduct retrospective after Phase 1 and adjust conventions/tooling before scaling.

## Done Criteria
- [x] A documented and reproducible unit-test framework is present.
- [x] `make unit-test` (or equivalent) runs green locally and in CI.
- [x] Initial high-value function-level tests are merged and stable.
