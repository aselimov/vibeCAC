# Function-Level Unit Testing Implementation Plan

## Objective
Build a repeatable, automated unit-testing workflow for individual Fortran procedures/modules in `src/`, while preserving existing scenario/integration tests under `src/tests/`.

## Scope and Constraints
- Focus on **function-level/module-level** verification (small, deterministic tests).
- Keep compatibility with existing compiler/toolchain assumptions (`mpif90`, GNU Fortran flags in current makefiles).
- Minimize disruption to production source layout.

## Step-by-Step Requirements

- [ ] **Inventory and classify testable units**
  - [ ] Enumerate candidate procedures in `src/*.f90` by module/file.
  - [ ] Tag each as: pure math/helper, stateful logic, MPI-dependent, file-I/O-dependent.
  - [ ] Prioritize units that are deterministic and side-effect light (e.g., math, parser/token logic, utility transforms).

- [ ] **Define unit-test architecture and conventions**
  - [ ] Create a dedicated unit-test tree (e.g., `src/unit_tests/`) separate from current integration-style tests in `src/tests/`.
  - [ ] Establish naming rules:
    - [ ] test files: `test_<module_or_function>.f90`
    - [ ] test executables: `ut_<module_or_function>`
  - [ ] Standardize assertion style (`assert_equal`, `assert_close`, `assert_true`, `assert_false`) and failure messaging format.
  - [ ] Require one test case = one clearly named procedure/subroutine with explicit Arrange/Act/Assert sections.

- [ ] **Choose and integrate a Fortran unit testing framework**
  - [ ] Evaluate `pFUnit` vs. `test-drive` (or equivalent lightweight harness) against project constraints:
    - [ ] compiler support,
    - [ ] MPI compatibility needs,
    - [ ] setup complexity,
    - [ ] CI friendliness.
  - [ ] Select one framework and document rationale in a short `docs/testing/unit-testing.md` design note.
  - [ ] Add bootstrap/install instructions for local developer environments.

- [ ] **Create build system support for unit tests**
  - [ ] Add a top-level build target for unit tests (e.g., `make unit-test`).
  - [ ] Add compile/link rules for unit-test binaries that link required production modules.
  - [ ] Ensure unit-test builds use debug-friendly flags (`-O0 -g -fcheck=all -fbacktrace` or equivalent).
  - [ ] Add selective execution support (single test file, single test case where framework permits).

- [ ] **Improve testability seams in production code**
  - [ ] Refactor tightly coupled routines to allow injection/mocking of external dependencies where practical (I/O, global state, MPI wrappers).
  - [ ] Isolate pure computation from orchestration code into smaller procedures.
  - [ ] Add explicit interfaces and reduce hidden dependencies to simplify direct procedure invocation in tests.
  - [ ] Preserve runtime behavior by validating refactors against existing `src/tests` workloads.

- [ ] **Implement first wave of baseline unit tests (high-value targets)**
  - [ ] Add initial tests for deterministic numeric helpers (e.g., in `math.f90`, related utility modules).
  - [ ] Add tests for parser/input validation edge cases where deterministic.
  - [ ] Add boundary-condition tests:
    - [ ] zero/empty inputs,
    - [ ] min/max valid ranges,
    - [ ] invalid argument handling.
  - [ ] Add floating-point tolerance policy and verify tests use consistent epsilons.

- [ ] **Define coverage and quality gates**
  - [ ] Add line/function coverage reporting workflow (e.g., `gcov`/`lcov` compatible with Fortran toolchain).
  - [ ] Set initial minimum coverage targets for unit-test scope (starting threshold, then ratchet policy).
  - [ ] Add pass/fail gating rules for new/modified code (must include unit tests or explicit exemption note).

- [ ] **Integrate unit tests into CI and developer workflow**
  - [ ] Add CI job(s) to compile and run unit tests on every PR.
  - [ ] Add separate CI stage for existing heavier integration tests to keep signal clear.
  - [ ] Enforce fast-feedback runtime target for unit suite (goal: a few minutes).
  - [ ] Publish test and coverage artifacts in CI for debugging failures.

- [ ] **Document developer usage and maintenance policy**
  - [ ] Add a contributor guide section: how to write, run, and debug unit tests locally.
  - [ ] Provide templates/examples for adding a new unit test file.
  - [ ] Define policy for flaky tests (identification, quarantine, fix SLA).
  - [ ] Define ownership/rotation for maintaining test infrastructure.

- [ ] **Rollout plan**
  - [ ] Phase 1: framework + build plumbing + 5–10 representative unit tests.
  - [ ] Phase 2: cover core deterministic modules and enforce test requirement for touched functions.
  - [ ] Phase 3: extend seams/mocks for currently hard-to-test logic and raise coverage gates.
  - [ ] Conduct retrospective after Phase 1 and adjust conventions/tooling before scaling.

## Done Criteria
- [ ] A documented and reproducible unit-test framework is present.
- [ ] `make unit-test` (or equivalent) runs green locally and in CI.
- [ ] Initial high-value function-level tests are merged and stable.
- [ ] Contribution workflow requires tests (or justified exemption) for functional changes.
