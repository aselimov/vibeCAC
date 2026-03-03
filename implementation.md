# Function-Level Unit Testing Continuation Plan (Phase 2+)

## Objective
Continue expanding and hardening unit-test coverage across deterministic and semi-isolated logic in `src/`, with emphasis on maintainability, defect-prevention, and incremental testability improvements for harder-to-isolate routines.

## Context
- A baseline unit-test framework and CI workflow already exist.
- Initial seams were extracted and first-wave tests are in place.
- The next phase should focus on sustainable scale: broader coverage, higher signal quality, and explicit tracking of untested/high-risk areas.

## Success Metrics (Rolling)
- Unit-test count grows steadily each sprint (target: +10 to +20 meaningful tests per cycle where feasible).
- New/refactored deterministic logic lands with unit tests by default.
- Unit test failures are actionable (low flake rate, clear assertions, minimal false positives).
- Coverage increases in prioritized modules without adding brittle tests.

## Guiding Principles
- Prefer **small, deterministic, behavior-focused** tests over implementation-coupled tests.
- Test **decision logic first** (branch-heavy validation, boundary checks, parser behavior).
- Use seams and helper extraction to isolate logic from I/O, MPI, and global-state coupling.
- Keep test runtime fast and predictable; split heavier cases to dedicated suites when needed.

## Step-by-Step Continuation Requirements

- [ ] **1) Re-baseline current test inventory and gaps**
  - [ ] Generate up-to-date inventory of procedures/modules in `src/*.f90`.
  - [ ] Map each module to current unit-test status: covered / partially covered / untested.
  - [ ] Produce a gap table with columns:
    - [ ] module,
    - [ ] key procedures,
    - [ ] risk level (high/medium/low),
    - [ ] blocker type (MPI, I/O, global state, orchestration complexity),
    - [ ] estimated effort.

- [ ] **2) Prioritize a targeted backlog for next increments**
  - [ ] Define a top-10 candidate list of untested or weakly-tested routines.
  - [ ] Prioritize by defect impact and determinism (not just ease).
  - [ ] Label each backlog item as:
    - [ ] direct unit-testable now,
    - [ ] needs seam extraction,
    - [ ] better suited for integration tests.

- [ ] **3) Expand assertion and test helper utilities**
  - [ ] Add/standardize helper assertions for:
    - [ ] array comparisons,
    - [ ] structured floating-point tolerance checks,
    - [ ] expected-failure message patterns.
  - [ ] Ensure failure output includes enough context for quick triage (expected vs actual + case label).
  - [ ] Consolidate common fixtures/builders to reduce repeated setup code.

- [ ] **4) Implement second-wave deterministic module tests**
  - [ ] Add tests for branch-heavy validation and parsing helpers not yet covered.
  - [ ] Add tests for numeric/state transition helpers with boundary and invalid-input focus.
  - [ ] For each selected module, include:
    - [ ] happy-path behavior,
    - [ ] boundary values,
    - [ ] invalid/guard-rail behavior,
    - [ ] regression test(s) for known historical bugs when applicable.

- [ ] **5) Systematically introduce testability seams in hard modules**
  - [ ] Identify routines blocked by mixed responsibilities (compute + I/O + orchestration).
  - [ ] Extract pure decision helpers from stateful paths where behavior can be preserved.
  - [ ] Introduce narrow interfaces/adapters for external dependencies to enable controlled tests.
  - [ ] Verify each seam refactor by running both unit and existing integration tests.

- [ ] **6) Strengthen coverage quality (not just percentage)**
  - [ ] Track coverage deltas per PR for touched modules.
  - [ ] Add branch/condition-focused checks for highly conditional routines where tooling allows.
  - [ ] Define and enforce “meaningful coverage” review rules:
    - [ ] no assertion-free tests,
    - [ ] no pure smoke tests for logic-heavy routines,
    - [ ] required negative-case coverage for validators/parsers.

- [ ] **7) Add regression-test intake workflow**
  - [ ] For each production bug fix in deterministic logic, require a reproducer unit test.
  - [ ] Tag regression tests with issue/bug references in comments for traceability.
  - [ ] Periodically review regression set for overlaps and cleanup opportunities.

- [ ] **8) Improve CI feedback and suite ergonomics**
  - [ ] Ensure fast unit suite remains default PR gate.
  - [ ] Keep optional/long-running tests separated and clearly labeled.
  - [ ] Publish per-test timing and identify outliers to keep runtime budget controlled.
  - [ ] Add quick local commands/documentation for:
    - [ ] running a single module test,
    - [ ] running failed tests only (where supported),
    - [ ] collecting local coverage output.

- [ ] **9) Define module-by-module completion criteria**
  - [ ] For each prioritized module, define a “done” checklist:
    - [ ] core branches covered,
    - [ ] boundary and invalid cases present,
    - [ ] no flaky behavior under repeated runs,
    - [ ] coverage meets module threshold.
  - [ ] Track completion in a living table with owner and last-update date.

- [ ] **10) Plan incremental execution timeline**
  - [ ] Iteration A (1–2 weeks): inventory refresh, top-10 backlog, helper utility upgrades.
  - [ ] Iteration B (1–2 weeks): second-wave tests for easiest high-impact modules.
  - [ ] Iteration C (2+ weeks): seam extraction in harder modules + accompanying tests.
  - [ ] Iteration D (ongoing): regression intake + threshold ratcheting + periodic test debt cleanup.

## Risks and Mitigations
- **Risk:** Coverage rises via shallow tests.  
  **Mitigation:** Enforce branch/negative-case expectations and review heuristics.

- **Risk:** Refactors for testability alter runtime behavior.  
  **Mitigation:** Pair seam changes with existing integration suite runs and focused regression cases.

- **Risk:** Unit suite slows down over time.  
  **Mitigation:** Track per-test timing, split slow tests, and keep fast suite as default gate.

## Definition of Done (This Continuation Plan)
- [ ] Archived prior completed implementation plan with descriptive naming.
- [ ] New continuation plan committed at `implementation.md`.
- [ ] First prioritized backlog slice selected and linked to concrete test tasks.
- [ ] At least one newly targeted module receives expanded unit-test coverage under this plan.
