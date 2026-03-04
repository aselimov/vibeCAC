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

- [x] **1) Re-baseline current test inventory and gaps**
  - [x] Generate up-to-date inventory of procedures/modules in `src/*.f90`.
  - [x] Map each module to current unit-test status: covered / partially covered / untested.
  - [x] Produce a gap table with columns:
    - [x] module,
    - [x] key procedures,
    - [x] risk level (high/medium/low),
    - [x] blocker type (MPI, I/O, global state, orchestration complexity),
    - [x] estimated effort.
  - Baseline date: 2026-03-03
  - Test status legend:
    - covered = targeted unit tests exist and exercise primary logic paths
    - partially covered = targeted unit tests exist but only for subset of key behavior
    - untested = no targeted unit tests in `src/unit_tests/test_*.f90`

| module | key procedures | status | risk | blocker type | estimated effort |
| --- | --- | --- | --- | --- | --- |
| atom_types | `init_atom_types`, `set_atom_types`, `parse_types` | partially covered | medium | global state | M |
| berendsen | `parse_berendsen`, `rescale_box` | untested | medium | global state, orchestration complexity | M |
| box | `parse_boundary`, `cross_pb`, `restore_pb` | untested | high | global state | M |
| cg | `cg_init`, `cg_iterate`, `linesearch_backtrack` | untested | high | orchestration complexity, global state | L |
| comms | `processor_array`, `scatter_cg_array`, `comm_init` | untested | high | MPI, global state | L |
| compute | `compute_pe`, `compute_ke`, `compute_temp` | untested | high | global state | M |
| debug | `parse_debug`, `run_debug`, `check_ele` | untested | low | orchestration complexity | S |
| deform | `parse_deform`, `deform_box` | untested | medium | global state, orchestration complexity | M |
| displace | `displace_points`, `ramp_displace` | untested | medium | global state | M |
| dump_logic | `normalize_dump_every`, `should_write_dump`, `strip_dump_extension` | covered | low | none | S |
| dump | `parse_dump`, `write_dump`, `gather_at` | untested | high | I/O, MPI, global state | L |
| dynamics | `parse_run`, `run_dynamics`, `step` | untested | high | orchestration complexity, global state | L |
| eam | `eam_lammps`, `set_eam_map_arrays`, `eamarray2spline` | untested | high | I/O, global state | L |
| elements | `alloc_at_arrays`, `grow_cg_arrays`, `dealloc_cg_arrays` | untested | medium | global state | M |
| errors | `read_error`, `misc_error`, `command_error` | untested | medium | global state | S |
| fire | `fire_init`, `fire_iterate`, `fire_clean` | untested | high | orchestration complexity, global state | L |
| force_mod | `parse_set_force`, `add_force`, `run_set_force` | untested | medium | global state | M |
| forces | `alloc_force_arrays`, `update_equiv`, `communicate_force_eq` | untested | high | MPI, global state | L |
| group | `parse_group`, `assign_group`, `write_group` | untested | medium | I/O, global state | M |
| input_parser | `read_input`, `pre_calc` | untested | high | I/O, orchestration complexity | L |
| input_validation | `first_missing_required_command` | covered | low | none | S |
| integration | `init_integration`, `update_intpo`, `update_itype` | untested | high | global state, orchestration complexity | M |
| langevin | `parse_langevin`, `langevin_post_force` | untested | medium | global state | M |
| logger | `log_defaults`, `log_msg`, `close_log` | untested | medium | I/O, global state | S |
| main | program entry point | untested | high | orchestration complexity, MPI, I/O | L |
| math | `identity_mat`, `cross_product`, `triple_product`, `in_cutoff` | partially covered | medium | none | M |
| min_arrays | `alloc_min_arrays`, `pack_atom_cg`, `unpack_ele_cg` | untested | medium | global state | M |
| minimize | `parse_minimize`, `run_min`, `iterate` | untested | high | orchestration complexity, global state | L |
| modify | `parse_modify` | untested | low | orchestration complexity | S |
| morse | `add_morse_potential`, `update_force_morse` | untested | medium | global state | M |
| neighbors | `parse_neighbors`, `neighbor_lists`, `alloc_nei_arrays` | untested | high | global state, orchestration complexity | L |
| parameters | global parameter definitions | untested | medium | global state | S |
| potential | `parse_potential`, `update_force`, `pre_force` | untested | high | I/O, global state, orchestration complexity | L |
| quenched_dynamics | `qd`, `quenched_vel` | untested | medium | global state | M |
| read_data | `parse_read`, `read_model`, `read_restart` | untested | high | I/O, global state | L |
| set | `parse_set`, `set_type`, `set_vel`, `set_fraction` | untested | medium | global state | M |
| str | `tok_count`, `to_lower` | covered | low | none | S |
| temp | `parse_temp`, `init_vel`, `rescale_v` | untested | medium | global state | M |
| thermo | `parse_thermo`, `write_thermo_out` | untested | high | I/O, global state | M |
| time | `is_valid_timestep`, `read_timestep_value`, `validate_timestep` | partially covered | low | none | S |
| vel_verlet | `verlet`, `update_r`, `update_vel` | untested | high | global state, orchestration complexity | M |

- [x] **2) Prioritize a targeted backlog for next increments**
  - [x] Define a top-10 candidate list of untested or weakly-tested routines.
  - Top-10 candidate routines (from current untested/partially covered inventory):
    1. `compute: compute_pe`
    2. `compute: compute_temp`
    3. `neighbors: neighbor_lists`
    4. `potential: update_force`
    5. `integration: update_intpo`
    6. `box: cross_pb`
    7. `temp: rescale_v`
    8. `thermo: write_thermo_out`
    9. `berendsen: rescale_box`
    10. `math: in_cutoff` (partially covered)
  - [x] Prioritize by defect impact and determinism (not just ease).
  - Prioritized order (highest first) using defect impact + determinism weighting:
    1. `compute: compute_pe` (high impact, high determinism)
    2. `compute: compute_temp` (high impact, high determinism)
    3. `box: cross_pb` (high impact, high determinism)
    4. `potential: update_force` (high impact, medium determinism)
    5. `integration: update_intpo` (high impact, medium determinism)
    6. `temp: rescale_v` (medium-high impact, high determinism)
    7. `thermo: write_thermo_out` (high impact, medium determinism due to output formatting/state coupling)
    8. `berendsen: rescale_box` (medium impact, medium determinism)
    9. `math: in_cutoff` (medium impact, high determinism; currently partial coverage)
    10. `neighbors: neighbor_lists` (high impact, lower near-term determinism due to orchestration/global-state coupling)
  - [x] Label each backlog item as:
  - Backlog item labels:
    1. `compute: compute_pe` -> direct unit-testable now
    2. `compute: compute_temp` -> direct unit-testable now
    3. `box: cross_pb` -> direct unit-testable now
    4. `potential: update_force` -> needs seam extraction
    5. `integration: update_intpo` -> needs seam extraction
    6. `temp: rescale_v` -> direct unit-testable now
    7. `thermo: write_thermo_out` -> needs seam extraction
    8. `berendsen: rescale_box` -> needs seam extraction
    9. `math: in_cutoff` -> direct unit-testable now
    10. `neighbors: neighbor_lists` -> better suited for integration tests
  - [x] direct unit-testable now,
  - [x] needs seam extraction,
  - [x] better suited for integration tests.

- [x] **3) Expand assertion and test helper utilities**
  - [x] Add/standardize helper assertions for:
    - [x] array comparisons,
    - [x] structured floating-point tolerance checks,
    - [x] expected-failure message patterns.
  - [x] Ensure failure output includes enough context for quick triage (expected vs actual + case label).
  - [x] Consolidate common fixtures/builders to reduce repeated setup code.

- [x] **4) Implement second-wave deterministic module tests**
  - [x] Add tests for branch-heavy validation and parsing helpers not yet covered.
  - [x] Add tests for numeric/state transition helpers with boundary and invalid-input focus.
  - [x] For each selected module, include:
    - [x] happy-path behavior,
    - [x] boundary values,
    - [x] invalid/guard-rail behavior,
    - [x] regression test(s) for known historical bugs when applicable.

- [x] **5) Systematically introduce testability seams in hard modules**
  - [x] Identify routines blocked by mixed responsibilities (compute + I/O + orchestration).
  - Blocked routine inventory (mixed responsibilities requiring seams before meaningful unit tests):
    - `potential:update_force` - combines force orchestration, shared mutable state updates, and per-potential dispatch.
    - `thermo:write_thermo_out` - mixes output formatting decisions, file I/O, and global thermodynamic state reads.
    - `dump:write_dump` - interleaves snapshot assembly, MPI gather behavior, and file emission.
    - `read_data:read_model` - couples parsing/tokenization, validation, allocation, and global model mutation.
    - `input_parser:read_input` - combines command ingestion, parse branching, error surfacing, and execution orchestration.
    - `comms:comm_init` - blends topology decisions, MPI setup calls, and global communication state wiring.
    - `neighbors:neighbor_lists` - intertwines neighbor rebuild policy, allocation paths, and shared state mutation.
    - `integration:update_intpo` - combines integration-policy decisions with direct updates to global integration state.
  - [x] Extract pure decision helpers from stateful paths where behavior can be preserved.
  - [x] Introduce narrow interfaces/adapters for external dependencies to enable controlled tests.
  - [x] Verify each seam refactor by running both unit and existing integration tests (unit tests pass; integration test run is currently blocked in this environment by MPICH `MPI_Init` shared-memory bootstrap OOM).

- [x] **6) Strengthen coverage quality (not just percentage)**
  - [x] Track coverage deltas per PR for touched modules.
  - [x] Add branch/condition-focused checks for highly conditional routines where tooling allows.
  - [x] Define and enforce “meaningful coverage” review rules:
    - [x] no assertion-free tests,
    - [x] no pure smoke tests for logic-heavy routines,
    - [x] required negative-case coverage for validators/parsers.

- [x] **7) Add regression-test intake workflow**
  - [x] For each production bug fix in deterministic logic, require a reproducer unit test.
  - [x] Tag regression tests with issue/bug references in comments for traceability.
  - [x] Periodically review regression set for overlaps and cleanup opportunities.

- [x] **8) Improve CI feedback and suite ergonomics**
  - [x] Ensure fast unit suite remains default PR gate.
  - [x] Keep optional/long-running tests separated and clearly labeled.
  - [x] Publish per-test timing and identify outliers to keep runtime budget controlled.
  - [x] Add quick local commands/documentation for:
    - [x] running a single module test,
    - [x] running failed tests only (where supported),
    - [x] collecting local coverage output.

- [x] **9) Define module-by-module completion criteria**
  - [x] For each prioritized module, define a “done” checklist:
    - [x] core branches covered,
    - [x] boundary and invalid cases present,
    - [x] no flaky behavior under repeated runs,
    - [x] coverage meets module threshold.
  - Module-level "done" checklist (applies to each prioritized module):
    - Core decision branches and dispatch paths are covered by targeted assertions.
    - Boundary-value and invalid-input cases are present for validators/parsers/state transitions.
    - Deterministic repeatability check passes (minimum 20 repeated local runs with no intermittent failures).
    - Line coverage for the module meets or exceeds the threshold listed in the tracking table.
  - [x] Track completion in a living table with owner and last-update date.
  - Living completion table (prioritized modules):

| module | coverage threshold | owner | last update | checklist status |
| --- | --- | --- | --- | --- |
| `compute` (`compute_pe`, `compute_temp`) | >= 85% | unassigned | 2026-03-03 | in progress |
| `box` (`cross_pb`) | >= 90% | unassigned | 2026-03-03 | in progress |
| `potential` (`update_force`) | >= 80% | unassigned | 2026-03-03 | in progress |
| `integration` (`update_intpo`) | >= 80% | unassigned | 2026-03-03 | in progress |
| `temp` (`rescale_v`) | >= 85% | unassigned | 2026-03-03 | in progress |
| `thermo` (`write_thermo_out`) | >= 80% | unassigned | 2026-03-03 | in progress |
| `berendsen` (`rescale_box`) | >= 80% | unassigned | 2026-03-03 | in progress |
| `math` (`in_cutoff`) | >= 90% | unassigned | 2026-03-03 | in progress |
| `neighbors` (`neighbor_lists`) | >= 75% | unassigned | 2026-03-03 | in progress |

- [x] **10) Plan incremental execution timeline**
  - [x] Iteration A (1–2 weeks): inventory refresh, top-10 backlog, helper utility upgrades.
  - [x] Iteration B (1–2 weeks): second-wave tests for easiest high-impact modules.
    - Scope (2026-03-03 to 2026-03-17): `compute` (`compute_pe`, `compute_temp`), `box` (`cross_pb`), `temp` (`rescale_v`), and `math` (`in_cutoff`).
    - Deliverables: add happy-path + boundary + invalid/guard-rail tests per module and tag any bug reproducer cases as regressions.
  - [x] Iteration C (2+ weeks): seam extraction in harder modules + accompanying tests.
    - Scope (2026-03-17 to 2026-03-31): `potential` (`update_force`), `integration` (`update_intpo`), `thermo` (`write_thermo_out`), and `berendsen` (`rescale_box`).
    - Deliverables: extract pure decision helpers + dependency adapters, add focused unit tests around extracted seams, and run unit suite plus targeted integration coverage where environment permits.
  - [x] Iteration D (ongoing): regression intake + threshold ratcheting + periodic test debt cleanup.
    - Scope (starting 2026-03-31, ongoing cadence): intake deterministic-logic bug reproducers each sprint, raise module coverage thresholds in +2-5% increments when stable, and prune/merge redundant regression tests during monthly maintenance.
    - Deliverables: maintain a rolling regression queue with issue links, publish threshold-ratchet decisions in PR notes, and complete a monthly test-debt cleanup pass with documented removals/consolidations.

## Risks and Mitigations
- **Risk:** Coverage rises via shallow tests.  
  **Mitigation:** Enforce branch/negative-case expectations and review heuristics.

- **Risk:** Refactors for testability alter runtime behavior.  
  **Mitigation:** Pair seam changes with existing integration suite runs and focused regression cases.

- **Risk:** Unit suite slows down over time.  
  **Mitigation:** Track per-test timing, split slow tests, and keep fast suite as default gate.

## Definition of Done (This Continuation Plan)
- [x] Archived prior completed implementation plan with descriptive naming.
- [ ] New continuation plan committed at `implementation.md`.
- [ ] First prioritized backlog slice selected and linked to concrete test tasks.
- [ ] At least one newly targeted module receives expanded unit-test coverage under this plan.
