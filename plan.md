# C++ Conversion Plan for `vibeCAC`

## 1. Executive Summary

This repository is a medium-sized scientific simulation codebase written primarily in Fortran 90, built with `mpif90`, organized around global-state modules, and split across three major concerns:

1. **Simulation core**: MPI domain decomposition, neighbor lists, force evaluation, integration, minimization, and thermodynamic output.
2. **Input/runtime orchestration**: command parsing, defaults initialization, simulation loop control, and data loading.
3. **Verification and docs**: lightweight unit tests, physics/integration tests, and MkDocs documentation.

A successful conversion to C++ should **not** be treated as a literal line-by-line language translation. The Fortran implementation relies heavily on:

- shared mutable module state,
- implicit coupling through `use parameters` and other central modules,
- MPI wrappers and domain decomposition logic,
- array-oriented numerics,
- Makefile-based build/test workflows.

Because of that, the recommended path is a **phased architectural migration** that preserves physics correctness and command compatibility while gradually replacing Fortran subsystems with C++ equivalents. The best target end state is a **modern C++17 or C++20** codebase with:

- **CMake** as the build system,
- **MPI** retained through standard C MPI bindings or a thin wrapper,
- optional BLAS/LAPACK only where demonstrably valuable,
- explicit ownership and data flow,
- unit-tested pure logic extracted from orchestration-heavy runtime code,
- compatibility layers for existing input files and regression tests.

## 2. Current-State Assessment

### 2.1 High-level structure

The repository centers on `src/*.f90`, with a single `main` program initializing MPI, defaults, logging, input parsing, and runtime execution. The top-level build delegates to `src/Makefile`, which compiles all Fortran sources using `mpif90` and links a single `CAC` executable. There is also a lightweight unit-test harness under `src/unit_tests/` and broader simulation/integration tests under `src/tests/`.

### 2.2 Architecture characteristics that matter for migration

#### A. Global-state-heavy design

The code strongly depends on module-level state, especially via `parameters`, but also through modules such as `elements`, `forces`, `neighbors`, `integration`, and `comms`. This is the single biggest migration driver because C++ will work best if this state is reorganized into explicit data structures and services.

#### B. MPI is a first-class runtime dependency

MPI is used in startup, communication, domain decomposition, error handling, timing, and portions of compute/data distribution. That means the conversion plan must preserve a distributed-memory execution model rather than first converting to serial code.

#### C. Large “core engine” modules exist

A few files are especially large and likely to dominate migration cost:

- `comms.f90`
- `neighbors.f90`
- `eam.f90`
- `dump.f90`
- `elements.f90`
- `math.f90`
- `cg.f90`
- `integration.f90`
- `temp.f90`

These should be treated as independent migration workstreams rather than bundled together.

#### D. Two test layers already exist

The repository already distinguishes between:

- **unit tests** for deterministic helpers, and
- **integration/physics tests** for simulation behavior.

That is an advantage: the migration can preserve existing expected behavior instead of relying on informal manual checks.

### 2.3 Likely domain model from the current Fortran organization

A practical C++ design should start from the existing conceptual decomposition rather than inventing a new one from scratch.

| Current Fortran area | Likely C++ responsibility |
| --- | --- |
| `parameters` | configuration/constants/runtime options structs |
| `elements`, `atom_types`, `group` | atom/element storage and metadata models |
| `box` | simulation box, boundary conditions, transforms |
| `neighbors` | neighbor-list builder and cell lists |
| `forces`, `potential`, `eam`, `morse`, `force_mod` | force pipeline and potential models |
| `integration`, `vel_verlet`, `dynamics`, `temp`, `langevin`, `berendsen`, `deform` | time integration and thermostats/barostats |
| `minimize`, `fire`, `cg`, `min_arrays` | minimization subsystem |
| `comms` | MPI topology, halo exchange, scatter/gather services |
| `read_data`, `input_parser`, `set`, `modify`, `displace` | input model, command parser, mutating commands |
| `dump`, `thermo`, `logger`, `time`, `errors` | output/logging/error/reporting services |

## 3. Migration Goals

The conversion should optimize for the following goals, in this order:

1. **Behavioral equivalence** for supported physics and commands.
2. **Deterministic reproducibility** where currently expected.
3. **Improved maintainability** through explicit ownership and smaller interfaces.
4. **Preserved MPI scalability characteristics** as much as practical.
5. **Incremental deliverability** so the code remains runnable during the transition.

Non-goals for the first migration wave:

- redesigning the scientific method itself,
- broad feature expansion,
- GPU offload,
- replacing MPI with a different parallel runtime,
- chasing micro-optimizations before correctness parity.

## 4. Recommended Target Architecture in C++

### 4.1 Language and tooling choices

- **Language**: C++20 preferred, C++17 acceptable if toolchain constraints require it.
- **Build system**: CMake.
- **Testing**:
  - GoogleTest or Catch2 for unit tests,
  - existing integration tests preserved initially via adapter scripts,
  - optional Python-based regression harness for golden-file comparisons.
- **Formatting/linting**:
  - `clang-format`,
  - `clang-tidy` for selected checks,
  - compiler warnings as errors in CI for new C++ targets.
- **Packaging/deps**:
  - avoid heavy dependencies unless they clearly replace bespoke logic,
  - keep MPI as an external toolchain dependency.

### 4.2 Data model

A good end-state decomposition would likely include:

- `SimulationConfig`
- `RuntimeContext`
- `MpiContext`
- `SimulationBox`
- `AtomData`
- `ElementData`
- `NeighborList`
- `Potential` base interface with concrete `EamPotential`, `MorsePotential`, etc.
- `Integrator` base interface with `VelocityVerletIntegrator`
- `Minimizer` base interface with `FireMinimizer` and `CgMinimizer`
- `Command` parser/executor layer
- `Logger`, `ThermoWriter`, `DumpWriter`

### 4.3 Data layout guidance

This code appears numerics-heavy and array-oriented. For performance and MPI interoperability, prefer:

- contiguous containers (`std::vector<T>` or PMR-backed containers where useful),
- structure-of-arrays layouts for particle and force data,
- explicit views/spans (`std::span`) for non-owning access,
- minimal virtual dispatch in hot loops,
- isolating polymorphism at setup/configuration boundaries rather than inside inner force loops.

### 4.4 MPI strategy

Use the **C MPI API from C++** rather than deprecated C++ MPI bindings. Wrap it in a thin RAII layer for:

- initialization/finalization,
- communicator ownership,
- rank/size queries,
- cartesian topology creation,
- collective wrappers,
- error reporting.

Do **not** over-abstract MPI early; a thin wrapper is safer than a bespoke messaging framework.

## 5. Major Risks and How to Control Them

### 5.1 Physics drift during translation

Risk: small indexing, initialization, or floating-point-order differences create silent simulation divergence.

Mitigation:

- preserve existing tests,
- create golden outputs for representative inputs,
- compare energies/forces/thermo output within tolerances,
- migrate pure helper logic first to establish validation patterns.

### 5.2 Over-big-bang rewrite failure

Risk: the team spends months rewriting without a runnable artifact.

Mitigation:

- keep a mixed-language period,
- migrate subsystem-by-subsystem,
- maintain a passing executable at each milestone,
- use compatibility shims where necessary.

### 5.3 Performance regressions

Risk: naïve C++ abstractions degrade memory locality or MPI exchange efficiency.

Mitigation:

- preserve data-oriented layouts,
- benchmark representative workloads after each major subsystem port,
- avoid object-per-atom designs,
- review hot loops before finalizing abstraction layers.

### 5.4 Hidden coupling from Fortran modules

Risk: implicit shared state in Fortran is not fully discovered until late.

Mitigation:

- inventory globals first,
- explicitly map module state into structs/classes,
- require new C++ APIs to take dependencies as arguments or constructor-injected collaborators.

## 6. Phased Migration Strategy

## Phase 0 — Discovery and Baseline Capture

### Objectives

- fully understand the current behavior surface,
- inventory architecture, dependencies, and global state,
- create a measurable migration baseline.

### Tasks

- [ ] **Inventory public behavior**
   - enumerate supported commands from docs and parser code,
   - map required command ordering and default behaviors,
   - identify all file formats read/written.

- [ ] **Inventory runtime modules**
   - produce a dependency graph for `src/*.f90`,
   - classify modules as pure, stateful, MPI-heavy, I/O-heavy, or orchestration-heavy,
   - identify modules with the most inbound dependencies.

- [ ] **Capture correctness baselines**
   - choose a small set of representative input decks,
   - record expected thermo output, dumps, energies, and minimization convergence behavior,
   - preserve exact or tolerance-based golden references.

- [ ] **Capture performance baselines**
   - measure startup time, step time, and memory for representative cases,
   - measure 1-rank and multi-rank behavior if the environment allows.

### Deliverables

- architecture inventory document,
- command/input compatibility matrix,
- baseline regression corpus,
- performance benchmark report.

## Phase 1 — Introduce a Mixed-Language Build and Test Harness

### Objectives

- make C++ code a first-class citizen without breaking the current Fortran build,
- allow incremental porting.

### Tasks

- [ ] Add a top-level **CMake** build that can:
   - build the existing Fortran executable,
   - build new C++ libraries/tests,
   - link mixed Fortran/C++ targets where needed.

- [ ] Establish compiler/toolchain policy:
   - MPI toolchain selection,
   - optimization/debug presets,
   - sanitizer-friendly debug build for C++.

- [ ] Add C++ unit test harness:
   - GoogleTest or Catch2,
   - one sample migrated helper module to prove build/test wiring.

- [ ] Add regression test runner:
   - shell or Python wrapper that runs the old executable and future C++ executable on the same inputs,
   - compare structured outputs.

### Deliverables

- `CMakeLists.txt` and supporting modules,
- CI path building Fortran + C++,
- first C++ unit tests,
- executable comparison harness.

## Phase 2 — Port Low-Risk, Deterministic Utility Modules First

### Objectives

- de-risk coding conventions, data types, and test methodology,
- create reusable utility patterns before touching core physics.

### Good early candidates

- `str`
- `input_validation`
- selected portions of `time`
- selected portions of `math`
- `dump_logic`

### Tasks

- [ ] Port pure helper functions to header/implementation pairs or small libraries.
- [ ] Match edge-case behavior with unit tests copied from current Fortran expectations.
- [ ] Define common numeric utilities:
   - tolerance comparison,
   - vector/matrix helpers,
   - error/status handling.
- [ ] Establish project conventions for:
   - namespaces,
   - file naming,
   - exception policy vs error returns,
   - `enum class` usage for command/options state.

### Deliverables

- `libcac_util` or equivalent,
- first passing C++ test suite,
- coding standards document.

## Phase 3 — Extract Explicit State Models Before Porting Heavy Logic

### Objectives

- replace implicit module globals with explicit C++ state containers,
- make later subsystem ports tractable.

### Tasks

- [ ] Create a **state inventory** for `parameters`, `elements`, `forces`, `integration`, `neighbors`, and `comms`.
- [ ] Group Fortran globals into coherent C++ structs such as:
   - `PhysicalConstants`,
   - `RuntimeFlags`,
   - `DomainState`,
   - `AtomState`,
   - `ForceState`,
   - `NeighborState`,
   - `IntegrationState`.
- [ ] Define ownership/lifetime rules.
- [ ] Decide which state is:
   - immutable configuration,
   - mutable per-step simulation state,
   - MPI/process-local metadata,
   - derived/cache data.
- [ ] Add serialization/debug dump helpers for state comparison against Fortran runs.

### Deliverables

- state model design doc,
- C++ types for core simulation state,
- unit tests for initialization/default behavior.

## Phase 4 — Port Input, Configuration, and Non-Hot-Path I/O

### Objectives

- create a C++ front end around the simulation before touching the deepest compute kernels.

### Candidate modules

- `logger`
- `errors`
- `time`
- `input_parser`
- `read_data`
- `set`
- `modify`
- `group`
- parts of `dump` and `thermo`

### Strategy

Split this phase into two sub-layers:

1. **Parsing/model layer**: tokenization, command recognition, validation.
2. **Execution layer**: invoking simulation services with parsed commands.

### Tasks

- [ ] Define a command AST or structured command model.
- [ ] Preserve command names and ordering requirements from the existing input system.
- [ ] Keep output formatting stable enough for regression diffs.
- [ ] Implement error messages with source line context.
- [ ] Add parser tests from real input snippets.

### Deliverables

- C++ parser library,
- C++ logging/error framework,
- compatibility-tested input execution path.

## Phase 5 — Port Geometry, Topology, and Neighbor Infrastructure

### Objectives

- establish the data movement and geometric primitives needed by force evaluation.

### Candidate modules

- `box`
- `elements`
- `atom_types`
- `integration` state pieces
- `neighbors`
- selected `comms` pieces

### Tasks

- [ ] Port simulation box and periodic-boundary operations.
- [ ] Port atom/element storage and resizing logic.
- [ ] Port neighbor-list data structures.
- [ ] Port cell decomposition and update logic.
- [ ] Port MPI scatter/gather and halo exchanges needed by neighbor rebuilds.
- [ ] Create deterministic geometry and neighbor regression tests.

### Deliverables

- C++ geometry/domain library,
- C++ particle storage library,
- passing neighbor-list regression suite.

## Phase 6 — Port Force Evaluation Pipeline

### Objectives

- migrate the scientific core while preserving result fidelity.

### Candidate modules

- `forces`
- `potential`
- `eam`
- `morse`
- `force_mod`
- `compute`

### Tasks

- [ ] Define a `Potential` interface that separates:
   - setup/loading,
   - per-step precomputation,
   - force/energy evaluation,
   - postprocessing.
- [ ] Port simple force paths first:
   - `compute`,
   - `morse`,
   - common force accumulation logic.
- [ ] Port EAM carefully with dedicated regression coverage:
   - spline generation,
   - map arrays,
   - electron density accumulation,
   - force update logic.
- [ ] Validate against Fortran outputs on small and medium systems.
- [ ] Benchmark hot loops and revisit data layout if needed.

### Deliverables

- C++ force/potential subsystem,
- tolerance-based energy/force parity reports,
- performance comparison report.

## Phase 7 — Port Integration, Dynamics, and Thermostat/Barostat Logic

### Objectives

- restore full time-evolution capabilities once force calculation is stable.

### Candidate modules

- `vel_verlet`
- `dynamics`
- `temp`
- `langevin`
- `berendsen`
- `deform`
- `quenched_dynamics`

### Tasks

- [ ] Define a clean time-step pipeline:
   - pre-step,
   - force update,
   - integration update,
   - post-step hooks,
   - thermo/dump checkpoints.
- [ ] Port velocity/position update kernels.
- [ ] Port thermostats/barostats with reproducibility controls for random components.
- [ ] Port deformation and ramp/displacement behaviors.
- [ ] Validate long-ish runs against thermodynamic trend baselines.

### Deliverables

- C++ runtime stepping engine,
- dynamics parity suite,
- documented random-seed policy.

## Phase 8 — Port Minimization Subsystem

### Objectives

- preserve static relaxation workflows after the dynamic runtime is stable.

### Candidate modules

- `minimize`
- `fire`
- `cg`
- `min_arrays`

### Tasks

- [ ] Port minimization state containers.
- [ ] Port FIRE before CG if it is operationally simpler.
- [ ] Port line search and convergence criteria with explicit tests.
- [ ] Compare convergence paths and final energies against Fortran references.

### Deliverables

- C++ minimization subsystem,
- minimizer regression suite,
- convergence tolerance policy.

## Phase 9 — Port Remaining MPI/Output Glue and Remove Fortran Runtime Dependency

### Objectives

- finish any bridging code,
- retire the Fortran executable once all required features are migrated.

### Tasks

- [ ] Replace remaining mixed-language interfaces.
- [ ] Ensure docs/examples run through the C++ executable.
- [ ] Remove now-obsolete Fortran build targets only after parity is confirmed.
- [ ] Preserve the old Fortran code on a maintenance branch/tag for scientific traceability.

### Deliverables

- standalone C++ executable,
- retired mixed-language bridge layer,
- migration completion report.

## 7. Suggested Subsystem Port Order

A practical order, balancing dependency risk and validation ease, is:

1. build/test harness,
2. pure utility modules,
3. state model extraction,
4. parser/logging/error layers,
5. box/elements/atom metadata,
6. neighbor infrastructure,
7. compute/common force accumulation,
8. morse potential,
9. EAM potential,
10. dynamics/integration,
11. thermo/dump output,
12. minimization,
13. cleanup and Fortran retirement.

This order is intentionally **not** the same as the current compile order. It is optimized for migration safety.

## 8. Detailed Work Breakdown by Existing Fortran Module

| Fortran file | Migration recommendation | Notes |
| --- | --- | --- |
| `parameters.f90` | redesign, not direct port | replace with explicit config/constants/state types |
| `str.f90` | direct/near-direct port | low risk |
| `input_validation.f90` | direct/near-direct port | low risk |
| `time.f90` | split port | timing/logging concerns separable from validation |
| `logger.f90` | redesign light | use streams/files + severity levels |
| `errors.f90` | redesign | prefer typed errors/exceptions or status objects |
| `math.f90` | selective port | separate pure numerics from MPI-dependent helpers |
| `box.f90` | redesign light | map well to `SimulationBox` |
| `atom_types.f90` | redesign light | map to metadata/lookup tables |
| `elements.f90` | redesign medium | central storage layer |
| `group.f90` | redesign medium | may become selection/group service |
| `comms.f90` | redesign heavy | likely one of the hardest ports |
| `neighbors.f90` | redesign heavy | major performance-sensitive subsystem |
| `forces.f90` | redesign medium | shared accumulation layer |
| `compute.f90` | port after state cleanup | useful for parity metrics |
| `potential.f90` | redesign medium | become force pipeline coordinator |
| `eam.f90` | redesign heavy | likely highest scientific risk |
| `morse.f90` | moderate port | good first concrete potential |
| `force_mod.f90` | moderate port | inject external/set forces cleanly |
| `integration.f90` | redesign medium | should become explicit runtime state |
| `vel_verlet.f90` | direct-ish port after state cleanup | hot loop, keep lean |
| `dynamics.f90` | redesign medium | orchestration layer |
| `temp.f90` | split port | RNG and rescaling need care |
| `langevin.f90` | split port | stochastic behavior must be controlled |
| `berendsen.f90` | moderate port | after box + compute state exist |
| `deform.f90` | moderate port | tied to box and neighbor rebuild logic |
| `read_data.f90` | redesign medium | input/model loading layer |
| `input_parser.f90` | redesign heavy | good AST/executor candidate |
| `dump.f90` | redesign heavy | separate formatting, gather, and file I/O |
| `thermo.f90` | moderate redesign | separate data collection from rendering |
| `modify.f90` | fold into command execution layer | low-medium risk |
| `set.f90` | fold into command execution layer | low-medium risk |
| `displace.f90` | fold into command execution layer | low-medium risk |
| `debug.f90` | postpone | lower priority |
| `min_arrays.f90` | port with minimizers | dependent on minimization design |
| `fire.f90` | port late | after force/dynamics state exists |
| `cg.f90` | port late | complex and coupled |
| `minimize.f90` | port late | orchestration layer |
| `main.f90` | rewrite last | becomes thin program entry point |

## 9. Validation Strategy

The migration should use three validation layers.

### 9.1 Unit-level parity

For deterministic helpers and pure math/state logic:

- exact comparisons for discrete outputs,
- tolerance comparisons for floating-point outputs,
- edge-case coverage for parser and boundary logic.

### 9.2 Subsystem regression

For larger pieces such as neighbor lists, EAM spline setup, dump formatting, and thermo output:

- golden files,
- structured diffs,
- tolerance-based numerical comparisons.

### 9.3 End-to-end simulation parity

For representative models:

- same input deck,
- same MPI rank counts where practical,
- compare final energies, selected forces, thermo traces, and convergence behavior.

Where exact bitwise equality is unrealistic, define acceptance tolerances per subsystem.

## 10. Build, CI, and Developer Workflow Plan

### Build-system milestones

1. **Milestone A**: CMake builds Fortran executable and C++ tests.
2. **Milestone B**: mixed-language bridge target exists.
3. **Milestone C**: C++ executable runs a minimal case.
4. **Milestone D**: Fortran build becomes optional.

### CI stages

1. formatting/lint,
2. C++ unit tests,
3. existing Fortran unit tests,
4. representative integration/regression tests,
5. optional performance smoke benchmarks.

### Developer conventions

- all new functionality lands in C++ once migration starts,
- Fortran changes should be limited to bug fixes or temporary bridge support,
- every subsystem port needs a validation note describing what was compared and how.

## 11. Documentation Plan

Update documentation in parallel with the migration.

### Needed docs

- architecture overview of the C++ design,
- build instructions for mixed-language and full-C++ modes,
- command compatibility matrix,
- migration status tracker by module,
- contributor guide for writing new C++ tests and benchmarks.

### Documentation rule

Whenever a Fortran module is superseded, document:

- what replaced it,
- what behaviors remain equivalent,
- what known gaps still exist.

## 12. Team and Resourcing Guidance

If multiple contributors are involved, split work by subsystem boundaries rather than by individual files where possible.

### Suggested ownership slices

- **Infrastructure owner**: CMake, CI, formatting, test harness.
- **Runtime owner**: parser, command execution, logging, errors.
- **Geometry/data owner**: box, elements, groups, neighbor state.
- **MPI owner**: comms and distributed-memory correctness.
- **Physics owner**: potentials, forces, compute, integration.
- **Minimization owner**: FIRE/CG/minimization flows.

A designated reviewer should own cross-cutting concerns:

- floating-point stability,
- MPI semantics,
- regression policy.

## 13. Rough Timeline Estimate

This is a substantial migration. A realistic timeline depends on team size and how much scientific verification is required.

### Solo or part-time effort

- likely **many months to over a year** for a safe full migration.

### Small focused team

- **3–6 months** for a credible mixed-language and partial-C++ system,
- **6–12 months** for a high-confidence full transition, depending heavily on EAM/comms complexity and validation depth.

The highest schedule risk is not syntax translation; it is validation of scientific equivalence in `comms`, `neighbors`, `eam`, and minimization behavior.

## 14. Immediate Next Steps for This Repository

The most useful next actions are:

- [ ] **Create a module/state inventory document**
   - enumerate every global variable in `parameters`, `elements`, `forces`, `neighbors`, `integration`, and `comms`.

- [ ] **Stand up a CMake skeleton**
   - build the existing Fortran sources first,
   - add one C++ utility library and one C++ test target.

- [ ] **Build a regression corpus**
   - pick 5–10 representative inputs covering:
     - read + potential setup,
     - neighbor generation,
     - thermo output,
     - dynamics,
     - minimization,
     - dump generation.

- [ ] **Port the lowest-risk modules first**
   - `str`, `input_validation`, selected `math`, `dump_logic`.

- [ ] **Design explicit core state types**
   - before porting `neighbors`, `forces`, or `potential`.

- [ ] **Choose a first medium-complexity pilot subsystem**
   - recommended: `box` + selected `elements` support logic,
   - not recommended as first pilot: `comms`, `neighbors`, or `eam`.

## 15. Definition of Success

The conversion should be considered successful when:

- the primary documented workflows run through a C++ executable,
- representative integration cases match Fortran behavior within agreed tolerances,
- the C++ build and tests are the default developer path,
- Fortran is no longer required for supported features,
- architecture is clearer than the current module-global design rather than merely being the same design rewritten in another syntax.

## 16. Final Recommendation

The repository is a **good candidate for conversion to C++**, but only if the work is approached as an **incremental architectural migration** instead of a direct transliteration. The biggest determinants of success will be:

- preserving behavior with regression tests,
- making global state explicit,
- porting MPI/data-layout-sensitive systems deliberately,
- delaying the hardest scientific kernels until the supporting architecture is ready.

If this plan is followed, the codebase can move toward a more maintainable, testable, and extensible modern C++ implementation without sacrificing the scientific behavior already embodied in the existing Fortran code.
