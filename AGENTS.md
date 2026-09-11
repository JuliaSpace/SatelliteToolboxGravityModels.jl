# Repository Guide

## Package Structure

- Single package (`SatelliteToolboxGravityModels`). Supports Julia 1.10 up to 1.13 (`[compat] julia = "1.10, 1.11, 1.12, 1.13"`).
- `src/SatelliteToolboxGravityModels.jl` is the module entrypoint and controls the include order: the `GravityModels` submodule (`src/GravityModels/GravityModels.jl`) is included and re-exported first, then `src/types.jl`, then the ICGEM files (`src/icgem/{api,compute,fetch,parse,show}.jl`).
- The `GravityModels` submodule defines the model API (`src/GravityModels/api.jl`), the `Workspace` with the Legendre buffers and precomputed coefficients plus the shared input processing (`src/GravityModels/workspace.jl`), the evaluation functions (`accelerations.jl`, `gravitational_field_derivative.jl`, `potential.jl`), and the time conversion helpers (`time.jl`). The submodule must be included before `src/types.jl` because `IcgemFile` subtypes `GravityModels.AbstractGravityModel`.
- `src/types.jl` defines `IcgemParseError`, the ICGEM coefficient types (`IcgemGfcCoefficient`, `IcgemPeriodicTerm`, `IcgemTimeVariableCoefficient`), and `IcgemFile{T, N}`, whose type parameter `N` is the `Val` with the Legendre normalization. Time-variable coefficients are stored sparsely in a vector plus an index (`time_variable_index`) defined only up to `max_time_variable_degree`.
- `ext/SatelliteToolboxGravityModelsZygoteExt.jl` is a package extension that loads only when **both** ForwardDiff and Zygote are loaded (see `[weakdeps]`/`[extensions]` in `Project.toml`); exercising it requires loading both trigger packages. Its pullbacks must not forward the user `workspace`, since ForwardDiff dual numbers cannot be stored in its buffers.
- Printing uses the public helpers of SatelliteToolboxBase.jl v2.1 (`print_tree`, `format_value`, `type_name`, `PrintedField`, `PrintedSection`), imported in the entrypoint.
- Tests do not mirror `src/` files. `test/runtests.jl` always runs `test/icgem.jl` and `test/gravity_models.jl`; on non-prerelease Julia it additionally `Pkg.add`s DifferentiationInterface, ForwardDiff, Zygote, JET, AllocCheck, and Aqua at runtime and runs `test/zygote_ext.jl` plus `test/allocations.jl` (allocation tests are skipped on macOS with Julia 1.12+ by `runtests.jl`, although AllocCheck works there with Julia 1.13).
- `test/icgem_test_files/` holds small hand-written ICGEM files, including `icgem2_time_variable.gfc` in the format 2.0 with validity intervals; `test/test_results/` holds reference grids computed by the ICGEM calculation service.
- Test-only dependencies are declared in `[extras]` + `[targets]` in `Project.toml` (Test, DelimitedFiles, LinearAlgebra, Pkg, SatelliteToolboxTransformations). SatelliteToolboxTransformations must be a release compatible with SatelliteToolboxBase v2 (v1.3.0 or later).
- Tests clear the package scratch space first, so ICGEM model files are re-downloaded during the run — network access is required, and the extra `Pkg.add` calls also hit the network.

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'` — the first run precompiles for minutes while printing little; use generous timeouts.
- Focused test file: test-only deps (e.g. SatelliteToolboxTransformations) are not available in a plain `--project=.` session, so use TestEnv.jl (must be installed in the default env): `julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/icgem.jl")'`. The test files use paths relative to `test/`, so run them from that directory.
- There is no test-name selector.
- CI (`.github/workflows/ci.yml`) builds via `julia-actions/julia-buildpkg` before testing and covers Julia 1.10 and latest stable on Ubuntu (x64), macOS (arm64), and Windows (x64), with coverage uploaded to Codecov; a separate workflow tests nightly on the same OSes.
- Build docs: `julia --project=docs docs/make.jl` (first run: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`). The build downloads EGM96 for the `@repl` blocks and fails on unresolved `@ref` links, so docstrings in the `GravityModels` submodule must not `@ref` names defined only in the main module. CI deploys docs from Julia stable via `julia-actions/julia-docdeploy`.
- Format: `julia -e 'using JuliaFormatter; format(".")'` — no `--project=.` (JuliaFormatter is not a package dependency; it must be available in the default env). Verify with `git diff --exit-code`.

## Code Style

- Formatting is configured in `.JuliaFormatter.toml` (Blue style with overrides: alignment options enabled, `whitespace_in_kwargs`, `whitespace_typedefs`, and several transformations like `pipe_to_function_call` disabled). The config file is the source of truth; CI does **not** run a format check, so apply the formatter manually after changes.
- Every function, macro, and structure — including private `_`-prefixed helpers — carries a docstring with the signature line ending in `-> <return value>`, followed by `# Arguments`, `# Keywords`, and `# Returns` sections when applicable; methods that differ only in the time argument type stack their signatures in one docstring. Exceptions: `Base` interface overloads (`show`, `showerror`, `zero`, `broadcastable`) and `ChainRulesCore.rrule` methods in the extension are undocumented by convention.
- Section separator comments (`# == Name ===...`) and file header blocks (`## Description ###...`) follow a fixed 92-column width, as do the dotted end-of-line comments; match the surrounding pattern. Docstring signature lines, URLs, and error message strings may exceed the width.
- Parser errors use `IcgemParseError` (with the line number when applicable); argument validation uses `ArgumentError`; invalid data lines are skipped with a `@warn` that starts with `[Line N]`.

## Behavioral Constraints

- New tests follow the `@testset "Name" verbose = true begin ... end` pattern used in `test/runtests.jl`; match it when adding coverage.
- Changes to differentiable code paths (potential, accelerations, field derivatives) must keep the ForwardDiff/Zygote extension tests (`test/zygote_ext.jl`) passing, and performance-sensitive paths are covered by JET/AllocCheck in `test/allocations.jl`, which check zero allocations when a `Workspace` is passed.
- The evaluation kernels evaluate the Legendre functions at the angle from the polar axis folded to the northern hemisphere and recover the southern hemisphere through parity; `test/gravity_models.jl` compares the results at the poles component-wise against `BigFloat` evaluations, so polar accuracy must be preserved.
- The printed representations are checked verbatim in `test/icgem.jl`; changing a label or the formatting requires updating those expectations.

## Not Configured

- No format-check CI job exists; formatting is applied manually (see Code Style).
- No linter or pre-commit hooks are configured.
- `Manifest.toml` is gitignored; do not commit it.
