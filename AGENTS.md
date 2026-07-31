# Repository Guide

## Package Structure

- Single package (`SatelliteToolboxGravityModels`). Supports Julia 1.10 up to 1.12 (`[compat] julia = "1.10, 1.11, 1.12"`).
- `src/SatelliteToolboxGravityModels.jl` is the module entrypoint and controls the include order: the `GravityModels` submodule (`src/GravityModels/GravityModels.jl`) is included and re-exported first, then `src/types.jl`, then the ICGEM files (`src/icgem/{api,compute,fetch,parse,show}.jl`).
- `ext/SatelliteToolboxGravityModelsZygoteExt.jl` is a package extension that loads only when **both** ForwardDiff and Zygote are loaded (see `[weakdeps]`/`[extensions]` in `Project.toml`); exercising it requires loading both trigger packages.
- Tests do not mirror `src/` files. `test/runtests.jl` always runs `test/icgem.jl` and `test/gravity_models.jl`; on non-prerelease Julia it additionally `Pkg.add`s DifferentiationInterface, ForwardDiff, Zygote, JET, AllocCheck, and Aqua at runtime and runs `test/zygote_ext.jl` plus `test/allocations.jl` (allocation tests are skipped on macOS with Julia 1.12+).
- Test-only dependencies are declared in `[extras]` + `[targets]` in `Project.toml` (Test, DelimitedFiles, LinearAlgebra, Pkg, SatelliteToolboxTransformations).
- Tests clear the package scratch space first, so ICGEM model files are re-downloaded during the run — network access is required, and the extra `Pkg.add` calls also hit the network.

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'` — the first run precompiles for minutes while printing little; use generous timeouts.
- Focused test file: test-only deps (e.g. SatelliteToolboxTransformations) are not available in a plain `--project=.` session, so use TestEnv.jl (must be installed in the default env): `julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/icgem.jl")'`
- There is no test-name selector.
- CI (`.github/workflows/ci.yml`) builds via `julia-actions/julia-buildpkg` before testing and covers Julia 1.10 and latest stable on Ubuntu (x64), macOS (arm64), and Windows (x64), with coverage uploaded to Codecov; a separate workflow tests nightly on the same OSes.
- Build docs: `julia --project=docs docs/make.jl` (first run: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'`). CI deploys docs from Julia stable via `julia-actions/julia-docdeploy`.
- Format: `julia -e 'using JuliaFormatter; format(".")'` — no `--project=.` (JuliaFormatter is not a package dependency; it must be available in the default env). Verify with `git diff --exit-code`.

## Code Style

- Formatting is configured in `.JuliaFormatter.toml` (Blue style with overrides: alignment options enabled, `whitespace_in_kwargs`, `whitespace_typedefs`, and several transformations like `pipe_to_function_call` disabled). The config file is the source of truth; CI does **not** run a format check, so apply the formatter manually after changes.
- Every function, macro, and structure — including private `_`-prefixed helpers — carries a docstring with the signature line ending in `-> <return value>`, followed by `# Arguments`, `# Keywords`, and `# Returns` sections when applicable; match this template. Exceptions: `Base` interface overloads (`show`, `zero`) and `ChainRulesCore.rrule` methods in the extension are undocumented by convention.
- Section separator comments (`# == Name ===...`) and file header blocks (`## Description ###...`) follow a fixed 92-column width; match the surrounding pattern.

## Behavioral Constraints

- New tests follow the `@testset "Name" verbose = true begin ... end` pattern used in `test/runtests.jl`; match it when adding coverage.
- Changes to differentiable code paths (potential, accelerations, field derivatives) must keep the ForwardDiff/Zygote extension tests (`test/zygote_ext.jl`) passing, and performance-sensitive paths are covered by JET/AllocCheck in `test/allocations.jl`.

## Not Configured

- No format-check CI job exists; formatting is applied manually (see Code Style).
- No linter or pre-commit hooks are configured.
- `Manifest.toml` is gitignored; do not commit it.
