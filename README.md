<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

SatelliteToolboxGravityModels.jl
================================

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteToolboxGravityModels.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteToolboxGravityModels.jl?token=47G4OLV6PD&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxGravityModels.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-16A34A?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-D97706?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&logo=julia&logoColor=white&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteToolboxGravityModels.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/blob/main/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.10396094-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.10396094)

This package implements the support of gravity models for the **SatelliteToolbox.jl**
ecosystem. We can use it to obtain, for example, the gravitational acceleration to
implement highly accurate numerical orbit propagators.

Currently, we have the following functionalities:

- Compute the gravity field derivative in spherical coordinates;
- Compute the gravitational acceleration; and
- Compute the gravity acceleration, taking into account the body's rotation rate.

The package contains an API to allow the user to access any gravity model while computing
the previous functions. We also provide in-built support for
[ICGEM](http://icgem.gfz-potsdam.de/home) files.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxGravityModels")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/stable
