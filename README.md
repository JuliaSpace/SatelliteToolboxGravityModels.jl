SatelliteToolboxGravityModels.jl
================================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolboxGravityModels.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxGravityModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package implements the support of gravity models for the **SatelliteToolbox.jl**
ecosystem. We can use it to obtain, for example, the gravitational acceleration to
implement highly accurate numerical orbit propagators.

Currently, we have the following functionalities:

- Compute the gravity field derivative in spherical coordinates;
- Compute the gravitational acceleration; and
- Compute the gravity acceleration, taking into account the Earth's rotation rate.

The package contains an API to allow the user to access any gravity model while computing
the previous functions. We also provide in-built support for
[ICGEM](http://icgem.gfz-potsdam.de/home) files.

## Installation

```julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxGravityModels")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/stable
