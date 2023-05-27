SatelliteToolboxGravityModels.jl
================================

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
