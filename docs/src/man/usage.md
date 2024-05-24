# Usage

```@meta
CurrentModule = SatelliteToolboxGravityModels
```

```@repl usage
using SatelliteToolboxGravityModels
```

## Initialization

We can initialize a gravity model using the function:

```julia
load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T
```

where the arguments and keywords depend on the gravity model type `T`. For ICGEM files, we
must use `T = IcgemFile` and the following signature:

```julia
GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::DataType = Float64)
```

where it loads the ICGEM file in the path `filename` converting the coefficients to the type
`T`.

We also provide a function to help downloading the ICGEM files:

```julia
fetch_icgem_file(url::AbstractString; kwargs...)
fetch_icgem_file(model::Symbol; kwargs...)
```

It fetches a ICGEM file from the `url` and return its file path to be parsed with the
function [`GravityModels.load`](@ref). If the file already exists, it will not be
re-downloaded unless the keyword `force = true` is passed.

Notice that the functions downloads the files to a [scratch
space](https://github.com/JuliaPackaging/Scratch.jl).

A symbol can be passed instead the URL to fetch pre-configured gravity field models. The
supported values are:

- `:EGM96`: Earth Gravitational Model from 1996.
- `:EGM2008`: Earth Gravitational Model from 2008.
- `:JGM2`: Joint Gravity Model 2.
- `:JGM3`: Joint Gravity Model 3.

Finally, we can initialize, for example, the EGM96 model using:

```@repl usage
egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
```

## Gravitational Field Derivative

The following function:

```julia
gravitational_field_derivative(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime("2000-01-01"); kwargs...) where T<:Number -> NTuple{3, T}
```

computes the gravitational field derivative [SI] with respect to the spherical coordinates:

```math
\frac{\partial U}{\partial r},~ \frac{\partial U}{\partial \phi},~ \frac{\partial U}{\partial \lambda},~
```

using the `model` in the position `r` [m], represented in ITRF, at instant `time`. If the
latter argument is omitted, the J2000.0 epoch is used.

!!! info

    In this case, $$\phi$$ is the geocentric latitude and $$\lambda$$ is the longitude.

The following keywords are available:

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available.
    (**Default** = -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default** = -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default** = `nothing`)

```@repl usage
GravityModels.gravitational_field_derivative(egm96, [6378.137e3, 0, 0])
```

## Gravitational Acceleration

The gravitational acceleration is the acceleration caused by the central body mass only,
i.e., without considering the centrifugal potential. We can compute it using the function:

```julia
gravitational_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime("2000-01-01"); kwargs...) where T<:Number -> NTuple{3, T}
```

where it returns the gravitational field acceleration [m / s²] represented in ITRF using the
`model` in the position `r` [m], also represented in ITRF, at instant `time`. If the latter
argument is omitted, the J2000.0 epoch is used.

The following keywords are available:

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available.
    (**Default** = -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default** = -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default** = `nothing`)

```@repl usage
GravityModels.gravitational_acceleration(egm96, [6378.137e3, 0, 0])
```

## Gravity Acceleration

The gravity acceleration is the compound acceleration caused by the central body mass and
the centrifugal force due to the planet's rotation. We can compute it using the function:

```julia
gravity_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime("2000-01-01"); kwargs...) where T<:Number -> NTuple{3, T}
```

where it computes the gravity acceleration [m / s²] represented in ITRF using the `model` in
the position `r` [m], also represented in ITRF, at instant `time`. If the latter argument is
omitted, the J2000.0 epoch is used.

The following keywords are available:

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available.
    (**Default** = -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default** = -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default** = `nothing`)

Thus, we can compute the gravity acceleration in Equator using the EGM96 model by:

```@repl usage
egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96));
```

Whereas we can obtain the gravity acceleration at the poles by:

```@repl usage
GravityModels.gravitational_acceleration(egm96, [0, 0, 6356.7523e3])
```
