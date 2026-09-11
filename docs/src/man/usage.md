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
GravityModels.load(::Type{T}, args...; kwargs...) where {T <: AbstractGravityModel} -> T
```

where the arguments and keywords depend on the gravity model type `T`. For ICGEM files, we
must use `T = IcgemFile` and the following signature:

```julia
GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::Type = Float64; kwargs...)
```

where it loads the ICGEM file in the path `filename` converting the coefficients to the type
`T`. The ICGEM format does not define the angular speed of the central body, which is
required to compute the gravity acceleration. Hence, it must be provided using the keyword
`angular_speed` [rad/s] for bodies other than Earth (**Default**: `EARTH_ANGULAR_SPEED`).
Both the ICGEM formats 1.0 and 2.0 are supported, including time-variable coefficients with
validity intervals.

We also provide a function to help downloading the ICGEM files:

```julia
fetch_icgem_file(url::AbstractString; kwargs...)
fetch_icgem_file(model::Symbol; kwargs...)
```

It fetches an ICGEM file from the `url` and returns its file path to be parsed with the
function [`GravityModels.load`](@ref). If the file already exists, it will not be
re-downloaded unless the keyword `force = true` is passed.

Notice that the function downloads the files to a [scratch
space](https://github.com/JuliaPackaging/Scratch.jl).

A symbol can be passed instead of the URL to fetch pre-configured gravity field models. The
supported values are:

- `:EGM96`: Earth Gravitational Model from 1996.
- `:EGM2008`: Earth Gravitational Model from 2008.
- `:JGM2`: Joint Gravity Model 2.
- `:JGM3`: Joint Gravity Model 3.

Finally, we can initialize, for example, the EGM96 model using:

```@repl usage
egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
```

## Workspace

All the functions that evaluate a model accept the keyword `workspace`, which receives an
object created by:

```julia
GravityModels.Workspace(model::AbstractGravityModel; kwargs...) -> Workspace
```

The workspace holds the buffers used to compute the associated Legendre functions and their
derivatives, together with precomputed recursion coefficients. Hence, when the functions
are called many times for the same model, e.g. in a numerical orbit propagator, the
workspace avoids allocations and largely improves the performance. The following keywords
are available:

- `max_degree::Int`: Maximum degree supported by the workspace. If it is higher than the
    maximum degree of the model, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree of the model.
    (**Default**: -1)
- `max_order::Int`: Maximum order supported by the workspace. If it is higher than
    `max_degree`, it will be clamped. If it is lower than 0, it will be set to the same
    value as `max_degree`.
    (**Default**: -1)
- `T::Type{<:AbstractFloat}`: Element type of the workspace, which must be the type obtained
    by promoting the type of the model coefficients, the element type of the position, and
    the type of the time used in the evaluations.
    (**Default**: type of the model coefficients)

The evaluation functions throw an `ArgumentError` if the workspace element type or the
supported degree and order do not match the computation.

!!! warning

    The workspace holds mutable buffers. Hence, it must not be shared among threads that
    evaluate the model concurrently. Create one workspace per thread instead.

```@repl usage
workspace = GravityModels.Workspace(egm96)
```

## Common Keywords

The functions described in the following sections accept the keywords:

- `max_degree::Int`: Maximum degree used in the spherical harmonics. If it is higher than
    the available number of coefficients in the model, it will be clamped. If it is lower
    than 0, it will be set to the maximum degree available.
    (**Default**: -1)
- `max_order::Int`: Maximum order used in the spherical harmonics. If it is higher than
    `max_degree`, it will be clamped. If it is lower than 0, it will be set to the same
    value as `max_degree`.
    (**Default**: -1)
- `workspace::Union{Nothing, Workspace}`: Workspace created for the model as described in
    the previous section. If it is `nothing`, the buffers are allocated at every call.
    (**Default**: `nothing`)

The time can be passed as a `DateTime` object or as the number of elapsed seconds [s] from
the J2000.0 epoch (2000-01-01T12:00:00). If it is omitted, the J2000.0 epoch is used.

## Gravitational Potential

The following function:

```julia
GravityModels.gravitational_potential(model::AbstractGravityModel, r::AbstractVector, time = 0; kwargs...) -> RT
```

computes the gravitational potential [m²/s²] using the `model` in the position `r` [m],
represented in the body-fixed frame (ITRF for Earth), at instant `time`. The gravitational
potential is the potential caused by the central body mass only, i.e., without considering
the centrifugal potential.

```@repl usage
GravityModels.gravitational_potential(egm96, [6378.137e3, 0, 0]; workspace)
```

## Gravitational Field Derivative

The following function:

```julia
GravityModels.gravitational_field_derivative(model::AbstractGravityModel, r::AbstractVector, time = 0; kwargs...) -> RT, RT, RT
```

computes the gravitational field derivative with respect to the spherical coordinates:

```math
\frac{\partial U}{\partial r},~ \frac{\partial U}{\partial \phi},~ \frac{\partial U}{\partial \lambda},~
```

using the `model` in the position `r` [m], represented in the body-fixed frame (ITRF for
Earth), at instant `time`. The derivatives have units [m/s²], [m²/s²], and [m²/s²],
respectively.

!!! info

    In this case, $$\phi$$ is the geocentric latitude and $$\lambda$$ is the longitude.

```@repl usage
GravityModels.gravitational_field_derivative(egm96, [6378.137e3, 0, 0]; workspace)
```

## Gravitational Acceleration

The gravitational acceleration is the acceleration caused by the central body mass only,
i.e., without considering the centrifugal potential. We can compute it using the function:

```julia
GravityModels.gravitational_acceleration(model::AbstractGravityModel, r::AbstractVector, time = 0; kwargs...) -> SVector{3, RT}
```

where it returns the gravitational acceleration [m/s²] represented in the body-fixed frame
(ITRF for Earth) using the `model` in the position `r` [m], also represented in the
body-fixed frame, at instant `time`.

```@repl usage
GravityModels.gravitational_acceleration(egm96, [6378.137e3, 0, 0]; workspace)
```

The algorithm is accurate at the poles, including positions exactly on the polar axis:

```@repl usage
GravityModels.gravitational_acceleration(egm96, [0, 0, 6356.7523e3]; workspace)
```

## Gravity Acceleration

The gravity acceleration is the compound acceleration caused by the central body mass and
the centrifugal force due to the body's rotation. We can compute it using the function:

```julia
GravityModels.gravity_acceleration(model::AbstractGravityModel, r::AbstractVector, time = 0; kwargs...) -> SVector{3, RT}
```

where it computes the gravity acceleration [m/s²] represented in the body-fixed frame (ITRF
for Earth) using the `model` in the position `r` [m], also represented in the body-fixed
frame, at instant `time`. Besides the common keywords, this function accepts:

- `ω::Number`: Angular speed of the body [rad/s], which defaults to the value stored in the
    model (see [`GravityModels.angular_speed`](@ref)).
    (**Default**: `GravityModels.angular_speed(model)`)

Thus, we can compute the gravity acceleration in the Equator using the EGM96 model by:

```@repl usage
GravityModels.gravity_acceleration(egm96, [6378.137e3, 0, 0]; workspace)
```

Whereas we can obtain the gravity acceleration at the North pole by:

```@repl usage
GravityModels.gravity_acceleration(egm96, [0, 0, 6356.7523e3]; workspace)
```

## Automatic Differentiation

The evaluation functions can be differentiated with respect to the position and the time
using **ForwardDiff.jl**. If both **ForwardDiff.jl** and **Zygote.jl** are loaded, a package
extension provides the reverse rules required by **Zygote.jl**, whose pullbacks compute the
Jacobians with **ForwardDiff.jl**. Notice that the workspace cannot be used in this case
because its buffers cannot store the dual numbers.
