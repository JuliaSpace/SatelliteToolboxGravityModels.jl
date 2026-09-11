# GravityModels API

This document describes the API required for a gravity model. The user can add new models by
overloading the functions listed here.

## Structure

All models require a structure with supertype `AbstractGravityModel{T <: Number}`, where `T`
is the type of the coefficients in the model.

## API Functions

```julia
function coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::Number) where {T <: Number} -> RT, RT
```

This function must return the coefficients `Clm` and `Slm` of the gravity `model` for the
specified `degree`, `order`, and `time`, expressed as the number of elapsed seconds from the
J2000.0 epoch (2000-01-01T12:00:00). Hence:

```julia
coefficients(model, 10, 8, 0.0)
```

must return a tuple with the `Clm` and `Slm`, respectively, for the degree 10, order 8, and
computed at the J2000.0 epoch. The return type `RT` is `T` or its promotion with the type of
`time` if the model has time-variable coefficients.

> **Note**
> If the model has constant coefficients, the function must still accept the positional
> argument `time`, but it will be neglected. The package already defines the methods that
> receive a `DateTime` object and that omit the `time` for the sake of usage simplification.

---

```julia
function angular_speed(model::AbstractGravityModel{T}) where {T <: Number} -> T
```

This function must return the angular speed [rad/s] of the central body, which is used to
compute the centrifugal acceleration in `gravity_acceleration`.

---

```julia
function coefficient_norm(model::AbstractGravityModel) -> Val
```

This function must return the normalization we must use in the spherical harmonics when
computing the Legendre associated functions, wrapped in a `Val`. The accepted values are:

- `Val(:full)`: Use full normalization.
- `Val(:schmidt)`: Use Schmidt quasi-normalization.
- `Val(:unnormalized)`: Do not perform normalization.

The return type must be inferable from the type of `model`, e.g. by storing the `Val` in a
type parameter, so that the evaluation functions are type stable.

---

```julia
function gravity_constant(model::AbstractGravityModel{T}) where {T <: Number} -> T
```

This function must return the gravity constant [m³/s²] for the gravity model.

---

```julia
function load(::Type{T}, args...; kwargs...) where {T <: AbstractGravityModel} -> T
```

This function must return the gravity model structure, which is loaded using the arguments
`args...` and keywords `kwargs...`.

---

```julia
function maximum_degree(model::AbstractGravityModel) -> Int
```

This function must return the maximum degree of the gravity `model`.

---

```julia
function radius(model::AbstractGravityModel{T}) where {T <: Number} -> T
```

This function must return the reference radius [m] for the gravity model.
