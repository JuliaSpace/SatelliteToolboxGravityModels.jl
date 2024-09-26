# GravityModels API

This document describes the API required for a gravity model. The user can add new models by
overloading the functions listed here.

## Structure

All models require a structure with supertype `AbstractGravityModel{T<:Number}`, where `T`
is the type of the coefficients in the model. For the function that passes the seconds since
J2000.0, the return type `RT` is the promoted type of the set {`T`, `W`}, where `W` is the
type of the time variable.

## API Functions

```julia
function coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::Number) where {T<:Number}> -> T, T
```

This function must return the coefficients `Clm` and `Slm` of the gravity `model` for the
specified `degree`, `order`, and `time`. The latter is described as the number of ellapsed
seconds since J2000.0 epoch (2000-01-01T12:00:00.000). Hence:

```julia
coefficients(model, 10, 8, 7.404048e8)
```

must return a `Tuple{T, T}` with the `Clm` and `Slm`, respectively, for the degree 10, order
8, and computed at 2023-06-19T00:00:00.000.

> **Note**
> If the model has constant coefficients, the function must still accept the positional
> argument `time`, but it will be neglected. The package already defines the function
> without the `time` for the sake of usage simplification.

> **Note**
> The time can also be specified as a `DateTime` object. The package already defines the
> function with this signature.

---

```julia
function coefficient_norm(model::AbstractGravityModel) -> Symbol
```

This function must return the normalization we must use in the spherical harmonics when
computing the Legendre associated functions. The accepted values are:

- `:full`: Use full normalization.
- `:schmidt`: Use Schmidt quasi-normalization.
- `:unnormalized`: Do not perform normalization.

---

```julia
function gravity_constant(model::AbstractGravityModel{T}) where T<:Number -> T
```

This function must return the gravity constant [m³ / s²] for the gravity model.

---

```julia
function load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T
```

This function must return gravity model structure, which is loaded using the arguments
`args...` and keywords `kwargs...`.

---

```julia
function maximum_degree(model::AbstractGravityModel) -> Int
```

This function must return the maximum degree of the gravity `model`.

---

```julia
function radius(model::AbstractGravityModel{T}) where T<:Number -> T
```

This function must return the radius [m] for the gravity model.
