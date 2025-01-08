# GravityModels API

This document describes the API required for a gravity model. The user can add new models by
overloading the functions listed here.

## Structure

All models require a structure with supertype `AbstractGravityModel{T<:Number, NT}`, where `T`
is the type of the coefficients in the model and `NT` is the norm type.

## API Functions

```julia
function coefficients(model::AbstractGravityModel{T, NT}, degree::Int, order::Int, time::DataTime = DateTime("2000-01-01")) where {T<:Number, NT} -> T, T
```

This function must return the coefficients `Clm` and `Slm` of the gravity `model` for the
specified `degree`, `order`, and `time`. If the latter argument is omitted, the J2000.0
epoch is used. Hence:

```julia
coefficients(model, 10, 8, DateTime("2023-06-19"))
```

must return a `Tuple{T, T}` with the `Clm` and `Slm`, respectively, for the degree 10, order
8, and computed at day 2023-06-19.

> **Note**
> If the model has constant coefficients, the function must still accept the positional
> argument `time`, but it will be neglected. The package already defines the function
> without the `time` for the sake of usage simplification.

---

```julia
function coefficient_norm(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> Symbol
```

This function must return the normalization we must use in the spherical harmonics when
computing the Legendre associated functions. The accepted values are:

- `:full`: Use full normalization.
- `:schmidt`: Use Schmidt quasi-normalization.
- `:unnormalized`: Do not perform normalization.

---

```julia
function gravity_constant(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> T
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
function radius(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> T
```

This function must return the radius [m] for the gravity model.
