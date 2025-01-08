## Description #############################################################################
#
# Define the API functions for the gravity models.
#
############################################################################################

"""
    coefficients(model::AbstractGravityModel{T, NT}, degree::Int, order::Int[, time::Union{Number, DateTime}]) where {T<:Number, NT} -> T, T

Return the `Clm` and `Slm` coefficients of the gravity `model` for the specified `degree`,
`order`, and `time`. If the latter argument is omitted, the J2000.0 epoch is used.

`time` can be expressed using a `DateTime` object or the number of ellapsed seconds from
J2000.0 epoch.
"""
function coefficients end

function coefficients(model::AbstractGravityModel, degree::Int, order::Int)
    return coefficients(model, degree, order, 0)
end

function coefficients(model::AbstractGravityModel, degree::Int, order::Int, time::DateTime)
    t = Dates.value(time - _DT_J2000) / 1000
    return coefficients(model, degree, order, t)
end

"""
    coefficient_norm(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> Symbol

Return the normalization we must use in the spherical harmonics when computing the Legendre
associated functions. The accepted values are:

- `:full`: Use full normalization.
- `:schmidt`: Use Schmidt quasi-normalization.
- `:unnormalized`: Do not perform normalization.
"""
function coefficient_norm end

"""
    gravity_constant(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> T

Return the gravity constant [m³ / s²] for the gravity model.
"""
function gravity_constant end

"""
    load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T

Load a gravity model of type `T` using the arguments `args...` and keywords `kwargs...`.
"""
function load end

"""
    maximum_degree(model::AbstractGravityModel) -> Int

Return the maximum degree of the gravity `model`.
"""
function maximum_degree end

"""
    radius(model::AbstractGravityModel{T, NT}) where {T<:Number, NT} -> T

Return the radius [m] for the gravity model.
"""
function radius end
