## Description #############################################################################
#
# Define the API functions for the gravity models.
#
############################################################################################

"""
    coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::Number = -43200) where T<:Number -> T, T

Return the `Clm` and `Slm` coefficients of the gravity `model` for the specified `degree`,
`order`, and `time`. If the latter argument is omitted, the J2000.0 epoch is used.
"""
function coefficients end

function coefficients(model::AbstractGravityModel, degree::Int, order::Int, time::Number = -43200)
    return coefficients(model, degree, order, time)
end

"""
    coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::DateTime = DateTime("2000-01-01")) where T<:Number -> T, T

Return the `Clm` and `Slm` coefficients of the gravity `model` for the specified `degree`,
`order`, and `time`. If the latter argument is omitted, the J2000.0 epoch is used.
"""
function coefficients(model::AbstractGravityModel, degree::Int, order::Int, time::DateTime = DateTime("2000-01-01"))
    return coefficients(model, degree, order, time)
end

"""
    coefficient_norm(model::AbstractGravityModel) where T<:Number -> Symbol

Return the normalization we must use in the spherical harmonics when computing the Legendre
associated functions. The accepted values are:

- `:full`: Use full normalization.
- `:schmidt`: Use Schmidt quasi-normalization.
- `:unnormalized`: Do not perform normalization.
"""
function coefficient_norm end

"""
    gravity_constant(model::AbstractGravityModel{T}) where T<:Number -> T

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
    radius(model::AbstractGravityModel{T}) where T<:Number -> T

Return the radius [m] for the gravity model.
"""
function radius end
