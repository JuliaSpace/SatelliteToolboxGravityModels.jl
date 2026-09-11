## Description #############################################################################
#
# Define the API functions for the gravity models.
#
############################################################################################

"""
    coefficients(model::AbstractGravityModel, degree::Int, order::Int[, time]) -> T, T

Return the `Clm` and `Slm` coefficients [-] of the gravity `model` for the specified
`degree`, `order`, and `time`. If the latter argument is omitted, the J2000.0 epoch
(2000-01-01T12:00:00) is used.

# Arguments

- `model::AbstractGravityModel{T}`: Gravity model.
- `degree::Int`: Degree of the coefficients.
- `order::Int`: Order of the coefficients.
- `time::Union{Number, DateTime}`: Time at which the coefficients are computed, expressed
    as a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.
    (**Default**: J2000.0 epoch)

# Returns

- `T`: Coefficient `Clm` [-] for the specified `degree`, `order`, and `time`.
- `T`: Coefficient `Slm` [-] for the specified `degree`, `order`, and `time`.
"""
function coefficients end

function coefficients(model::AbstractGravityModel, degree::Int, order::Int)
    return coefficients(model, degree, order, 0)
end

function coefficients(model::AbstractGravityModel, degree::Int, order::Int, time::DateTime)
    return coefficients(model, degree, order, _to_j2000_seconds(time))
end

"""
    angular_speed(model::AbstractGravityModel{T}) -> T

Return the angular speed [rad/s] of the central body of the gravity `model`, which is used
to compute the centrifugal acceleration in [`gravity_acceleration`](@ref).
"""
function angular_speed end

"""
    coefficient_norm(model::AbstractGravityModel) -> Val

Return the normalization we must use in the spherical harmonics when computing the
Legendre associated functions for the gravity `model`, wrapped in a `Val` so that the
Legendre functions can be dispatched on it. The accepted values are:

- `Val(:full)`: Use full normalization.
- `Val(:schmidt)`: Use Schmidt quasi-normalization.
- `Val(:unnormalized)`: Do not perform normalization.
"""
function coefficient_norm end

"""
    gravity_constant(model::AbstractGravityModel{T}) -> T

Return the gravity constant [m³/s²] of the gravity `model`.
"""
function gravity_constant end

"""
    load(::Type{T}, args...; kwargs...) -> T

Load a gravity model of type `T` using the arguments `args...` and keywords `kwargs...`.
"""
function load end

"""
    maximum_degree(model::AbstractGravityModel) -> Int

Return the maximum degree of the gravity `model`.
"""
function maximum_degree end

"""
    radius(model::AbstractGravityModel{T}) -> T

Return the reference radius [m] of the gravity `model`.
"""
function radius end
