## Description #############################################################################
#
# Function to compute accelerations.
#
## References ##############################################################################
#
# [1] Barthelmes, F (2013). Definition of Functions of the Geopotential and Their
#     Calculation from Spherical Harmonic Models. Scientific Technical Report STR09/02.
#     GeoForschungsZentrum (GFZ).
#
############################################################################################

"""
    gravitational_acceleration(model::AbstractGravityModel, r::AbstractVector[, time]; kwargs...) -> SVector{3, RT}

Compute the gravitational acceleration [m/s²] represented in the body-fixed frame (ITRF
for Earth) using the `model` in the position `r` [m], also represented in the body-fixed
frame, at instant `time`. If the latter argument is omitted, the J2000.0 epoch
(2000-01-01T12:00:00) is used.

The return element type `RT` is obtained by promoting the type of the `model`
coefficients, the element type of `r`, and the type of `time`.

!!! note

    Gravitational acceleration is the acceleration caused by the central body mass only,
    i.e., without considering the centrifugal potential.

!!! note

    The matrices `P` and `dP` are lower triangular. Hence, the algorithm performance for
    large models can be improved if they are created using the `LowerTriangularStorage`
    (defined in SatelliteToolboxBase.jl) with a row-major ordering. If those matrices are
    not provided by the user, they will be created using that type of storage.

See also: [`gravity_acceleration`](@ref)

# Arguments

- `model::AbstractGravityModel{T}`: Gravity model.
- `r::AbstractVector`: Position [m] in the body-fixed frame (ITRF for Earth) at which the
    acceleration is computed.
- `time::Union{Number, DateTime}`: Time at which the acceleration is computed, expressed
    as a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.
    (**Default**: J2000.0 epoch)

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available.
    (**Default**: -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default**: -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default**: `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default**: `nothing`)

# Returns

- `SVector{3, RT}`: Gravitational acceleration [m/s²] represented in the body-fixed frame
    (ITRF for Earth).

# References

- **[1]** Barthelmes, F (2013). *Definition of Functions of the Geopotential and Their
    Calculation from Spherical Harmonic Models*. Scientific Technical Report STR09/02.
    GeoForschungsZentrum (GFZ), p. 22.
"""
function gravitational_acceleration(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::Number = 0;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number}
    RT = promote_type(T, V, typeof(time))

    # == Partial Derivatives of the Gravitational Field ====================================

    n_max, m_max, n_max_P, m_max_P, n_max_dP, m_max_dP, P, dP = _prepare_field_derivative_inputs(
        model, RT, max_degree, max_order, P, dP
    )

    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ, ∂U_∂λ_over_cosϕ_pole = _gravitational_field_derivative_kernel(
        model, r, time, n_max, m_max, n_max_P, m_max_P, n_max_dP, m_max_dP, P, dP
    )

    # == Acceleration Represented in the UEN Frame =========================================

    r_gc, ρ_gc, _, sin_λ, cos_λ, _ = _spherical_coordinates(r, RT)

    sin_ϕ = r[3] / r_gc
    cos_ϕ = ρ_gc / r_gc

    # Compute the partial derivatives in spherical coordinate systems [1, p. 22] (eq. 120):
    #
    #     ∂U        1      ∂U     1   ∂U
    #    ---- , ---------.---- , ---.----
    #     ∂r     r.cos ϕ   ∂λ     r   ∂ϕ
    #
    # Notice that r ⋅ cos(ϕ_gc) = ρ_gc. On the polar axis, both `∂U / ∂λ` and `ρ_gc` are
    # 0 and we must use the limit computed by the kernel.
    a_u = ∂U_∂r
    a_e = (ρ_gc > 0) ? ∂U_∂λ / ρ_gc : ∂U_∂λ_over_cosϕ_pole / r_gc
    a_n = ∂U_∂ϕ / r_gc

    # == Acceleration Represented in the ITRF ==============================================

    # The unit vectors of the local UEN (Up-East-North) reference frame represented in the
    # body-fixed frame (ITRF for Earth) are:
    #
    #   up    = [ cos(ϕ) cos(λ),  cos(ϕ) sin(λ), sin(ϕ) ],
    #   east  = [        -sin(λ),         cos(λ),      0 ],
    #   north = [-sin(ϕ) cos(λ), -sin(ϕ) sin(λ), cos(ϕ) ].
    a_itrf = @SVector [
        a_u * cos_ϕ * cos_λ - a_e * sin_λ - a_n * sin_ϕ * cos_λ,
        a_u * cos_ϕ * sin_λ + a_e * cos_λ - a_n * sin_ϕ * sin_λ,
        a_u * sin_ϕ + a_n * cos_ϕ,
    ]

    return a_itrf
end

function gravitational_acceleration(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number}
    return gravitational_acceleration(
        model,
        r,
        _to_j2000_seconds(time);
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
    )
end

"""
    gravity_acceleration(model::AbstractGravityModel, r::AbstractVector[, time]; kwargs...) -> SVector{3, RT}

Compute the gravity acceleration [m/s²] represented in the body-fixed frame (ITRF for
Earth) using the `model` in the position `r` [m], also represented in the body-fixed
frame, at instant `time`. If the latter argument is omitted, the J2000.0 epoch
(2000-01-01T12:00:00) is used.

The return element type `RT` is obtained by promoting the type of the `model`
coefficients, the element type of `r`, and the type of `time`.

!!! note

    Gravity acceleration is the compound acceleration caused by the central body mass and
    the centrifugal force due to the planet's rotation.

    For non-Earth bodies, the body's rotation rate can be provided using the `ω` keyword.

!!! note

    The matrices `P` and `dP` are lower triangular. Hence, the algorithm performance for
    large models can be improved if they are created using the `LowerTriangularStorage`
    (defined in SatelliteToolboxBase.jl) with a row-major ordering. If those matrices are
    not provided by the user, they will be created using that type of storage.

See also: [`gravitational_acceleration`](@ref)

# Arguments

- `model::AbstractGravityModel{T}`: Gravity model.
- `r::AbstractVector`: Position [m] in the body-fixed frame (ITRF for Earth) at which the
    acceleration is computed.
- `time::Union{Number, DateTime}`: Time at which the acceleration is computed, expressed
    as a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.
    (**Default**: J2000.0 epoch)

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available.
    (**Default**: -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than `max_degree`, it will be clamped.
    If it is lower than 0, it will be set to the same value as `max_degree`.
    (**Default**: -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default**: `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default**: `nothing`)
- `ω::Number`: Rotation rate of the body [rad/s]. For non-Earth bodies, provide the
    appropriate rotation rate for the celestial body.
    (**Default**: `EARTH_ANGULAR_SPEED`)

# Returns

- `SVector{3, RT}`: Gravity acceleration [m/s²] represented in the body-fixed frame (ITRF
    for Earth).

# References

- **[1]** Barthelmes, F (2013). *Definition of Functions of the Geopotential and Their
    Calculation from Spherical Harmonic Models*. Scientific Technical Report STR09/02.
    GeoForschungsZentrum (GFZ), pp. 22-23.
"""
function gravity_acceleration(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::Number = 0;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
    ω::Number = EARTH_ANGULAR_SPEED,
) where {T <: Number, V <: Number}

    # == Gravitational Acceleration ========================================================

    grav_itrf = gravitational_acceleration(
        model, r, time; max_degree = max_degree, max_order = max_order, P = P, dP = dP
    )

    # == Centripetal acceleration ==========================================================
    #
    # The centripetal acceleration at the desired position lies in a plane parallel to the
    # Equatorial plane and points toward the Earth's rotation axes.
    #
    # Applying the centrifugal potential in eq. 124 [1, p. 23] into eq. 120 [1, p. 22] and
    # then converting from the UEN reference frame to ITRF, one gets:
    #
    #                   ┌                          ┐                  ┌           ┐
    #                   │ ω² r cos(ϕ_gc) cos(λ_gc) │                  │ cos(λ_gc) │
    #   α_centrifugal = │ ω² r cos(ϕ_gc) sin(λ_gc) │ = ω² r cos(ϕ_gc) │ sin(λ_gc) │,
    #                   │            0             │                  │    0      │
    #                   └                          ┘                  └           ┘
    #
    # where ω is the Earth rotation rate, cos(λ_gc) = r_x / √(r_x² + r_y²), and
    # sin(λ_gc) = r_y / √(r_x² + r_y²).
    #
    # Since r ⋅ cos(ϕ_gc) = ρ_gc = √(r_x² + r_y²), the expression above simplifies to:
    #
    #   α_centrifugal = [ω² ⋅ r_x, ω² ⋅ r_y, 0].
    #
    # This form avoids the division by ρ_gc, which is 0 at the poles.
    #
    # NOTE: The Earth rotation axes is not aligned with the Z-axis of ITRF. However, this
    # approximation provides a sufficient accuracy for most applications.
    ω² = ω^2

    centrifugal_accel_itrf = @SVector [ω² * r[1], ω² * r[2], zero(T)]

    # Finally, compute the gravity acceleration.
    g_itrf = grav_itrf + centrifugal_accel_itrf

    return g_itrf
end

function gravity_acceleration(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
    ω::Number = EARTH_ANGULAR_SPEED,
) where {T <: Number, V <: Number}
    return gravity_acceleration(
        model,
        r,
        _to_j2000_seconds(time);
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
        ω = ω,
    )
end
