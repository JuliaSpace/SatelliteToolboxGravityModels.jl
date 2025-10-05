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
    gravitational_acceleration(model::AbstractGravityModel{Number, NormType}, r::AbstractVector{Number}[, time::Union{Number, DateTime}]; kwargs...) -> NTuple{3, RT}

Compute the gravitational acceleration [m / s²] represented in ITRF using the `model` in the
position `r` [m], also represented in ITRF, at instant `time`. If the latter argument is
omitted, the J2000.0 epoch is used (2000-01-01T12:00:00).

`time` can be expressed using a `DateTime` object or the number of ellapsed seconds from
J2000.0 epoch.

!!! note

    Gravitational acceleration is the acceleration caused by the central body mass only,
    i.e., without considering the centrifugal potential.

# Keywords

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

# Returns

- `RT`: The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
- `RT`: The derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`).
- `RT`: The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).
"""
function gravitational_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V};
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}
    return gravitational_acceleration(model, r, 0; max_degree, max_order, P, dP)
end

function gravitational_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}

    # Compute the partial derivatives of the gravitational field w.r.t. the spherical
    # coordinates.
    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ = gravitational_field_derivative(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP
    )

    # Auxiliary variables.
    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc  + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc
    ϕ_gc  = atan(r[3], ρ_gc)
    λ_gc  = atan(r[2], r[1])

    # == Acceleration Represented in the ITRF ==============================================

    # Compute the partial derivatives in spherical coordinate systems [1, p. 22] (eq. 120):
    #
    #     ∂U        1      ∂U     1   ∂U
    #    ---- , ---------.---- , ---.----
    #     ∂r     r.cos ϕ   ∂λ     r   ∂ϕ
    #
    # Notice that the singularity is not a problem here. When computing `cos(π / 2)` a very
    # small number will be returned and `∂U / ∂λ` is 0. Hence, the 2nd component will be 0.

    a_uen = @SVector [
        ∂U_∂r,
        ∂U_∂λ / (r_gc * cos(ϕ_gc)),
        ∂U_∂ϕ / r_gc
    ]

    # The vector `a_uen` is represented in the local UEN (Up-Earth-North) reference frame.
    # Hence, we need to describe the unitary vectors of this frame in the ECEF reference
    # frame. This can be accomplished by the following rotations matrix.
    D_itrf_uen = angle_to_dcm(ϕ_gc, -λ_gc, :YZ)
    a_itrf = D_itrf_uen * a_uen

    return a_itrf
end

function gravitational_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}

    t = Dates.value(time - _DT_J2000) / 1000

    return gravitational_acceleration(
        model,
        r,
        t;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
    )
end

"""
    gravity_acceleration(model::AbstractGravityModel{Number, NormType}, r::AbstractVector{Number}[, time::Union{Number, DataTime}]; kwargs...) -> NTuple{3, RT}

Compute the gravity acceleration [m / s²] represented in ITRF using the `model` in the
position `r` [m], also represented in ITRF, at instant `time`. If the latter argument is
omitted, the J2000.0 epoch is used.

`time` can be expressed using a `DateTime` object or the number of ellapsed seconds from
J2000.0 epoch.

!!! note

    Gravity acceleration is the compound acceleration caused by the central body mass and
    the centrifugal force due to the planet's rotation.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational field derivative. If it is higher than the available number of
    coefficients in the `model`, it will be clamped. If it is lower than 0, it will be set
    to the maximum degree available. (**Default** = -1)
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

# Returns

- `RT`: The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
- `RT`: The derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`).
- `RT`: The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).
"""
function gravity_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V};
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number,V<:Number, NT<:Val}
    return gravity_acceleration(model, r, 0; max_degree, max_order, P, dP)
end

function gravity_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}

    # == Gravitational Acceleration ========================================================

    grav_itrf = gravitational_acceleration(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP
    )

    # == Centripetal acceleration ==========================================================
    #
    # The centripetal acceleration has the following value:
    #
    #   cp_accel = ω² r cos(ϕ_gc),
    #
    # where `ω` is the Earth's rotation rate, `r` is the distance from the Earth's center,
    # and `ϕ_gc` is the geocentric latitude. Notice that:
    #
    #   cos(ϕ_gc) = √(r_x² + r_y²) / √(r_x² + r_y² + r_z²),
    #
    # and:
    #
    #   r = √(r_x² + r_y² + r_z²).
    #
    # Hence, letting `ρ_gc = √(r_x² + r_y²)` one gets:
    #
    #   r ⋅ cos(ϕ_gc) = ρ_gc
    #
    # Finally:
    #
    #   cp_accel = ω² ρ_gc

    ρ²_gc    = r[1]^2 + r[2]^2
    ρ_gc     = √ρ²_gc
    cp_accel = EARTH_ANGULAR_SPEED^2 * ρ_gc

    # The centripetal acceleration at the desired position lies in a plane parallel to the
    # Equatorial plane and points toward the Earth's rotation axes.
    #
    # Applying the centrifugal potential in eq. 124 [1, p. 23] into eq. 120 [1, p. 22] and
    # then converting from the UEN reference frame to ITRF, one gets:
    #
    #                   ┌                          ┐                  ┌           ┐
    #                   │ ω² r cos(ϕ_gc) cos(λ_gc) │                  │ cos(λ_gc) │
    #   α_centrifugal = │            0             │ = ω² r cos(ϕ_gc) │    0      │,
    #                   │ ω² r cos(ϕ_gc) sin(λ_gc) │                  │ sin(λ_gc) │
    #                   └                          ┘                  └           ┘
    #
    # where ω is the Earth rotation rate, cos(λ_gc) = r_x / √(r_x² + r_y²), and
    # sin(λ_gc) = r_y / √(r_x² + r_y²).
    #
    # NOTE: The Earth rotation axes is not aligned with the Z-axis of ITRF. However, this
    # approximation provides a sufficient accuracy for most applications.
    centrifugal_accel_itrf = @SVector [
        r[1] / ρ_gc * cp_accel,
        r[2] / ρ_gc * cp_accel,
        T(0)
    ]

    # Finally, compute the gravity acceleration.
    g_itrf = grav_itrf + centrifugal_accel_itrf

    return g_itrf
end

function gravity_acceleration(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}

    t = Dates.value(time - _DT_J2000) / 1000

    return gravity_acceleration(
        model,
        r,
        t;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
    )
end
