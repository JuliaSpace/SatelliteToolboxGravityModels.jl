## Description #############################################################################
#
# Function to compute the gravitational potential.
#
## References ##############################################################################
#
# [1] Barthelmes, F (2013). Definition of Functions of the Geopotential and Their
#     Calculation from Spherical Harmonic Models. Scientific Technical Report STR09/02.
#     GeoForschungsZentrum (GFZ).
#
############################################################################################

"""
    gravitational_potential(model::AbstractGravityModel, r::AbstractVector[, time]; kwargs...) -> RT

Compute the gravitational potential `U` [m²/s²] using the `model` in the position `r` [m],
represented in the body-fixed frame (ITRF for Earth), at instant `time`. If the latter
argument is omitted, the J2000.0 epoch (2000-01-01T12:00:00) is used.

The return type `RT` is obtained by promoting the type of the `model` coefficients, the
element type of `r`, and the type of `time`.

!!! note

    Gravitational potential is the potential caused by the central body mass only, i.e.,
    without considering the centrifugal potential.

!!! note

    The matrix `P` is lower triangular. Hence, the algorithm performance for large models
    can be improved if it is created using the `LowerTriangularStorage` (defined in
    SatelliteToolboxBase.jl) with a row-major ordering. If this matrix is not provided by
    the user, it will be created using that type of storage.

# Arguments

- `model::AbstractGravityModel{T}`: Gravity model.
- `r::AbstractVector`: Position [m] in the body-fixed frame (ITRF for Earth) at which the
    potential is computed.
- `time::Union{Number, DateTime}`: Time at which the potential is computed, expressed as a
    `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.
    (**Default**: J2000.0 epoch)

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational potential. If it is higher than the available number of coefficients in
    the `model`, it will be clamped. If it is lower than 0, it will be set to the maximum
    degree available.
    (**Default**: -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational potential. If it is higher than `max_degree`, it will be clamped. If it
    is lower than 0, it will be set to the same value as `max_degree`.
    (**Default**: -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default**: `nothing`)

# References

- **[1]** Barthelmes, F (2013). *Definition of Functions of the Geopotential and Their
    Calculation from Spherical Harmonic Models*. Scientific Technical Report STR09/02.
    GeoForschungsZentrum (GFZ), p. 19.
"""
function gravitational_potential(
    model::AbstractGravityModel{T},
    r::AbstractVector{V};
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number}
    return gravitational_potential(model, r, 0; max_degree, max_order, P)
end

function gravitational_potential(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::W;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number, W <: Number}
    RT = promote_type(T, V, W)

    n_max, m_max = _process_degree_and_order(model, max_degree, max_order)

    # Check if the matrix related to Legendre must be allocated.
    if isnothing(P)
        P = LowerTriangularStorage{RowMajor, RT}(n_max + 1)
    else
        _check_legendre_matrix(P, "P", n_max, m_max)
    end

    # Call the kernel through a function barrier. Hence, the hot loop is always compiled
    # with a concrete type for `P`, even when it is allocated here.
    return _gravitational_potential_kernel(model, r, time, n_max, m_max, P)
end

function gravitational_potential(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number}
    t = Dates.value(time - _DT_J2000) / 1000

    return gravitational_potential(
        model, r, t; max_degree = max_degree, max_order = max_order, P = P
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _gravitational_potential_kernel(
        model::AbstractGravityModel,
        r::AbstractVector,
        time::Number,
        n_max::Int,
        m_max::Int,
        P::AbstractMatrix
    ) -> RT

Compute the gravitational potential [m²/s²] of `model` at the position `r` [m],
represented in the body-fixed frame (ITRF for Earth), and instant `time`, expressed as the
number of elapsed seconds [s] from the J2000.0 epoch (2000-01-01T12:00:00), using the
spherical harmonics up to degree `n_max` and order `m_max`.

This function is the kernel of [`gravitational_potential`](@ref), called through a
function barrier so the hot loop is compiled with a concrete type for `P`. It assumes all
inputs were already processed: `n_max` and `m_max` must be valid for `model`, and `P` must
have at least `n_max + 1 × m_max + 1` elements, which are overwritten with the associated
Legendre function values.
"""
function _gravitational_potential_kernel(
    model::AbstractGravityModel{T},
    r::AbstractVector{V},
    time::W,
    n_max::Int,
    m_max::Int,
    P::AbstractMatrix,
) where {T <: Number, V <: Number, W <: Number}
    RT = promote_type(T, V, W)

    # == Unpack Gravity Model Data =========================================================

    μ = gravity_constant(model)
    R₀ = radius(model)
    norm = coefficient_norm(model)

    # == Geocentric Spherical Coordinates ==================================================

    r_gc, _, θ, sin_λ, cos_λ, south = _spherical_coordinates(r, RT)

    # == Auxiliary Variables ===============================================================

    # The Legendre functions are evaluated at the angle θ between the position and the
    # polar axis, which is the geocentric colatitude θ_gc = π / 2 - ϕ_gc in the northern
    # hemisphere and π - θ_gc in the southern one (see `_spherical_coordinates`). In the
    # latter case, the parity P_n,m[cos(π - θ)] = (-1)^(n + m) P_n,m[cos(θ)] is used to
    # recover the values at the colatitude: the factor (-1)^m is absorbed by shifting the
    # longitude by π, which changes the sign of its sine and cosine, and the factor (-1)^n
    # is absorbed by changing the sign of the ratio R₀ / r.
    if south
        sin_λ = -sin_λ
        cos_λ = -cos_λ
        ratio = -R₀ / r_gc
    else
        ratio = +R₀ / r_gc
    end

    fact = ratio

    # Sine and cosine of twice the geocentric longitude, which are used to initialize the
    # recursion that computes `sin(m * λ_gc)` and `cos(m * λ_gc)`.
    sin_2λ = 2sin_λ * cos_λ
    cos_2λ = cos_λ^2 - sin_λ^2

    # == Gravitational Potential ===========================================================

    U = RT(1)  # Gravitational potential

    # Compute the associated Legendre functions `P_n,m[cos(θ)]` with the required
    # normalization. Since θ ∈ [0, π / 2], the functions are well-defined and no sign
    # adjustments are needed.
    legendre!(norm, P, θ, n_max, m_max; ph_term = false)

    # Compute the potential.
    @inbounds for n in 2:n_max
        aux_U = RT(0)

        # == Sine and Cosine with m = 1 ====================================================
        #
        # These values will be used to update recursively `sin(m * λ_gc)` and
        # `cos(m * λ_gc)`, reducing the computational burden.
        sin_mλ   = RT(0)      # sin( 0 * λ_gc)
        sin_m_1λ = -sin_λ    # sin(-1 * λ_gc)
        sin_m_2λ = -sin_2λ   # sin(-2 * λ_gc)
        cos_mλ   = RT(1)      # cos( 0 * λ_gc)
        cos_m_1λ = +cos_λ    # cos(-1 * λ_gc)
        cos_m_2λ = +cos_2λ   # cos(-2 * λ_gc)

        # == Compute the Contributions When `m ∈ [0, min(n, m_max)]` =======================

        for m in 0:min(n, m_max)
            # Compute recursively `sin(m * λ_gc)` and `cos(m * λ_gc)`.
            sin_mλ = 2cos_λ * sin_m_1λ - sin_m_2λ
            cos_mλ = 2cos_λ * cos_m_1λ - cos_m_2λ

            # == Get the Spherical Harmonics Coefficients ==================================

            clm, slm = coefficients(model, n, m, time)

            CcSs_nm = clm * cos_mλ + slm * sin_mλ

            # == Compute the Contributions for `m` =========================================

            P_nm = P[n + 1, m + 1]

            aux_U += P_nm * CcSs_nm

            # == Update the Values for the Next Step =======================================

            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        # fact = (R₀ / r)^n
        fact *= ratio

        # aux_U *= (R₀ / r)^n
        aux_U *= fact

        U += aux_U
    end

    U *= μ / r_gc

    return U
end
