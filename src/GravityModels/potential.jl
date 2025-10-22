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
    gravitational_potential(model::AbstractGravityModel{Number, NormType}, r::AbstractVector{Number}[, time::Union{Number, DateTime}]; kwargs...) -> RT

Compute the gravitational potential [J / kg] or [m² / s²] using the `model` in the position
`r` [m], represented in ITRF, at instant `time`. If the latter argument is omitted, the
J2000.0 epoch is used (2000-01-01T12:00:00).

`time` can be expressed using a `DateTime` object or the number of ellapsed seconds from
J2000.0 epoch.

!!! note

    Gravitational potential is the potential caused by the central body mass only, i.e.,
    without considering the centrifugal potential.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    gravitational potential. If it is higher than the available number of coefficients in
    the `model`, it will be clamped. If it is lower than 0, it will be set to the maximum
    degree available.
    (**Default** = -1)
- `max_order::Int`: Maximum order used in the spherical harmonics when computing the
    gravitational potential. If it is higher than `max_degree`, it will be clamped. If it
    is lower than 0, it will be set to the same value as `max_degree`.
    (**Default** = -1)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)

!!! note

    The matrix `P` is lower triangular. Hence, the algorithm peformance for large models
    can be improved if it is created using the `LowerTriangularStorage` (defined in
    SatelliteToolboxBase.jl) with a row-major ordering. If this matrix is not provided by
    the user, it will be created using that type of storage.

# Returns

- `RT`: The gravitational potential `U`.
"""
function gravitational_potential(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V};
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}
    return gravitational_potential(model, r, 0; max_degree, max_order, P)
end

function gravitational_potential(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::W;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, W<:Number, NT<:Val}

    RT = promote_type(T, V, W)

    # == Unpack Gravity Model Data =========================================================

    μ  = gravity_constant(model)
    R₀ = radius(model)
    norm_type = coefficient_norm(model)
    model_max_degree = maximum_degree(model)

    # == Process the Inputs ================================================================

    # Check maximum degree value.
    if (max_degree < 0) || (max_degree > model_max_degree)
        max_degree = model_max_degree
    end

    # Check maximum order value.
    if (max_order < 0) || (max_order > max_degree)
        max_order = max_degree
    end

    # Obtain the degree and order used for the computation.
    n_max = max_degree
    m_max = max_order

    # Check if the matrix related to Legendre must be computed.
    if isnothing(P)
        P = LowerTriangularStorage{RowMajor, RT}(n_max + 1)
    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(P)

        if (rows < n_max + 1) || (cols < m_max + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max + 1) rows and $(m_max + 1) columns."))
        end
    end

    # == Geocentric Latitude and Longitude =================================================

    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc  + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc
    ϕ_gc  = atan(r[3], ρ_gc)
    λ_gc  = atan(r[2], r[1])

    # == Auxiliary Variables ===============================================================

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational burden.

    sin_λ,  cos_λ  = sincos(λ_gc)
    sin_2λ, cos_2λ = sincos(2λ_gc)

    # == Gravitational Potential ===========================================================

    U = RT(1)  # Gravitational potential

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the required
    # normalization.
    #
    # Notice that cos(ϕ_gc - π / 2) = sin(ϕ_gc).
    legendre!(Val(norm_type), P, ϕ_gc - RT(π / 2), n_max, m_max; ph_term = false)

    # Auxiliary variables.
    ratio = R₀ / r_gc
    fact  = ratio

    # Compute the potential.
    @inbounds for n in 2:n_max
        aux_U = RT(0)

        # == Sine and Cosine with m = 1 ====================================================
        #
        # This values will be used to update recursively `sin(m * λ_gc)` and
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

function gravitational_potential(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT<:Val}

    t = Dates.value(time - _DT_J2000) / 1000

    return gravitational_potential(
        model,
        r,
        t;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
    )
end

# Shorter alias for gravitational_potential
const potential = gravitational_potential
