## Description #############################################################################
#
# Function to compute the gravitational field derivative.
#
############################################################################################

"""
    gravitational_field_derivative(model::AbstractGravityModel, r::AbstractVector[, time]; kwargs...) -> NTuple{3, RT}

Compute the gravitational field derivative with respect to the spherical coordinates
(`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) using the `model` in the position `r` [m], represented in the
body-fixed frame (ITRF for Earth), at instant `time`. If the latter argument is omitted,
the J2000.0 epoch (2000-01-01T12:00:00) is used.

The return element type `RT` is obtained by promoting the type of the `model`
coefficients, the element type of `r`, and the type of `time`.

!!! info

    In this case, `ϕ` is the geocentric latitude and `λ` is the longitude.

!!! note

    The matrices `P` and `dP` are lower triangular. Hence, the algorithm performance for
    large models can be improved if they are created using the `LowerTriangularStorage`
    (defined in SatelliteToolboxBase.jl) with a row-major ordering. If those matrices are
    not provided by the user, they will be created using that type of storage.

# Arguments

- `model::AbstractGravityModel{T, NT}`: Gravity model.
- `r::AbstractVector`: Position [m] in the body-fixed frame (ITRF for Earth) at which the
    derivative is computed.
- `time::Union{Number, DateTime}`: Time at which the derivative is computed, expressed as
    a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.
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

- `RT`: Derivative of the gravitational field w.r.t. the radius (`∂U/∂r`) [m/s²].
- `RT`: Derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`)
    [m²/s²].
- `RT`: Derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`) [m²/s²].
"""
function gravitational_field_derivative(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V};
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number, NT <: Val}
    return gravitational_field_derivative(model, r, 0; max_degree, max_order, P, dP)
end

function gravitational_field_derivative(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::W;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number, W <: Number, NT <: Val}
    RT = promote_type(T, V, W)

    n_max, m_max, n_max_P, m_max_P, n_max_dP, m_max_dP, P, dP = _prepare_field_derivative_inputs(
        model, RT, max_degree, max_order, P, dP
    )

    # Call the kernel through a function barrier. Hence, the hot loop is always compiled
    # with concrete types for `P` and `dP`, even when they are allocated here.
    ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ, _ = _gravitational_field_derivative_kernel(
        model, r, time, n_max, m_max, n_max_P, m_max_P, n_max_dP, m_max_dP, P, dP
    )

    return ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ
end

function gravitational_field_derivative(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::DateTime;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T <: Number, V <: Number, NT <: Val}
    t = Dates.value(time - _DT_J2000) / 1000

    return gravitational_field_derivative(
        model, r, t; max_degree = max_degree, max_order = max_order, P = P, dP = dP
    )
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _process_degree_and_order(model::AbstractGravityModel, max_degree::Int, max_order::Int) -> Int, Int

Return the maximum degree and order used in the spherical harmonics computation of `model`
given the requested `max_degree` and `max_order`. If `max_degree` is negative or higher
than the maximum degree of `model`, it is clamped to the latter. If `max_order` is negative
or higher than the selected degree, it is set to the selected degree.

# Returns

- `Int`: Maximum degree used in the computation.
- `Int`: Maximum order used in the computation.
"""
function _process_degree_and_order(
    model::AbstractGravityModel, max_degree::Int, max_order::Int
)
    model_max_degree = maximum_degree(model)

    n_max =
        ((max_degree < 0) || (max_degree > model_max_degree)) ? model_max_degree :
        max_degree
    m_max = ((max_order < 0) || (max_order > n_max)) ? n_max : max_order

    return n_max, m_max
end

"""
    _check_legendre_matrix(M::AbstractMatrix, name::String, n_max::Int, m_max::Int) -> Nothing

Check if the matrix `M`, called `name` in the error messages, has at least `n_max + 1` rows
and `m_max + 1` columns, throwing an `ArgumentError` otherwise.
"""
function _check_legendre_matrix(M::AbstractMatrix, name::String, n_max::Int, m_max::Int)
    rows, cols = size(M)

    if (rows < n_max + 1) || (cols < m_max + 1)
        throw(
            ArgumentError(
                "Matrix `$name` must have at least $(n_max + 1) rows and $(m_max + 1) columns.",
            ),
        )
    end

    return nothing
end

"""
    _prepare_field_derivative_inputs(model::AbstractGravityModel, RT::Type, max_degree::Int, max_order::Int, P::Union{Nothing, AbstractMatrix}, dP::Union{Nothing, AbstractMatrix}) -> Int, Int, Int, Int, Int, Int, AbstractMatrix, AbstractMatrix

Process the inputs of the gravitational field derivative computation of `model` with
element type `RT`, returning the degrees and orders used in the computation and the
matrices `P` and `dP` to store the associated Legendre functions and their derivatives.

The requested `max_degree` and `max_order` are clamped as described in
[`_process_degree_and_order`](@ref). If `P` or `dP` is `nothing`, the corresponding matrix
is allocated using a `LowerTriangularStorage` with row-major ordering and element type
`RT`. Otherwise, the function throws an `ArgumentError` if the matrix is too small.

# Returns

- `Int`: Maximum degree `n_max` used in the computation.
- `Int`: Maximum order `m_max` used in the computation.
- `Int`: Maximum degree computed in `P`.
- `Int`: Maximum order computed in `P`, which is `m_max + 1` if `m_max < n_max` because
    the derivative computation requires one additional order.
- `Int`: Maximum degree computed in `dP`.
- `Int`: Maximum order computed in `dP`.
- `AbstractMatrix`: Matrix `P`.
- `AbstractMatrix`: Matrix `dP`.
"""
function _prepare_field_derivative_inputs(
    model::AbstractGravityModel,
    ::Type{RT},
    max_degree::Int,
    max_order::Int,
    P::Union{Nothing, AbstractMatrix},
    dP::Union{Nothing, AbstractMatrix},
) where {RT}
    n_max, m_max = _process_degree_and_order(model, max_degree, max_order)

    # Obtain the required sizes for the matrices P and dP.
    #
    # Notice that, to compute the derivative if `m_max < n_max`, we need that `P` has an
    # order at least one time higher than `dP`. Otherwise, we will access regions with
    # undefined numbers.
    n_max_P  = n_max
    m_max_P  = (n_max == m_max) ? m_max : m_max + 1
    n_max_dP = n_max
    m_max_dP = m_max

    # Check if the matrices related to Legendre must be allocated.
    if isnothing(P)
        P = LowerTriangularStorage{RowMajor, RT}(n_max_P + 1)
    else
        _check_legendre_matrix(P, "P", n_max_P, m_max_P)
    end

    if isnothing(dP)
        dP = LowerTriangularStorage{RowMajor, RT}(n_max_dP + 1)
    else
        _check_legendre_matrix(dP, "dP", n_max_dP, m_max_dP)
    end

    return n_max, m_max, n_max_P, m_max_P, n_max_dP, m_max_dP, P, dP
end

"""
    _spherical_coordinates(r::AbstractVector, ::Type{RT}) -> RT, RT, RT, RT, RT, Bool

Compute the geocentric spherical coordinates of the position `r` [m] using the element type
`RT`.

The sine and cosine of the longitude are obtained directly from the coordinates, avoiding
trigonometric functions. On the polar axis, where the longitude is undefined, the longitude
is defined as 0.

The returned angle is the geocentric colatitude folded to the northern hemisphere, i.e. the
angle between the position and the polar axis in the interval `[0, π / 2]`. The last
returned value indicates whether the position lies in the southern hemisphere, in which
case the geocentric colatitude is `π` minus the returned angle. The angle is computed
directly from the coordinates instead of subtracting the latitude from `π / 2`, avoiding
the catastrophic cancellation near the poles, where the latitude rounds to `π / 2`. The
folding is required because the floating-point numbers near `π` cannot represent the
colatitude near the south pole as accurately as those near `0` represent it near the north
pole.

# Returns

- `RT`: Distance from the origin [m].
- `RT`: Distance from the polar axis [m].
- `RT`: Geocentric colatitude folded to the northern hemisphere [rad] in the interval
    `[0, π / 2]`.
- `RT`: Sine of the longitude [-].
- `RT`: Cosine of the longitude [-].
- `Bool`: `true` if the position lies in the southern hemisphere, `false` otherwise.
"""
function _spherical_coordinates(r::AbstractVector, ::Type{RT}) where {RT}
    ρ²_gc = r[1]^2 + r[2]^2
    r²_gc = ρ²_gc + r[3]^2
    r_gc  = √r²_gc
    ρ_gc  = √ρ²_gc

    south = r[3] < 0
    θ = atan(ρ_gc, abs(r[3]))

    if ρ_gc > 0
        sin_λ = RT(r[2] / ρ_gc)
        cos_λ = RT(r[1] / ρ_gc)
    else
        sin_λ = zero(RT)
        cos_λ = one(RT)
    end

    return RT(r_gc), RT(ρ_gc), RT(θ), sin_λ, cos_λ, south
end

"""
    _gravitational_field_derivative_kernel(
        model::AbstractGravityModel,
        r::AbstractVector,
        time::Number,
        n_max::Int,
        m_max::Int,
        n_max_P::Int,
        m_max_P::Int,
        n_max_dP::Int,
        m_max_dP::Int,
        P::AbstractMatrix,
        dP::AbstractMatrix
    ) -> NTuple{4, RT}

Compute the derivative of the gravitational field of `model` with respect to the spherical
coordinates at the position `r` [m], represented in the body-fixed frame (ITRF for Earth),
and instant `time`, expressed as the number of elapsed seconds [s] from the J2000.0 epoch
(2000-01-01T12:00:00), using the spherical harmonics up to degree `n_max` and order
`m_max`.

This function is the kernel of [`gravitational_field_derivative`](@ref), called through a
function barrier so the hot loop is compiled with concrete types for `P` and `dP`. It
assumes all inputs were already processed: `n_max` and `m_max` must be valid for `model`,
and `P` and `dP` must have at least `n_max_P + 1 × m_max_P + 1` and
`n_max_dP + 1 × m_max_dP + 1` elements, respectively, which are overwritten with the
associated Legendre function values and their derivatives.

# Returns

- `RT`: Derivative of the gravitational field w.r.t. the radius (`∂U/∂r`) [m/s²].
- `RT`: Derivative of the gravitational field w.r.t. the geocentric latitude (`∂U/∂ϕ`)
    [m²/s²].
- `RT`: Derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`) [m²/s²].
- `RT`: Limit of `(∂U/∂λ) / cos(ϕ)` [m²/s²] at the poles, which is required to compute the
    east component of the acceleration on the polar axis, where both `∂U/∂λ` and `cos(ϕ)`
    are 0. This value is meaningful only if `r` lies on the polar axis.
"""
function _gravitational_field_derivative_kernel(
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::W,
    n_max::Int,
    m_max::Int,
    n_max_P::Int,
    m_max_P::Int,
    n_max_dP::Int,
    m_max_dP::Int,
    P::AbstractMatrix,
    dP::AbstractMatrix,
) where {T <: Number, V <: Number, W <: Number, NT <: Val}
    RT = promote_type(T, V, W)

    # == Unpack Gravity Model Data =========================================================

    μ = gravity_constant(model)
    R₀ = radius(model)
    norm_type = coefficient_norm(model)

    # == Geocentric Spherical Coordinates ==================================================

    r_gc, _, θ, sin_λ, cos_λ, south = _spherical_coordinates(r, RT)

    # == Auxiliary Variables ===============================================================

    # The Legendre functions are evaluated at the angle θ between the position and the
    # polar axis, which is the geocentric colatitude θ_gc = π / 2 - ϕ_gc in the northern
    # hemisphere and π - θ_gc in the southern one (see `_spherical_coordinates`). In the
    # latter case, the parity P_n,m[cos(π - θ)] = (-1)^(n + m) P_n,m[cos(θ)] is used to
    # recover the values at the colatitude:
    #
    #   - The factor (-1)^m is absorbed by shifting the longitude by π, which changes the
    #     sign of its sine and cosine;
    #   - The factor (-1)^n is absorbed by changing the sign of the ratio R₀ / r; and
    #   - The derivatives w.r.t. θ_gc are the derivatives w.r.t. θ with the sign changed,
    #     which is applied at the end.
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

    # == First Derivative of the Non-Spherical Portion of the Gravitational Field ==========

    ∂U_∂r = RT(1)  # ........................................... Derivative w.r.t. the radius
    ∂U_∂ϕ = RT(0)  # .............................. Derivative w.r.t. the geocentric latitude
    ∂U_∂λ = RT(0)  # ............................. Derivative w.r.t. the geocentric longitude
    ∂U_∂λ_over_cosϕ_pole = RT(0)  # ........... Limit of (∂U/∂λ) / cos(ϕ_gc) at the poles

    # Compute the associated Legendre functions `P_n,m[cos(θ)]` with the required
    # normalization and their first-order derivatives w.r.t. θ. Since θ ∈ [0, π / 2], the
    # functions are well-defined and no sign adjustments are needed.
    legendre!(Val(norm_type), P, θ, n_max_P, m_max_P; ph_term = false)
    dlegendre!(Val(norm_type), dP, θ, P, n_max_dP, m_max_dP; ph_term = false)

    # Compute the derivatives.
    @inbounds for n in 2:n_max
        aux_∂U_∂r = RT(0)
        aux_∂U_∂ϕ = RT(0)
        aux_∂U_∂λ = RT(0)
        aux_∂U_∂λ_over_cosϕ_pole = RT(0)

        # == Sine and Cosine with m = 1 ====================================================
        #
        # These values will be used to update recursively `sin(m * λ_gc)` and
        # `cos(m * λ_gc)`, reducing the computational burden.
        #
        # TODO: Cache the computation.
        # We tried to compute those values only once using an external vector to store the
        # values. However, it leads to a worst performance. This behavior need further
        # investigation.
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
            ScCs_nm = slm * cos_mλ - clm * sin_mλ

            # == Compute the Contributions for `m` =========================================

            P_nm  = P[n + 1, m + 1]
            dP_nm = dP[n + 1, m + 1]

            aux_∂U_∂r += P_nm * CcSs_nm
            aux_∂U_∂ϕ += dP_nm * CcSs_nm
            aux_∂U_∂λ += m * P_nm * ScCs_nm

            # At the poles, `P_n,m / cos(ϕ_gc)` vanishes for `m > 1` and tends to the
            # derivative of `P_n,1` for `m = 1`.
            (m == 1) && (aux_∂U_∂λ_over_cosϕ_pole += dP_nm * ScCs_nm)

            # == Update the Values for the Next Step =======================================

            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        # fact = (a / r)^n
        fact *= ratio

        # aux_<> *= (a / r)^n
        aux_∂U_∂r *= fact
        aux_∂U_∂ϕ *= fact
        aux_∂U_∂λ *= fact
        aux_∂U_∂λ_over_cosϕ_pole *= fact

        ∂U_∂r += (n + 1) * aux_∂U_∂r
        ∂U_∂ϕ += aux_∂U_∂ϕ
        ∂U_∂λ += aux_∂U_∂λ
        ∂U_∂λ_over_cosϕ_pole += aux_∂U_∂λ_over_cosϕ_pole
    end

    # The term `∂U_∂ϕ` was computed with the derivatives w.r.t. θ, which is θ_gc in the
    # northern hemisphere and π - θ_gc in the southern one. Since ϕ_gc = π / 2 - θ_gc, the
    # derivative w.r.t. ϕ_gc has the sign changed in the northern hemisphere only.
    #
    # The limit of `P_n,1[cos(θ)] / sin(θ)` when θ → 0 is the derivative of `P_n,1` w.r.t. θ
    # at 0. Hence, the term `∂U_∂λ_over_cosϕ_pole` does not require a sign change in any
    # hemisphere, since the parity factors were already absorbed by the longitude shift and
    # by the sign of the ratio R₀ / r.
    ∂U_∂r *= -μ / r_gc^2
    ∂U_∂ϕ *= south ? +μ / r_gc : -μ / r_gc
    ∂U_∂λ *= +μ / r_gc
    ∂U_∂λ_over_cosϕ_pole *= +μ / r_gc

    return ∂U_∂r, ∂U_∂ϕ, ∂U_∂λ, ∂U_∂λ_over_cosϕ_pole
end
