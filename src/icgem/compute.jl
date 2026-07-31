## Description #############################################################################
#
# Functions to compute the gravity model coefficients of a ICGEM file.
#
############################################################################################

"""
    icgem_coefficients(model::IcgemFile, degree::Int, order::Int, time) -> RT, RT

Compute the coefficients `Clm` and `Slm` [-] of the ICGEM `model` for the specified
`degree` and `order` at the instant `time`, expressed as a `DateTime` object or the number
of elapsed seconds [s] from the J2000.0 epoch (2000-01-01T12:00:00).

The function throws an `ArgumentError` if `order` is higher than `degree` or if `degree`
is higher than the maximum degree available in `model`.

# Arguments

- `model::IcgemFile{T}`: ICGEM model.
- `degree::Int`: Degree of the coefficients.
- `order::Int`: Order of the coefficients.
- `time::Union{Number, DateTime}`: Time at which the coefficients are computed, expressed
    as a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.

# Returns

- `RT`: Coefficient `Clm` [-].
- `RT`: Coefficient `Slm` [-].

The return type `RT` is `T` for models with only constant coefficients, or
`float(promote_type(T, typeof(time)))` for models with time-variable coefficients.
"""
function icgem_coefficients(
    model::IcgemFile{T}, degree::Int, order::Int, time::Number
) where {T <: Number}
    # First let's check if the degree and order is inside the expected range.
    order > degree &&
        throw(ArgumentError("`order` must be lower than or equal to `degree`."))
    degree > model.max_degree && throw(
        ArgumentError("The maximum degree available in the model is $(model.max_degree)."),
    )

    # Get the data element related to the degree and order.
    coefficient = @inbounds model.data[degree + 1, order + 1]
    return _compute_icgem_coefficient(coefficient, time)
end

function icgem_coefficients(
    model::IcgemFile{T}, degree::Int, order::Int, time::DateTime
) where {T <: Number}
    t = Dates.value(time - _DT_J2000) / 1000

    return icgem_coefficients(model, degree, order, t)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _compute_icgem_coefficient(coefficient::IcgemGfcCoefficient{T}, t::Number) -> T, T

Return the constant coefficients `Clm` [-] and `Slm` [-] stored in `coefficient`. The time
`t` [s] is unused since the coefficient is constant.
"""
function _compute_icgem_coefficient(coefficient::IcgemGfcCoefficient, t::Number)
    return coefficient.clm, coefficient.slm
end

"""
    _compute_icgem_coefficient(coefficient::IcgemGfctCoefficient{T}, t::Number) -> RT, RT

Compute the coefficients `Clm` [-] and `Slm` [-] of the time-variable `coefficient` at the
instant `t`, expressed as the number of elapsed seconds [s] from the J2000.0 epoch
(2000-01-01T12:00:00).

The coefficients are obtained by adding the linear trend and the sine and cosine periodic
terms to the values at the coefficient epoch, as described in the ICGEM format
documentation [1]. The elapsed time from the epoch is converted to Julian years (365.25
days).

The return type `RT` is `float(promote_type(T, typeof(t)))`.

# References

- **[1]** Barthelmes, F., Förste, C (2011). *The ICGEM-format*. GFZ Potsdam, Department 1
    "Geodesy and Remote Sensing".
"""
function _compute_icgem_coefficient(
    coefficient::IcgemGfctCoefficient{T}, t::Number
) where {T <: Number}
    # Promote the coefficients beforehand so both return paths have the same type.
    RT = float(promote_type(T, typeof(t)))

    clm = RT(coefficient.clm)
    slm = RT(coefficient.slm)

    coefficient.is_time_varying || return clm, slm

    # Elapsed time from coefficients epoch [year], considering a Julian year with 365.25
    # days.
    Δt = RT((t - coefficient.time) / 86400 / 365.25)

    # == Trend =============================================================================

    if coefficient.has_trend
        clm += coefficient.trend_clm * Δt
        slm += coefficient.trend_slm * Δt
    end

    # == asin ==============================================================================

    for c in coefficient.asin_coefficients
        A_clm, A_slm, p = c

        aux = sin(RT(2π) / p * Δt)
        clm += A_clm * aux
        slm += A_slm * aux
    end

    # == acos ==============================================================================

    for c in coefficient.acos_coefficients
        A_clm, A_slm, p = c

        aux = cos(RT(2π) / p * Δt)
        clm += A_clm * aux
        slm += A_slm * aux
    end

    return clm, slm
end
