## Description #############################################################################
#
# Functions to compute the gravity model coefficients of a ICGEM file.
#
## References ##############################################################################
#
# [1] Barthelmes, F., Förste, C (2011). The ICGEM-format. GFZ Potsdam, Department 1
#     "Geodesy and Remote Sensing".
#
############################################################################################

"""
    icgem_coefficients(model::IcgemFile, degree::Int, order::Int, time) -> RT, RT

Compute the coefficients `Clm` and `Slm` [-] of the ICGEM `model` for the specified
`degree` and `order` at the instant `time`, expressed as a `DateTime` object or the number
of elapsed seconds [s] from the J2000.0 epoch (2000-01-01T12:00:00).

The function throws an `ArgumentError` if `order` is higher than `degree` or if `degree`
is higher than the maximum degree available in `model`.

The return type `RT` is `float(promote_type(T, typeof(time)))`, where `T` is the type of
the `model` coefficients.

# Arguments

- `model::IcgemFile{T}`: ICGEM model.
- `degree::Int`: Degree of the coefficients.
- `order::Int`: Order of the coefficients.
- `time::Union{Number, DateTime}`: Time at which the coefficients are computed, expressed
    as a `DateTime` object or the number of elapsed seconds [s] from the J2000.0 epoch.

# Returns

- `RT`: Coefficient `Clm` [-].
- `RT`: Coefficient `Slm` [-].
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

    # Promote the coefficients beforehand so both return paths have the same type.
    RT = float(promote_type(T, typeof(time)))

    # Check if the coefficient is time-variable.
    if degree <= model.max_time_variable_degree
        k = @inbounds model.time_variable_index[degree + 1, order + 1]

        (k != 0) && return _compute_icgem_time_variable_coefficient(
            model.time_variable_coefficients, Int(k), time, RT
        )
    end

    # Otherwise, return the constant coefficient.
    coefficient = @inbounds model.data[degree + 1, order + 1]

    return RT(coefficient.clm), RT(coefficient.slm)
end

function icgem_coefficients(
    model::IcgemFile{T}, degree::Int, order::Int, time::DateTime
) where {T <: Number}
    return icgem_coefficients(model, degree, order, _to_j2000_seconds(time))
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _compute_icgem_time_variable_coefficient(coefficients::Vector{IcgemTimeVariableCoefficient{T}}, k::Int, t::Number, ::Type{RT}) -> RT, RT

Compute the coefficients `Clm` [-] and `Slm` [-] of the time-variable coefficient at the
instant `t`, expressed as the number of elapsed seconds [s] from the J2000.0 epoch
(2000-01-01T12:00:00), using the element type `RT`.

`k` must be the index in `coefficients` of the first object related to the desired degree
and order. The function selects the object whose validity interval contains `t`, assuming
that the objects related to the same degree and order are stored consecutively and sorted
by epoch. If `t` is before the first epoch, the first object is used, and if `t` is after
the last epoch, the last object is used.

The coefficients are obtained by adding the linear trend and the periodic terms to the
values at the coefficient epoch, as described in the ICGEM format documentation [1]. The
elapsed time from the epoch is converted to Julian years (365.25 days).

# References

- **[1]** Barthelmes, F., Förste, C (2011). *The ICGEM-format*. GFZ Potsdam, Department 1
    "Geodesy and Remote Sensing".
"""
function _compute_icgem_time_variable_coefficient(
    coefficients::Vector{IcgemTimeVariableCoefficient{T}}, k::Int, t::Number, ::Type{RT}
) where {T <: Number, RT}
    @inbounds begin
        c = coefficients[k]

        # Select the last validity interval that starts before or at `t`.
        while (k < length(coefficients))
            c_next = coefficients[k + 1]

            ((c_next.degree != c.degree) || (c_next.order != c.order)) && break
            (t < c_next.t₀) && break

            k += 1
            c = c_next
        end
    end

    # Elapsed time from coefficients epoch [year], considering a Julian year with 365.25
    # days.
    Δt = RT((t - c.t₀) / (86400 * 365.25))

    # == Trend =============================================================================

    clm = RT(c.clm) + c.trend_clm * Δt
    slm = RT(c.slm) + c.trend_slm * Δt

    # == Periodic Terms ====================================================================

    for p in c.periodic_terms
        sin_ωt, cos_ωt = sincos(RT(2π) * Δt / p.period)

        clm += p.amplitude_sin_clm * sin_ωt + p.amplitude_cos_clm * cos_ωt
        slm += p.amplitude_sin_slm * sin_ωt + p.amplitude_cos_slm * cos_ωt
    end

    return clm, slm
end
