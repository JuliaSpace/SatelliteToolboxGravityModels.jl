# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to compute the gravity model coefficients of a ICGEM file.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    icgem_coefficients(model::IcgemFile{T}, degree::Int, order::Int, t::DateTime) where T<:Number -> T, T

Compute the ICGEM coefficients (`Clm` and `Slm`) of the `model` for the specified `degree`
and `order` in the instant `t`.
"""
function icgem_coefficients(
    model::IcgemFile{T},
    degree::Int,
    order::Int,
    t::DateTime
) where T<:Number
    # First let's check if the degree and order is inside the expected range.
    order > degree && throw(ArgumentError("`order` must be lower than or equal to `degree`."))
    degree > model.max_degree &&
        throw(ArgumentError("The maximum degree available in the model is $(model.max_degree)."))

    # Get the data element related to the degree and order.
    coefficient = @inbounds model.data[degree + 1, order + 1]

    # If `coefficient` is `nothing`, we should return 0 because it was not defined.
    isnothing(coefficient) && return T(0), T(0)

    # Compute the `Clm` and `Slm` coefficients.
    clm, slm = _compute_icgem_coefficient(coefficient, t)

    return clm, slm
end

############################################################################################
#                                    Private Functions
############################################################################################

# Compute the coefficients `Clm` and `Slm` for a coefficient of type `IcgemGfcCoefficient`.
function _compute_icgem_coefficient(coefficient::IcgemGfcCoefficient, t::DateTime)
    return coefficient.clm, coefficient.slm
end

# Compute the coefficients `Clm` and `Slm` for a coefficient of type `IcgemGfctCoefficient`.
function _compute_icgem_coefficient(
    coefficient::IcgemGfctCoefficient{T},
    t::DateTime
) where T<:Number
    clm = coefficient.clm
    slm = coefficient.slm

    # Elapsed time from coefficients epoch [year].
    Δt = T(Dates.value(t - coefficient.time) / 1000 / 86400 / 365)

    # Trend
    # ======================================================================================

    if coefficient.has_trend
        clm += coefficient.trend_clm * Δt
        slm += coefficient.trend_slm * Δt
    end

    # asin
    # ======================================================================================

    for c in coefficient.asin_coefficients
        A_clm, A_slm, p = c

        aux  = sin(T(2π) / p * Δt)
        clm += A_clm * aux
        slm += A_slm * aux
    end

    # acos
    # ======================================================================================

    for c in coefficient.acos_coefficients
        A_clm, A_slm, p = c

        aux  = cos(T(2π) / p * Δt)
        clm += A_clm * aux
        slm += A_slm * aux
    end

    return clm, slm
end
