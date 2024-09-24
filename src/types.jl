## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile

############################################################################################
#                                          ICGEM                                           #
############################################################################################

struct IcgemGfcCoefficient{T<:Number}
    clm::T
    slm::T
end

struct IcgemGfctCoefficient{T<:Number}
    clm::T
    slm::T
    time::Number # Seconds since JD_J2000

    # == Trend =============================================================================

    has_trend::Bool
    trend_clm::T
    trend_slm::T

    # == asin ==============================================================================

    asin_coefficients::Vector{NTuple{3, T}}

    # == acos ==============================================================================

    acos_coefficients::Vector{NTuple{3, T}}
end

struct IcgemFile{T<:Number} <: GravityModels.AbstractGravityModel{T}
    # Fields related to the header.
    product_type::Symbol
    model_name::String
    gravity_constant::T
    radius::T
    max_degree::Int
    errors::Symbol
    tide_system::Symbol
    norm::Symbol

    # Fields related to the data section.
    data::Matrix{Union{Nothing, IcgemGfcCoefficient{T}, IcgemGfctCoefficient{T}}}
end
