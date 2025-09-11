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
    time::T # .................................................. Seconds since J2000.0 epoch

    # == Trend =============================================================================

    has_trend::Bool
    trend_clm::T
    trend_slm::T

    # == asin ==============================================================================

    asin_coefficients::Vector{NTuple{3,T}}

    # == acos ==============================================================================

    acos_coefficients::Vector{NTuple{3,T}}
end

Base.length(c::IcgemGfcCoefficient) = 1
Base.length(c::IcgemGfctCoefficient) = 1
Base.iterate(c::IcgemGfcCoefficient) = (c, nothing)
Base.iterate(c::IcgemGfcCoefficient, state) = nothing
Base.iterate(c::IcgemGfctCoefficient) = (c, nothing)
Base.iterate(c::IcgemGfctCoefficient, state) = nothing

struct IcgemFile{T<:Number,NT<:Val} <: GravityModels.AbstractGravityModel{T,NT}
    # Fields related to the header.
    product_type::Symbol
    model_name::String
    gravity_constant::T
    radius::T
    max_degree::Int
    errors::Symbol
    tide_system::Symbol
    norm::NT

    # Fields related to the data section.
    data_static::Matrix{IcgemGfcCoefficient{T}}
    has_dynamic::Bool
    data_dynamic::Matrix{IcgemGfctCoefficient{T}}
end
