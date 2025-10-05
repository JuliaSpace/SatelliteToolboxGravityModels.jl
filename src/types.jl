## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile

abstract type AbstractIcgemCoefficient{T<:Number} end

############################################################################################
#                                          ICGEM                                           #
############################################################################################

struct IcgemGfcCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}
    clm::T
    slm::T
end

struct IcgemGfctCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}
    clm::T
    slm::T
    time::T # .................................................. Seconds since J2000.0 epoch

    # This variable indicates if the coefficient is time-varying. If false, the time
    # variable is ignored and the coefficient is treated as a regular `IcgemGfcCoefficient`.
    # This is useful to avoid unnecessary computations, leading to a huge performance boost.
    is_time_varying::Bool

    # == Trend =============================================================================

    has_trend::Bool
    trend_clm::T
    trend_slm::T

    # == asin ==============================================================================

    asin_coefficients::Vector{NTuple{3, T}}

    # == acos ==============================================================================

    acos_coefficients::Vector{NTuple{3, T}}
end

# Constructor for `IcgemGfctCoefficient` from `IcgemGfcCoefficient`.
IcgemGfctCoefficient(c::IcgemGfcCoefficient{T}) where T = IcgemGfctCoefficient(
    c.clm,
    c.slm,
    zero(T),
    false,
    false,
    zero(T),
    zero(T),
    NTuple{3,T}[],
    NTuple{3,T}[],
)

function Base.zero(::Type{IcgemGfcCoefficient{T}}) where T
    return IcgemGfcCoefficient(zero(T), zero(T))
end

function Base.zero(::Type{IcgemGfctCoefficient{T}}) where T
    return IcgemGfctCoefficient(
        zero(T),
        zero(T),
        zero(T),
        false,
        false,
        zero(T),
        zero(T),
        Vector{NTuple{3,T}}(),
        Vector{NTuple{3,T}}()
    )
end

struct IcgemFile{
    T<:Number,
    NT<:Val,
    Coeff<:AbstractIcgemCoefficient{T}
} <: GravityModels.AbstractGravityModel{T, NT}
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
    data::LowerTriangularStorage{RowMajor, Coeff}
end
