## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile

############################################################################################
#                                          ICGEM                                           #
############################################################################################

abstract type AbstractIcgemCoefficient{T<:Number} end

struct IcgemGfcCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}
    clm::T
    slm::T
end

struct IcgemGfctCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}
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

IcgemGfctCoefficient(c::IcgemGfcCoefficient{T}) where T = IcgemGfctCoefficient(
    c.clm,
    c.slm,
    zero(T),
    false,
    zero(T),
    zero(T),
    NTuple{3,T}[],
    NTuple{3,T}[],
)


Base.length(c::IcgemGfcCoefficient) = 1
Base.length(c::IcgemGfctCoefficient) = 1
Base.iterate(c::IcgemGfcCoefficient) = (c, nothing)
Base.iterate(c::IcgemGfcCoefficient, state) = nothing
Base.iterate(c::IcgemGfctCoefficient) = (c, nothing)
Base.iterate(c::IcgemGfctCoefficient, state) = nothing
Base.zero(::Type{IcgemGfcCoefficient{T}}) where T = IcgemGfcCoefficient(zero(T), zero(T))
Base.zero(::Type{IcgemGfctCoefficient{T}}) where T = IcgemGfctCoefficient(zero(T), zero(T), zero(T), false, zero(T), zero(T), Vector{NTuple{3,T}}(), Vector{NTuple{3,T}}())

# More compact storage for the coefficients to save on size and improve cache locality.
struct LowerTriangularStorage{T}
    data::Vector{T}
end

function LowerTriangularStorage{T}(degree, order) where T
    n = max(degree + 1, order + 1)
    len = (n * (n + 1)) ÷ 2
    return LowerTriangularStorage{T}(zeros(T, len))
end

_ij_to_lt_index(i::Int, j::Int) = (i * (i - 1)) ÷ 2 + j

Base.getindex(L::LowerTriangularStorage, i::Int, j::Int) = L.data[_ij_to_lt_index(i, j)]
Base.setindex!(L::LowerTriangularStorage, v, i::Int, j::Int) = (L.data[_ij_to_lt_index(i, j)] = v)

struct IcgemFile{T<:Number,NT<:Val,Coeff<:AbstractIcgemCoefficient{T}} <: GravityModels.AbstractGravityModel{T,NT}
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
    data::LowerTriangularStorage{Coeff}
end
