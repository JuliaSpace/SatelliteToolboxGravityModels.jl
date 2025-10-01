## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile

abstract type AbstractIcgemCoefficient{T<:Number} end

############################################################################################
#                                 Lower Triangular Storage                                 #
############################################################################################

# We create a custom storage type for the coefficients that only stores the lower triangular
# part of the matrix. This is because the spherical harmonic coefficients are defined only
# for the lower triangular part (i.e., where degree >= order). This saves memory and
# improves cache locality when accessing the coefficients.
struct LowerTriangularStorage{T}
    n::Int
    data::Vector{T}

    function LowerTriangularStorage{T}(num_rows::Int, num_columns::Int) where T
        n   = max(num_rows, num_columns)
        len = (n * (n + 1)) ÷ 2

        return new{T}(n, zeros(T, len))
    end
end

# Auxiliary function to convert (i, j) indices to the index in the 1D array.
_ij_to_lt_index(i::Int, j::Int) = (i * (i - 1)) ÷ 2 + j

# == Julia API =============================================================================

@propagate_inbounds function Base.getindex(
    L::LowerTriangularStorage{T},
    i::Int,
    j::Int
) where T<:AbstractIcgemCoefficient
    @boundscheck (i > L.n || j > L.n || i < 1 || j < 1) && throw(BoundsError(L, [i, j]))

    # For the upper triangular part, return zero.
    j > i && return zero(T)

    return L.data[_ij_to_lt_index(i, j)]
end

@propagate_inbounds function Base.setindex!(L::LowerTriangularStorage, v, i::Int, j::Int)
    # Notice that we also throw an error if trying to set an upper triangular element.
    @boundscheck (i > L.n || j > L.n || i < 1 || j < 1 || j > i) && throw(BoundsError(L, [i, j]))

    L.data[_ij_to_lt_index(i, j)] = v

    return v
end

function Base.summary(io::IO, L::LowerTriangularStorage{T}) where T<:AbstractIcgemCoefficient
    print(io, "$(L.n)×$(L.n) $(typeof(L))")
    return nothing
end

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

    # == Trend =============================================================================

    has_trend::Bool
    trend_clm::T
    trend_slm::T

    # == asin ==============================================================================

    asin_coefficients::Vector{NTuple{3, T}}

    # == acos ==============================================================================

    acos_coefficients::Vector{NTuple{3, T}}
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

function Base.zero(::Type{IcgemGfcCoefficient{T}}) where T
    return IcgemGfcCoefficient(zero(T), zero(T))
end

function Base.zero(::Type{IcgemGfctCoefficient{T}}) where T
    return IcgemGfctCoefficient(
        zero(T),
        zero(T),
        zero(T),
        false,
        zero(T), zero(T),
        Vector{NTuple{3,T}}(),
        Vector{NTuple{3,T}}()
    )
end

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
