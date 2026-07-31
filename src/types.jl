## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile

"""
    abstract type AbstractIcgemCoefficient{T<:Number}

Abstract type of all spherical harmonics coefficients stored in an ICGEM file.
"""
abstract type AbstractIcgemCoefficient{T <: Number} end

############################################################################################
#                                          ICGEM                                           #
############################################################################################

"""
    struct IcgemGfcCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}

Store a constant (`gfc`) spherical harmonics coefficient of an ICGEM file.

# Fields

- `clm::T`: Cosine coefficient `Clm` [-].
- `slm::T`: Sine coefficient `Slm` [-].
"""
struct IcgemGfcCoefficient{T <: Number} <: AbstractIcgemCoefficient{T}
    clm::T
    slm::T
end

"""
    struct IcgemGfctCoefficient{T<:Number} <: AbstractIcgemCoefficient{T}

Store a time-variable (`gfct`) spherical harmonics coefficient of an ICGEM file.

# Fields

- `clm::T`: Cosine coefficient `Clm` [-] at the epoch `time`.
- `slm::T`: Sine coefficient `Slm` [-] at the epoch `time`.
- `time::T`: Epoch of the coefficients, expressed as the number of elapsed seconds [s]
    since the J2000.0 epoch (2000-01-01T12:00:00).
- `is_time_varying::Bool`: Indicate whether the coefficient is time-varying. If `false`,
    the other time-related fields are ignored and the coefficient is treated as a regular
    [`IcgemGfcCoefficient`](@ref). This is useful to avoid unnecessary computations,
    leading to a huge performance boost.
- `has_trend::Bool`: Indicate whether the coefficient has a linear trend.
- `trend_clm::T`: Linear trend of `Clm` [year⁻¹].
- `trend_slm::T`: Linear trend of `Slm` [year⁻¹].
- `asin_coefficients::Vector{NTuple{3, T}}`: Sine periodic terms, in which each element
    contains the amplitude for `Clm` [-], the amplitude for `Slm` [-], and the period
    [year].
- `acos_coefficients::Vector{NTuple{3, T}}`: Cosine periodic terms, in which each element
    contains the amplitude for `Clm` [-], the amplitude for `Slm` [-], and the period
    [year].
"""
struct IcgemGfctCoefficient{T <: Number} <: AbstractIcgemCoefficient{T}
    clm::T
    slm::T
    time::T # .................................................. Seconds since J2000.0 epoch
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

"""
    IcgemGfctCoefficient(c::IcgemGfcCoefficient{T}) -> IcgemGfctCoefficient{T}

Create an [`IcgemGfctCoefficient`](@ref) from the constant coefficient `c`, keeping `Clm`
and `Slm` and marking the result as not time-varying.
"""
IcgemGfctCoefficient(c::IcgemGfcCoefficient{T}) where {T} = IcgemGfctCoefficient(
    c.clm,
    c.slm,
    zero(T),
    false,
    false,
    zero(T),
    zero(T),
    NTuple{3, T}[],
    NTuple{3, T}[],
)

function Base.zero(::Type{IcgemGfcCoefficient{T}}) where {T}
    return IcgemGfcCoefficient(zero(T), zero(T))
end

function Base.zero(::Type{IcgemGfctCoefficient{T}}) where {T}
    return IcgemGfctCoefficient(
        zero(T),
        zero(T),
        zero(T),
        false,
        false,
        zero(T),
        zero(T),
        Vector{NTuple{3, T}}(),
        Vector{NTuple{3, T}}(),
    )
end

"""
    struct IcgemFile{T<:Number, NT<:Val, Coeff<:AbstractIcgemCoefficient{T}} <: GravityModels.AbstractGravityModel{T, NT}

Store the information of a parsed ICGEM file.

# Fields

- `product_type::Symbol`: Product type of the model.
- `model_name::String`: Name of the gravity model.
- `gravity_constant::T`: Gravity constant [m³/s²] of the central body.
- `radius::T`: Reference radius [m] of the model.
- `max_degree::Int`: Maximum degree available in the model.
- `errors::Symbol`: Type of the errors described in the file (`:no`, `:calibrated`,
    `:calibrated_and_formal`, or `:formal`).
- `tide_system::Symbol`: Tide system of the model, or `:unknown` if the file does not
    specify it.
- `norm::NT`: Normalization of the model coefficients wrapped in a `Val`.
- `data::LowerTriangularStorage{RowMajor, Coeff}`: Spherical harmonics coefficients of the
    model, in which the element `[n + 1, m + 1]` is the coefficient of degree `n` and
    order `m`.
"""
struct IcgemFile{T <: Number, NT <: Val, Coeff <: AbstractIcgemCoefficient{T}} <:
       GravityModels.AbstractGravityModel{T, NT}
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
