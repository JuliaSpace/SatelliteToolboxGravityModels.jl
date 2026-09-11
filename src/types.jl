## Description #############################################################################
#
# Definition of types and structures.
#
############################################################################################

export IcgemFile, IcgemParseError

############################################################################################
#                                       Exceptions                                        #
############################################################################################

"""
    struct IcgemParseError <: Exception

Exception thrown when an ICGEM file does not conform to the ICGEM format.

# Fields

- `message::String`: Description of the problem found in the file.
- `line::Int`: Number of the line in which the problem was found, or 0 if the problem is
    not related to a specific line.
"""
struct IcgemParseError <: Exception
    message::String
    line::Int
end

"""
    IcgemParseError(message::String) -> IcgemParseError

Create an [`IcgemParseError`](@ref) with the `message` that is not related to a specific
line of the file.
"""
IcgemParseError(message::String) = IcgemParseError(message, 0)

function Base.showerror(io::IO, e::IcgemParseError)
    print(io, "IcgemParseError: ")
    (e.line > 0) && print(io, "[Line ", e.line, "] ")
    print(io, e.message)
    return nothing
end

############################################################################################
#                                          ICGEM                                           #
############################################################################################

"""
    struct IcgemGfcCoefficient{T <: Number}

Store a constant (`gfc`) spherical harmonics coefficient of an ICGEM file.

# Fields

- `clm::T`: Cosine coefficient `Clm` [-].
- `slm::T`: Sine coefficient `Slm` [-].
"""
struct IcgemGfcCoefficient{T <: Number}
    clm::T
    slm::T
end

# The coefficients are scalars in broadcasting operations, allowing, for example, to fill a
# storage with `zero(IcgemGfcCoefficient{T})`.
Base.broadcastable(c::IcgemGfcCoefficient) = Ref(c)

function Base.zero(::Type{IcgemGfcCoefficient{T}}) where {T}
    return IcgemGfcCoefficient(zero(T), zero(T))
end

"""
    struct IcgemPeriodicTerm{T <: Number}

Store a periodic term of a time-variable spherical harmonics coefficient of an ICGEM file.
The sine (`asin`) and cosine (`acos`) terms with the same period are stored together.

# Fields

- `amplitude_sin_clm::T`: Amplitude [-] of the sine term of `Clm`.
- `amplitude_sin_slm::T`: Amplitude [-] of the sine term of `Slm`.
- `amplitude_cos_clm::T`: Amplitude [-] of the cosine term of `Clm`.
- `amplitude_cos_slm::T`: Amplitude [-] of the cosine term of `Slm`.
- `period::T`: Period [year] of the term.
"""
struct IcgemPeriodicTerm{T <: Number}
    amplitude_sin_clm::T
    amplitude_sin_slm::T
    amplitude_cos_clm::T
    amplitude_cos_slm::T
    period::T
end

"""
    struct IcgemTimeVariableCoefficient{T <: Number}

Store a time-variable (`gfct`) spherical harmonics coefficient of an ICGEM file, together
with its linear trend (`trnd`) and periodic terms (`asin` and `acos`).

The coefficient is valid in the interval `[t₀, t₁)`. Files in the ICGEM format 1.0 do not
define the end of the validity interval, in which case `t₁` is `Inf`. Files in the ICGEM
format 2.0 can define several validity intervals for the same degree and order, each one
stored in a different object.

# Fields

- `degree::Int`: Degree of the coefficient.
- `order::Int`: Order of the coefficient.
- `clm::T`: Cosine coefficient `Clm` [-] at the epoch `t₀`.
- `slm::T`: Sine coefficient `Slm` [-] at the epoch `t₀`.
- `t₀::T`: Epoch of the coefficient, expressed as the number of elapsed seconds [s] since
    the J2000.0 epoch (2000-01-01T12:00:00).
- `t₁::T`: End of the validity interval, expressed as the number of elapsed seconds [s]
    since the J2000.0 epoch (2000-01-01T12:00:00), or `Inf` if the file does not define it.
- `trend_clm::T`: Linear trend of `Clm` [year⁻¹].
- `trend_slm::T`: Linear trend of `Slm` [year⁻¹].
- `periodic_terms::Vector{IcgemPeriodicTerm{T}}`: Periodic terms of the coefficient.
"""
struct IcgemTimeVariableCoefficient{T <: Number}
    degree::Int
    order::Int
    clm::T
    slm::T
    t₀::T
    t₁::T
    trend_clm::T
    trend_slm::T
    periodic_terms::Vector{IcgemPeriodicTerm{T}}
end

"""
    struct IcgemFile{T <: Number, N <: Val} <: GravityModels.AbstractGravityModel{T}

Store the information of a parsed ICGEM file.

The coefficients are stored in the field `data`. The time-variable coefficients are stored
in the field `time_variable_coefficients`, sorted by degree, order, and epoch, and the
field `time_variable_index` maps the degree and order to the index of the first object in
that vector related to them, or 0 if the coefficient is constant. This index is defined
only up to the degree `max_time_variable_degree`.

# Fields

- `product_type::Symbol`: Product type of the model.
- `model_name::String`: Name of the gravity model.
- `gravity_constant::T`: Gravity constant [m³/s²] of the central body.
- `radius::T`: Reference radius [m] of the model.
- `angular_speed::T`: Angular speed [rad/s] of the central body, used to compute the
    centrifugal acceleration. It is not defined in the ICGEM file and must be provided
    when the model is loaded.
- `max_degree::Int`: Maximum degree available in the model.
- `errors::Symbol`: Type of the errors described in the file (`:no`, `:calibrated`,
    `:calibrated_and_formal`, or `:formal`).
- `tide_system::Symbol`: Tide system of the model, or `:unknown` if the file does not
    specify it.
- `norm::N`: Normalization of the model coefficients wrapped in a `Val`, as returned by
    [`GravityModels.coefficient_norm`](@ref): `Val(:full)` for fully normalized
    coefficients or `Val(:unnormalized)` for unnormalized ones.
- `data::LowerTriangularStorage{RowMajor, IcgemGfcCoefficient{T}}`: Spherical harmonics
    coefficients of the model, in which the element `[n + 1, m + 1]` is the coefficient of
    degree `n` and order `m`. For time-variable coefficients, it contains the values at the
    epoch of the first validity interval.
- `max_time_variable_degree::Int`: Maximum degree of the time-variable coefficients, or -1
    if the model has only constant coefficients.
- `time_variable_index::LowerTriangularStorage{RowMajor, Int32}`: Index in which the
    element `[n + 1, m + 1]` is the position of the first time-variable coefficient of
    degree `n` and order `m` in `time_variable_coefficients`, or 0 if the coefficient is
    constant. It is defined only up to the degree `max_time_variable_degree`.
- `time_variable_coefficients::Vector{IcgemTimeVariableCoefficient{T}}`: Time-variable
    coefficients of the model sorted by degree, order, and epoch.
"""
struct IcgemFile{T <: Number, N <: Val} <: GravityModels.AbstractGravityModel{T}
    # Fields related to the header.
    product_type::Symbol
    model_name::String
    gravity_constant::T
    radius::T
    angular_speed::T
    max_degree::Int
    errors::Symbol
    tide_system::Symbol
    norm::N

    # Fields related to the data section.
    data::LowerTriangularStorage{RowMajor, IcgemGfcCoefficient{T}}
    max_time_variable_degree::Int
    time_variable_index::LowerTriangularStorage{RowMajor, Int32}
    time_variable_coefficients::Vector{IcgemTimeVariableCoefficient{T}}
end
