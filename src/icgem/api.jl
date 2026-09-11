## Description #############################################################################
#
# Functions related to the gravity model API.
#
############################################################################################

GravityModels.angular_speed(model::IcgemFile) = model.angular_speed

function GravityModels.coefficients(model::IcgemFile, degree::Int, order::Int, time::Number)
    return icgem_coefficients(model, degree, order, time)
end

GravityModels.coefficient_norm(model::IcgemFile) = model.norm

GravityModels.gravity_constant(model::IcgemFile) = model.gravity_constant

"""
    GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::Type = Float64; kwargs...) -> IcgemFile

Load the ICGEM file `filename` and return an [`IcgemFile`](@ref) object with its parsed
data. `T` is converted to float to obtain the type of the model coefficients. The function
throws an [`IcgemParseError`](@ref) if the file does not conform to the ICGEM format.

See also: [`parse_icgem`](@ref), [`fetch_icgem_file`](@ref)

# Keywords

- `angular_speed::Number`: Angular speed [rad/s] of the central body, which is not
    defined in the ICGEM file and is used to compute the centrifugal acceleration. It
    must be provided for models of bodies other than Earth.
    (**Default**: `EARTH_ANGULAR_SPEED`)
"""
function GravityModels.load(
    ::Type{IcgemFile},
    filename::AbstractString,
    ::Type{T} = Float64;
    angular_speed::Number = EARTH_ANGULAR_SPEED,
) where {T}
    return parse_icgem(filename, T; angular_speed = angular_speed)
end

GravityModels.maximum_degree(model::IcgemFile) = model.max_degree

GravityModels.radius(model::IcgemFile) = model.radius
