## Description #############################################################################
#
# Functions related to the gravity model API.
#
############################################################################################

function GravityModels.coefficients(model::IcgemFile, degree::Int, order::Int, time::Number)
    return icgem_coefficients(model, degree, order, time)
end

function GravityModels.coefficient_norm(
    model::IcgemFile{T, Val{NT}}
) where {T <: Number, NT}
    if NT === :unnormalized
        return :unnormalized
    else
        return :full
    end
end

GravityModels.gravity_constant(model::IcgemFile) = model.gravity_constant

"""
    GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::Type = Float64) -> IcgemFile

Load the ICGEM file `filename` and return an [`IcgemFile`](@ref) object with its parsed
data. `T` is converted to float to obtain the type of the model coefficients. The function
throws an [`IcgemParseError`](@ref) if the file does not conform to the ICGEM format.

See also: [`parse_icgem`](@ref), [`fetch_icgem_file`](@ref)
"""
function GravityModels.load(
    ::Type{IcgemFile}, filename::AbstractString, ::Type{T} = Float64
) where {T}
    return parse_icgem(filename, T)
end

GravityModels.maximum_degree(model::IcgemFile) = model.max_degree

GravityModels.radius(model::IcgemFile) = model.radius
