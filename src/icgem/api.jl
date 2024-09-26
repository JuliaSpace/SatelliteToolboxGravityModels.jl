## Description #############################################################################
#
# Functions related to the gravity model API.
#
############################################################################################
function GravityModels.coefficients(model::IcgemFile, degree::Int, order::Int, time::DateTime)
    time_JD = (datetime2julian(time) - JD_J2000) * 86400
    return icgem_coefficients(model, degree, order, time_JD)
end


function GravityModels.coefficients(model::IcgemFile, degree::Int, order::Int, time::Number)
    return icgem_coefficients(model, degree, order, time)
end

function GravityModels.coefficient_norm(model::IcgemFile)
    if model.norm === :unnormalized
        return :unnormalized
    else
        return :full
    end
end

GravityModels.gravity_constant(model::IcgemFile) = model.gravity_constant

function GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::DataType = Float64)
    return parse_icgem(filename, T)
end

GravityModels.maximum_degree(model::IcgemFile) = model.max_degree

GravityModels.radius(model::IcgemFile) = model.radius
