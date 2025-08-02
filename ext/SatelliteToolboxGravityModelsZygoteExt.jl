## Description #############################################################################
#
# Zygote Extension for the SatelliteToolboxGravityModels.jl package. Needed since gravitational_acceleration
# is mutating and changing that would require a large rework. Instead use ForwardDiff for this function instead.
#
############################################################################################

module SatelliteToolboxGravityModelsZygoteExt

using SatelliteToolboxGravityModels

using Zygote: NoTangent
using Zygote.ChainRulesCore: ChainRulesCore

using ForwardDiff

function ChainRulesCore.rrule(
    ::typeof(GravityModels.gravitational_acceleration),
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT}

    y = GravityModels.gravitational_acceleration(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
    )

    function _gravitational_acceleration_pullback(Δ)
        jac = ForwardDiff.jacobian(
            (x) -> GravityModels.gravitational_acceleration(
                model,
                x[1:3],
                x[4];
                max_degree = max_degree,
                max_order = max_order,
                P = P,
                dP = dP,
            ),
            [r; time]
        )

        vjp = Δ' * jac

        return (NoTangent(), NoTangent(), vjp[1:3], vjp[4], (NoTangent(), NoTangent(), NoTangent(), NoTangent()))

    end

    return y, _gravitational_acceleration_pullback

end 

end