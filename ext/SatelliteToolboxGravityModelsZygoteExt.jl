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

function ChainRulesCore.rrule(
    ::typeof(GravityModels.gravitational_potential),
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
) where {T<:Number, V<:Number, NT}

    y = GravityModels.gravitational_potential(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
    )

    function _gravitational_potential_pullback(Δ)
        grad = ForwardDiff.gradient(
            (x) -> GravityModels.gravitational_potential(
                model,
                x[1:3],
                x[4];
                max_degree = max_degree,
                max_order = max_order,
                P = P,
            ),
            [r; time]
        )

        vjp = Δ' * grad

        return (NoTangent(), NoTangent(), vjp[1:3], vjp[4], (NoTangent(), NoTangent(), NoTangent()))

    end

    return y, _gravitational_potential_pullback

end

function ChainRulesCore.rrule(
    ::typeof(GravityModels.gravitational_field_derivative),
    model::AbstractGravityModel{T, NT},
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = -1,
    max_order::Int = -1,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T<:Number, V<:Number, NT}

    y = GravityModels.gravitational_field_derivative(
        model,
        r,
        time;
        max_degree = max_degree,
        max_order = max_order,
        P = P,
        dP = dP,
    )

    function _gravitational_field_derivative_pullback(Δ)
        jac = ForwardDiff.jacobian(
            (x) -> begin
                result = GravityModels.gravitational_field_derivative(
                    model,
                    x[1:3],
                    x[4];
                    max_degree = max_degree,
                    max_order = max_order,
                    P = P,
                    dP = dP,
                )
                # Convert result to a vector for jacobian computation
                return collect(result)
            end,
            vcat(collect(r), time)
        )

        # Handle the tangent properly - extract the values from the tangent structure
        if Δ isa ChainRulesCore.Tangent
            # Extract the backing data from the tangent
            Δ_data = Δ.backing
            if Δ_data isa Tuple
                Δ_vec = map(Δ_data) do x
                    if x isa ChainRulesCore.ZeroTangent
                        0.0  # Convert ZeroTangent to actual zero
                    else
                        ChainRulesCore.unthunk(x)
                    end
                end
                Δ_vec = collect(Δ_vec)
            else
                Δ_vec = collect(ChainRulesCore.unthunk(Δ_data))
            end
        else
            Δ_vec = collect(ChainRulesCore.unthunk(Δ))
        end

        vjp = Δ_vec' * jac

        return (NoTangent(), NoTangent(), vjp[1:3], vjp[4], (NoTangent(), NoTangent(), NoTangent(), NoTangent()))

    end

    return y, _gravitational_field_derivative_pullback

end

end