## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxGravityModels;
        ambiguities = (recursive = false),
        deps_compat = (check_extras = false),
    )
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET Testing" begin
        rep = JET.test_package(
            SatelliteToolboxGravityModels;
            toplevel_logger = nothing,
            target_modules  = (@__MODULE__,),
        )
    end
end

@testset "Gravity Model Allocations" begin
    for norm in (:full, :unnormalized)
        model_type = IcgemFile{Float64, Val{norm}}

        @test length(
            check_allocs(
                (model, r_itrf, md, mo, w) -> GravityModels.gravitational_acceleration(
                    model, r_itrf; max_degree = md, max_order = mo, workspace = w
                ),
                (model_type, Vector{Float64}, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (model, r_itrf, t, md, mo, w) -> GravityModels.gravitational_acceleration(
                    model, r_itrf, t; max_degree = md, max_order = mo, workspace = w
                ),
                (model_type, Vector{Float64}, Float64, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (model, r_itrf, t, md, mo, w) -> GravityModels.gravity_acceleration(
                    model, r_itrf, t; max_degree = md, max_order = mo, workspace = w
                ),
                (model_type, Vector{Float64}, Float64, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (model, r_itrf, t, md, mo, w) ->
                    GravityModels.gravitational_field_derivative(
                        model, r_itrf, t; max_degree = md, max_order = mo, workspace = w
                    ),
                (model_type, Vector{Float64}, Float64, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0
    end
end

@testset "Gravitational Potential Allocations" begin
    for norm in (:full, :unnormalized)
        model_type = IcgemFile{Float64, Val{norm}}

        @test length(
            check_allocs(
                (model, r_itrf, md, mo, w) -> GravityModels.gravitational_potential(
                    model, r_itrf; max_degree = md, max_order = mo, workspace = w
                ),
                (model_type, Vector{Float64}, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (model, r_itrf, t, md, mo, w) -> GravityModels.gravitational_potential(
                    model, r_itrf, t; max_degree = md, max_order = mo, workspace = w
                ),
                (model_type, Vector{Float64}, Float64, Int, Int, Workspace{norm, Float64}),
            ),
        ) == 0
    end
end
