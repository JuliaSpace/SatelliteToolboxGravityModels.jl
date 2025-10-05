## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxGravityModels;
        ambiguities      = (recursive    = false),
        deps_compat      = (check_extras = false),
        persistent_tasks = false
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        SatelliteToolboxGravityModels;
        toplevel_logger = nothing,
        target_modules  = (@__MODULE__,)
    )
end

@testset "Gravity Model Allocations" begin
    import SatelliteToolboxGravityModels: IcgemGfcCoefficient, IcgemGfctCoefficient
    model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    max_deg = 10
    max_ord = 10
    P = zeros(max_deg, max_ord)
    dP = zeros(max_deg, max_ord)

    @test length(
        check_allocs(
            (model, x) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = max_deg,
                max_order = max_ord,
                P = P,
                dP = dP
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = max_deg,
                max_order = max_ord,
                P = P,
                dP = dP
            ),
            (IcgemFile{Float64, Val{:unnormalized}, IcgemGfcCoefficient{Float64}}, Vector{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (x) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = max_deg,
                max_order = max_ord,
                P = P,
                dP = dP
            ),
            (Vector{Float64},)
        )
    ) == 0

    @test length(
        check_allocs(
            (r_itrf, x) -> GravityModels.gravitational_acceleration(
                model,
                r_itrf,
                x;
                max_degree = max_deg,
                max_order = max_ord,
                P = P,
                dP = dP
            ),
            (Vector{Float64}, Float64)
        )
    ) == 0
end
