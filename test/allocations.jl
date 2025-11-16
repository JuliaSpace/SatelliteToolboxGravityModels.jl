## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(
        SatelliteToolboxGravityModels;
        ambiguities      = (recursive    = false),
        deps_compat      = (check_extras = false)
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        SatelliteToolboxGravityModels;
        toplevel_logger = nothing,
        target_modules  = (@__MODULE__,)
    )
end

@testset "Gravity Model Allocations (Matrix)" begin
    model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    max_deg = 10
    max_ord = 10
    P = zeros(max_deg, max_ord)
    dP = zeros(max_deg, max_ord)

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64}, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:unnormalized}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64}, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64}, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, r_itrf, t, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                r_itrf,
                t;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Float64, Int, Int, Matrix{Float64}, Matrix{Float64})
        )
    ) == 0
end

@testset "Gravitational Potential Allocations (Matrix)" begin
    model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    max_deg = 10
    max_ord = 10
    P = zeros(max_deg, max_ord)

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:unnormalized}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, Matrix{Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, r_itrf, t, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                r_itrf,
                t;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Float64, Int, Int, Matrix{Float64})
        )
    ) == 0
end

@testset "Gravity Model Allocations (LowerTriangularStorage)" begin
    model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    max_deg = 10
    max_ord = 10
    P = LowerTriangularStorage{RowMajor, Float64}(max_deg + 1)
    dP = LowerTriangularStorage{RowMajor, Float64}(max_deg + 1)

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:unnormalized}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, r_itrf, t, md, mo, p, dp) -> GravityModels.gravitational_acceleration(
                model,
                r_itrf,
                t;
                max_degree = md,
                max_order = mo,
                P = p,
                dP = dp
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Float64, Int, Int, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0
end

@testset "Gravitational Potential Allocations (LowerTriangularStorage)" begin
    model = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    max_deg = 10
    max_ord = 10
    P = LowerTriangularStorage{RowMajor, Float64}(max_deg + 1)

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:unnormalized}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, x, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                x;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Int, Int, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0

    @test length(
        check_allocs(
            (model, r_itrf, t, md, mo, p) -> GravityModels.gravitational_potential(
                model,
                r_itrf,
                t;
                max_degree = md,
                max_order = mo,
                P = p
            ),
            (IcgemFile{Float64, Val{:full}, IcgemGfcCoefficient{Float64}}, Vector{Float64}, Float64, Int, Int, LowerTriangularStorage{RowMajor, Float64})
        )
    ) == 0
end
