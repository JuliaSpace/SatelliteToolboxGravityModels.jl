using Test

using SatelliteToolboxGravityModels

using Dates
using DelimitedFiles
using LinearAlgebra
using SatelliteToolboxTransformations
using Scratch

# We must clear the scratch space before the tests to avoid errors.
clear_scratchspaces!(SatelliteToolboxGravityModels)

@testset "ICGEM" verbose = true begin
    include("./icgem.jl")
end

@testset "GravityModels API" verbose = true begin
    include("./gravity_models.jl")
end


if isempty(VERSION.prerelease)
    # Add Mooncake and Enzyme to the project if not the nightly version
    # Adding them via the Project.toml isn't working because it tries to compile them before reaching the gating
    using Pkg
    Pkg.add("DifferentiationInterface")
    Pkg.add("FiniteDiff")
    Pkg.add("ForwardDiff")
    Pkg.add("Mooncake")
    Pkg.add("PolyesterForwardDiff")
    Pkg.add("Zygote")

    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    #TODO: Remove this as Enzyme adds more support for Julia 1.12
    if (VERSION.major == 1 && VERSION.minor < 12)
        Pkg.add("Enzyme")
        using DifferentiationInterface
        using Enzyme, FiniteDiff, ForwardDiff, Mooncake, PolyesterForwardDiff, Zygote

        const _BACKENDS = (
            ("ForwardDiff", AutoForwardDiff()),
            ("Enzyme", AutoEnzyme()),
            ("Mooncake", AutoMooncake(;config=nothing)),
            ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
            ("Zygote", AutoZygote()),
        )
    else
        @warn "Enzyme is not fully supported on Julia 1.12, skipping tests"
        using DifferentiationInterface
        using FiniteDiff, ForwardDiff, Mooncake, PolyesterForwardDiff, Zygote

        const _BACKENDS = (
            ("ForwardDiff", AutoForwardDiff()),
            ("Mooncake", AutoMooncake(;config=nothing)),
            ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
            ("Zygote", AutoZygote()),
        )
    end

    using JET
    using AllocCheck
    using Aqua

    import SatelliteToolboxGravityModels: IcgemGfcCoefficient, IcgemGfctCoefficient
    import SatelliteToolboxBase: LowerTriangularStorage, RowMajor

    # Skip allocation tests on macOS with Julia 1.12+ due to AllocCheck detecting
    # platform-specific runtime calls (jl_get_pgcstack_static) as allocations
    if Sys.isapple() && (VERSION.major == 1 && VERSION.minor >= 12)
        @warn "Allocation tests skipped on macOS with Julia 1.12+ due to AllocCheck platform limitations"
    else
        @testset "Performance and Memory Allocations" verbose = true begin
            include("./allocations.jl")
        end
    end
    
    @testset "Differentiation Tests" verbose = true begin
        include("differentiability.jl")
    end
else
    @warn "Differentiation backends not guaranteed to work on julia-nightly, skipping tests"
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end
