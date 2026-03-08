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
    using Pkg
    Pkg.add("DifferentiationInterface")
    Pkg.add("ForwardDiff")
    Pkg.add("Zygote")
    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    @testset "Zygote Extension" verbose = true begin
        include("./zygote_ext.jl")
    end

    using JET
    using AllocCheck
    using Aqua

    import SatelliteToolboxGravityModels: IcgemGfcCoefficient, IcgemGfctCoefficient
    import SatelliteToolboxBase: LowerTriangularStorage, RowMajor

    if Sys.isapple() && (VERSION.major == 1 && VERSION.minor >= 12)
        @warn "Allocation tests skipped on macOS with Julia 1.12+ due to AllocCheck platform limitations"
    else
        @testset "Performance Tests" verbose = true begin
            include("./allocations.jl")
        end
    end
else
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end
