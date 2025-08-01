using Test

using SatelliteToolboxGravityModels

using Dates
using DelimitedFiles
using LinearAlgebra
using SatelliteToolboxTransformations
using Scratch

using DifferentiationInterface
using FiniteDiff, ForwardDiff, PolyesterForwardDiff, Zygote

if isempty(VERSION.prerelease)
    # Add Mooncake and Enzyme to the project if not the nightly version
    # Adding them via the Project.toml isn't working because it tries to compile them before reaching the gating
    using Pkg
    Pkg.add("Mooncake", "Enzyme")

    # Test with Mooncake and Enzyme along with the other backends
    using Mooncake, Enzyme
    const _BACKENDS = (
        ("ForwardDiff", AutoForwardDiff()),
        ("Enzyme", AutoEnzyme()),
        ("Mooncake", AutoMooncake(;config=nothing)),
        ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
        ("Zygote", AutoZygote()),
    )
else
    @warn "Mooncake.jl and Enzyme.jl not guaranteed to work on julia-nightly, skipping tests"
    const _BACKENDS = (
        ("ForwardDiff", AutoForwardDiff()),
        ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
        ("Zygote", AutoZygote()),
    )
end

using Aqua
using JET
using AllocCheck

# We must clear the scratch space before the tests to avoid errors.
clear_scratchspaces!(SatelliteToolboxGravityModels)

@testset "ICGEM" verbose = true begin
    include("./icgem.jl")
end

@testset "GravityModels API" verbose = true begin
    include("./gravity_models.jl")
end

@testset "Differentiation Tests" verbose = true begin
    include("differentiability.jl")
end

@testset "Performance and Memory Allocations" verbose = true begin
    include("./allocations.jl")
end
