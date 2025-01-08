using Test

using SatelliteToolboxGravityModels

using Dates
using DelimitedFiles
using LinearAlgebra
using SatelliteToolboxTransformations
using Scratch

using DifferentiationInterface
using FiniteDiff, ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote

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
