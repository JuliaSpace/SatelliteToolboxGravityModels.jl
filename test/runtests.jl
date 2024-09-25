using Test

using SatelliteToolboxGravityModels

using Dates
using DelimitedFiles
using LinearAlgebra
using SatelliteToolboxTransformations
using Scratch

using DifferentiationInterface
import FiniteDiff, ForwardDiff

# We must clear the scratch space before the tests to avoid errors.
clear_scratchspaces!(SatelliteToolboxGravityModels)

@testset "ICGEM" verbose = true begin
    include("./icgem.jl")
end

@testset "GravityModels API" verbose = true begin
    include("./gravity_models.jl")
end

#TODO: REMOVE VERSION CHECK WHEN 1.6 STOPS BEING SUPPORTED
if VERSION != v"1.6.0"
    @testset "Differentiation Tests" verbose = true begin
        include("differentiability.jl")
    end
end


