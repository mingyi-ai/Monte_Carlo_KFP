using Aqua
using MonteCarloKFP
using ProgressMeter
using Test

include("support/TestFixtures.jl")
using .TestFixtures

@testset "MonteCarloKFP" begin
    @testset "Package quality" begin
        Aqua.test_all(MonteCarloKFP)
    end

    @testset "Domains" begin
        include("domains/protocol.jl")
        include("domains/shapes.jl")
    end

    @testset "Dynamics" begin
        include("dynamics/protocol.jl")
    end

    @testset "Simulation" begin
        include("simulation/core.jl")
        include("simulation/dimensions.jl")
    end

    @testset "Examples" begin
        include("../examples/square_boundary.jl")
        include("examples/square_boundary.jl")
    end
end
