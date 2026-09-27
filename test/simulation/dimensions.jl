@testset "three-dimensional simulation" begin
    config = SimulationConfig(
        trajectories = 64,
        dt = 1.0f-3,
        max_steps = 4000,
        seed = 5,
    )
    dynamics = OverdampedLangevin(1.0f0)

    ball = BallDomain((0.0f0, 0.0f0, 0.0f0), 1.0f0)
    ball_result = simulate(dynamics, ball, (0.0f0, 0.0f0, 0.0f0), config)
    @test ball_result isa SimulationResult{3,Float32}
    @test all(
        index -> index == 0 || index == 1,
        ball_result.exit_component_indices,
    )
    @test exit_rate(ball_result) > 0.9

    exited = ball_result.exit_component_indices .!= 0
    exit_radii = [
        sqrt(sum(abs2, ball_result.exit_points[index])) for
        index in eachindex(ball_result.exit_points) if exited[index]
    ]
    @test all(isapprox.(exit_radii, 1.0f0; atol = 1.0f-5))

    box = BoxDomain((-1.0f0, -1.0f0, -1.0f0), (1.0f0, 1.0f0, 1.0f0))
    box_result = simulate(dynamics, box, (0.0f0, 0.0f0, 0.0f0), config)
    @test box_result isa SimulationResult{3,Float32}
    @test all(index -> 0 < index <= 6, box_result.exit_component_indices)

    # The kinetic layout is (v, x) with 2n components; a six-dimensional box
    # bounds both the momentum and the position.
    box6 = BoxDomain(ntuple(_ -> -1.0f0, 6), ntuple(_ -> 1.0f0, 6))
    kinetic_result = simulate(
        KineticBrownian(noise = 1.0f0),
        box6,
        ntuple(_ -> 0.0f0, 6),
        config,
    )
    @test kinetic_result isa SimulationResult{6,Float32}
    @test all(index -> 0 < index <= 6, kinetic_result.exit_component_indices)
end

@testset "one-dimensional custom-domain simulation" begin
    interval = IntervalDomain(1.0f0)
    config = SimulationConfig(
        trajectories = 64,
        dt = 1.0f-3,
        max_steps = 2000,
        seed = 7,
    )
    result = simulate(OverdampedLangevin(1.0f0), interval, (0.0f0,), config)

    @test result isa SimulationResult{1,Float32}
    @test all(index -> index <= 2, result.exit_component_indices)
    @test exit_rate(result) > 0.0
end
