@testset "simulation configuration" begin
    @test_throws ArgumentError SimulationConfig(
        trajectories = 0,
        dt = 0.1f0,
        max_steps = 1,
    )
    @test_throws ArgumentError SimulationConfig(
        trajectories = 1,
        dt = 0.0f0,
        max_steps = 1,
    )
end

@testset "deterministic square exit and interpolation" begin
    config =
        SimulationConfig(trajectories = 4, dt = 1.0f0, max_steps = 4, seed = 7)
    # (v, x) start: momentum 0.5 carries the position to the top wall at x = 1.
    result = simulate(
        KineticBrownian(noise = 0.0f0),
        SquareDomain(1.0f0),
        (0.5f0, 0.0f0),
        config,
    )

    @test result.exit_component_indices == fill(3, 4)
    @test result.step_counts == fill(2, 4)
    @test all(isapprox.(getindex.(result.exit_points, 2), 1.0f0))
    @test all(isapprox.(getindex.(result.exit_points, 1), 0.5f0))
    @test exit_rate(result) == 1.0
    @test exit_probability(result, 3) == 1.0
end

@testset "annulus exit components" begin
    config =
        SimulationConfig(trajectories = 3, dt = 1.0f0, max_steps = 4, seed = 11)
    result = simulate(
        KineticBrownian(noise = 0.0f0),
        AnnulusDomain(1.0f0, 3.0f0),
        (-0.5f0, 1.5f0),
        config,
    )

    @test result.exit_component_indices == fill(1, 3)
    radii = hypot.(
        getindex.(result.exit_points, 1),
        getindex.(result.exit_points, 2),
    )
    @test all(isapprox.(radii, 1.0f0; atol = 2.0f-6))
    @test conditional_exit_probability(result, 1) == 1.0
end

@testset "underdamped two-target exit" begin
    config =
        SimulationConfig(trajectories = 2, dt = 1.0f0, max_steps = 3, seed = 19)
    dynamics =
        UnderdampedLangevin(_ -> (0.0f0,); friction = 0.0f0, noise = 0.0f0)
    domain = TargetDisks(-2.0f0, 0.4f0, 1.0f0, 0.6f0)
    result = simulate(dynamics, domain, (0.5f0, 0.0f0), config)

    @test result.exit_component_indices == fill(2, 2)
    radii = hypot.(
        getindex.(result.exit_points, 1),
        getindex.(result.exit_points, 2) .- 1.0f0,
    )
    @test all(isapprox.(radii, 0.6f0; atol = 2.0f-6))
end

@testset "simulation result storage" begin
    config =
        SimulationConfig(trajectories = 3, dt = 0.1f0, max_steps = 2, seed = 3)
    result = simulate(
        KineticBrownian(noise = 1.0f0),
        SquareDomain(5.0f0),
        (0.0f0, 0.0f0),
        config,
    )

    @test result isa SimulationResult{2,Float32}
    @test result.final_points isa Vector{Point{2,Float32}}
    @test result.exit_points isa Vector{Point{2,Float32}}
    @test result.exit_component_indices isa Vector{Int}
    @test result.step_counts isa Vector{Int}
    @test length(result) == 3

    # The domain dimension and the start dimension must agree; start must be a
    # tuple so the state dimension is known statically.
    @test_throws ArgumentError simulate(
        KineticBrownian(noise = 1.0f0),
        SquareDomain(5.0f0),
        (0.0f0, 0.0f0, 0.0f0, 0.0f0),
        config,
    )
    @test_throws ArgumentError simulate(
        KineticBrownian(noise = 1.0f0),
        SquareDomain(5.0f0),
        [0.0f0, 0.0f0],
        config,
    )
end

@testset "optional progress reporting" begin
    config =
        SimulationConfig(trajectories = 4, dt = 1.0f0, max_steps = 4, seed = 7)
    dynamics = KineticBrownian(noise = 0.0f0)
    domain = SquareDomain(1.0f0)

    default_result = simulate(dynamics, domain, (0.5f0, 0.0f0), config)
    enabled_result = redirect_stderr(devnull) do
        simulate(dynamics, domain, (0.5f0, 0.0f0), config; progress = true)
    end
    @test enabled_result.final_points == default_result.final_points

    output = IOBuffer()
    meter = Progress(config.trajectories; dt = 0, output = output)
    custom_result =
        simulate(dynamics, domain, (0.5f0, 0.0f0), config; progress = meter)
    @test custom_result.final_points == default_result.final_points
    @test meter.counter == config.trajectories
    @test occursin("100%", String(take!(output)))
end

@testset "censoring and reproducibility" begin
    config = SimulationConfig(
        trajectories = 32,
        dt = 0.01f0,
        max_steps = 10,
        seed = 1234,
    )
    dynamics = KineticBrownian(noise = 1.0f0)
    domain = SquareDomain(10.0f0)
    serial_result =
        simulate(dynamics, domain, (0.0f0, 0.0f0), config; threaded = false)
    threaded_result =
        simulate(dynamics, domain, (0.0f0, 0.0f0), config; threaded = true)

    @test serial_result.final_points == threaded_result.final_points
    @test all(iszero, serial_result.exit_component_indices)
    @test all(point -> all(isnan, point), serial_result.exit_points)
    @test exit_rate(serial_result) == 0.0
    @test isnan(conditional_exit_probability(serial_result, 1))
end
