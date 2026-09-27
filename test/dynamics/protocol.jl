@testset "generic point and noise shape" begin
    # Point is a dimension-generic immutable vector.
    @test Point(1.0, 2.0) isa Point{2,Float64}
    @test Point(1.0, 2.0, 3.0) isa Point{3,Float64}
    @test Point([1, 2, 3]) isa Point{3,Float64}
    @test_throws ArgumentError Point()

    four = Point(1.0, 2.0, 3.0, 4.0)
    @test length(four) == 4
    @test four[3] == 3.0
    @test collect(four) == [1.0, 2.0, 3.0, 4.0]

    # Ordinary vector arithmetic preserves the fixed-size Point representation.
    left = Point(1.0f0, 2.0f0)
    right = Point(3.0f0, 5.0f0)
    @test left + right === Point(4.0f0, 7.0f0)
    @test right - left === Point(2.0f0, 3.0f0)
    @test 2 * left === Point(2.0f0, 4.0f0)
    @test left * 2 === Point(2.0f0, 4.0f0)
    @test right / 2 === Point(1.5f0, 2.5f0)
    @test 0.5 * left isa Point{2,Float64}

    # The kinetic interpretation is a validated accessor, not part of Point.
    # (v, x): the leading half is momentum, the trailing half is position.
    @test phase_momentum(four) == Point(1.0, 2.0)
    @test phase_position(four) == Point(3.0, 4.0)
    @test_throws ArgumentError phase_momentum(Point(1.0, 2.0, 3.0))
    @test_throws ArgumentError phase_position(Point(1.0))

    # The dynamics prescribes the noise shape.
    kinetic = KineticBrownian(noise = 0.5)
    @test noise_dimension(kinetic, Point(0.0f0, 0.0f0)) == 1
    @test noise_dimension(kinetic, Point(0.0f0, 0.0f0, 0.0f0, 0.0f0)) == 2
    @test_throws ArgumentError noise_dimension(
        kinetic,
        Point(0.0f0, 0.0f0, 0.0f0),
    )

    # A four-dimensional kinetic step: momentum noise, position drift.
    step = advance(kinetic, four, 0.1, (0.0, 0.0))
    @test step == Point(1.0, 2.0, 3.0 + 1.0 * 0.1, 4.0 + 2.0 * 0.1)
    noisy = advance(kinetic, Point(0.0, 0.0, 1.0, -1.0), 0.25, (0.3, -0.2))
    @test noisy[3] == 1.0 && noisy[4] == -1.0
    @test noisy[1] ≈ 0.5 * sqrt(0.25) * 0.3
    @test noisy[2] ≈ 0.5 * sqrt(0.25) * -0.2
    @test_throws ArgumentError advance(kinetic, four, 0.1, (0.0,))

    # The scalar protocol composes with call-site batching.
    states = [Point(randn(), randn(), randn(), randn()) for _ = 1:64]
    increments = [ntuple(_ -> randn(), 2) for _ = 1:64]
    batched = advance.(Ref(kinetic), states, 0.01, increments)
    @test length(batched) == 64
    @test eltype(batched) == Point{4,Float64}

    threaded = similar(states)
    Threads.@threads for index in eachindex(states)
        threaded[index] =
            advance(kinetic, states[index], 0.01, increments[index])
    end
    @test batched == threaded
end

@testset "custom dynamics protocol" begin
    config =
        SimulationConfig(trajectories = 3, dt = 1.0f0, max_steps = 2, seed = 5)
    result =
        simulate(FreeFall(2.0f0), SquareDomain(10.0f0), (0.0f0, 0.0f0), config)

    @test all(iszero, result.exit_component_indices)
    @test all(
        point -> point[1] == -4.0f0 && point[2] == -2.0f0,
        result.final_points,
    )
end
