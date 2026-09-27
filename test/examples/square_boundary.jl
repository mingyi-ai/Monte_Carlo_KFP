@testset "kinetic square boundary data" begin
    left_peak = prescribed_boundary_value(-1.0f0, -0.5f0)
    right_peak = prescribed_boundary_value(1.0f0, 0.5f0)
    @test left_peak ≈ exp(-1.0f0)
    @test right_peak ≈ exp(-1.0f0)

    # Diffusive velocity boundaries are regular everywhere.
    @test prescribed_boundary_value(0.25f0, -1.0f0) == 0.0f0
    @test prescribed_boundary_value(-0.25f0, 1.0f0) == 0.0f0
    @test isnothing(prescribed_boundary_value(0.0f0, 0.9999f0))
    @test isnothing(prescribed_boundary_value(0.9999f0, 0.5f0))

    # Transport-influx and grazing boundaries are not assigned Dirichlet data.
    @test isnothing(prescribed_boundary_value(-1.0f0, 0.5f0))
    @test isnothing(prescribed_boundary_value(-1.0f0, 0.0f0))
    @test isnothing(prescribed_boundary_value(1.0f0, -0.5f0))
    @test isnothing(prescribed_boundary_value(1.0f0, 0.0f0))

    square = SquareDomain(1.0f0)
    left_start = simulation_start(-1.0f0, 0.5f0)
    left_grazing_start = simulation_start(-1.0f0, 0.0f0)
    right_start = simulation_start(1.0f0, -0.5f0)
    @test point_location(square, left_start) == Interior
    @test point_location(square, left_grazing_start) == Interior
    @test point_location(square, right_start) == Interior
    @test simulation_start(0.25f0, -0.5f0) == (-0.5f0, 0.25f0)

    # Wrong-sign side exits caused by a finite Euler step are projected to the
    # zero grazing payoff, not used to define influx boundary grid values.
    @test boundary_payoff(-1.0f0, 0.05f0) == 0.0f0
    @test boundary_payoff(1.0f0, -0.05f0) == 0.0f0
    @test_throws ErrorException boundary_payoff(0.0f0, 0.0f0)
end
