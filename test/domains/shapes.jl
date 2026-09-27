@testset "domain construction validation" begin
    @test_throws ArgumentError AnnulusDomain(2.0f0, 1.0f0)
    @test_throws ArgumentError TargetDisks(-1.0f0, 2.0f0, 1.0f0, 2.0f0)
    @test_throws ArgumentError SquareDomain(0.0)
    @test_throws ArgumentError BoxDomain((0.0, 0.0), (1.0, 0.0))
    @test_throws ArgumentError LShapeDomain(-1.0)
    @test_throws ArgumentError BallDomain((0.0, 0.0), -1.0)

    @test SquareDomain(1) isa AbstractDomain
    @test AnnulusDomain(1, 2.0f0) isa AbstractDomain
    @test BallDomain((0.0f0, 0.0f0), 1.0f0) isa AbstractDomain
end

@testset "analytic domains" begin
    # Sphere-based domains implement the protocol directly.
    circle = BallDomain((0.0, 0.0), 1.0)
    @test circle isa AbstractDomain
    @test boundary_components(circle) == (:sphere,)
    @test MonteCarloKFP.domain_dimension(circle) == 2
    @test boundary_measure(circle, :sphere) ≈ 2π
    @test point_location(circle, (1.5, 0.0)) == Exterior
    @test point_location(circle, (1.0, 0.0)) == OnBoundary
    @test isapprox(boundary_point(circle, :sphere, 0.25)[1], 0.0; atol = 1e-12)
    @test all(
        point_location(circle, point) == OnBoundary for
        point in boundary_points(circle, 8)
    )

    # Three-dimensional box from planar faces.
    box = BoxDomain((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
    @test MonteCarloKFP.domain_dimension(box) == 3
    @test length(boundary_components(box)) == 6
    @test boundary_measure(box, :upper_1) == 4.0
    @test point_location(box, (0.0, 0.0, 0.0)) == Interior
    @test point_location(box, (1.0, 0.0, 0.0)) == OnBoundary
    @test point_location(box, (2.0, 0.0, 0.0)) == Exterior
    box_hit = first_boundary_hit(box, (0.0, 0.0, 0.0), (2.0, 0.5, 0.0))
    @test box_hit.point == Point(1.0, 0.25, 0.0)
    @test box_hit.component_index == 2
    @test_throws ArgumentError BoxDomain((1.0, 0.0), (0.0, 1.0))

    # Three-dimensional ball.
    ball = BallDomain((0.0, 0.0, 0.0), 1.0)
    @test MonteCarloKFP.domain_dimension(ball) == 3
    @test boundary_measure(ball, :sphere) ≈ 4π
    @test point_location(ball, (0.0, 0.0, 0.0)) == Interior
    @test point_location(ball, (1.0, 0.0, 0.0)) == OnBoundary
    @test point_location(ball, (2.0, 0.0, 0.0)) == Exterior
    ball_hit = first_boundary_hit(ball, (0.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    @test ball_hit.point == Point(1.0, 0.0, 0.0)
    @test ball_hit.segment_fraction == 0.5

    # Three-dimensional shell: outside the inner sphere, inside the outer one.
    shell = ShellDomain((0.0, 0.0, 0.0), 1.0, 2.0)
    @test boundary_components(shell) == (:inner, :outer)
    @test point_location(shell, (1.5, 0.0, 0.0)) == Interior
    @test point_location(shell, (0.5, 0.0, 0.0)) == Exterior
    @test point_location(shell, (2.5, 0.0, 0.0)) == Exterior
    shell_hit = first_boundary_hit(shell, (0.0, 0.0, 1.5), (0.0, 0.0, 3.0))
    @test shell_hit.component_index == 2
    @test shell_hit.point == Point(0.0, 0.0, 2.0)

    # Equal-arclength sampling is a one-dimensional procedure.
    @test_throws ArgumentError boundary_points(box, 8)

    # Concave square with its top-right quadrant removed.
    lshape = LShapeDomain(1.0f0)
    @test boundary_components(lshape) ==
          (:bottom, :right, :notch_bottom, :notch_left, :top, :left)
    @test point_location(lshape, (-0.5f0, 0.5f0)) == Interior
    @test point_location(lshape, (0.5f0, -0.5f0)) == Interior
    @test point_location(lshape, (0.5f0, 0.5f0)) == Exterior
    @test point_location(lshape, (0.0f0, 0.5f0)) == OnBoundary
    @test point_location(lshape, (0.5f0, 1.0f0)) == Exterior
    @test boundary_component(lshape, (0.5f0, 0.0f0)) == 3
    @test boundary_component(lshape, (0.0f0, 0.5f0)) == 4
    @test boundary_component(lshape, (-0.5f0, 0.0f0)) == 0
    @test sum(boundary_measure(lshape, index) for index = 1:6) == 8.0f0
    @test all(
        point_location(lshape, point) == OnBoundary for
        point in boundary_points(lshape, 16)
    )

    notch_hit = first_boundary_hit(lshape, (-0.5f0, 0.5f0), (0.5f0, 0.5f0))
    @test notch_hit.point == Point(0.0f0, 0.5f0)
    @test notch_hit.component_index == 4
    @test notch_hit.segment_fraction == 0.5f0
end
