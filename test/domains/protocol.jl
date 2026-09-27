@testset "boundary protocol" begin
    square = SquareDomain(1.0f0)
    @test square isa SquareDomain{Float32}
    @test boundary_components(square) == (:bottom, :right, :top, :left)
    @test point_location(square, (0.0f0, 0.0f0)) == Interior
    @test square((0.0f0, 0.0f0)) == Interior
    @test point_location(square, (1.0f0, 0.25f0)) == OnBoundary
    @test point_location(square, (1.01f0, 0.0f0)) == Exterior
    @test boundary_component(square, (1.0f0, 0.25f0)) == 2
    @test boundary_component(square, (0.0f0, 0.0f0)) == 0
    @test boundary_measure(square, :top) == 2.0f0

    square_points = boundary_points(square, 8)
    @test length(square_points) == 8
    @test all(
        point_location(square, point) == OnBoundary for point in square_points
    )
    closed_points = [square_points; [square_points[1]]]
    distances = [
        hypot(
            closed_points[index+1][1] - closed_points[index][1],
            closed_points[index+1][2] - closed_points[index][2],
        ) for index = 1:8
    ]
    @test all(isapprox.(distances, 1.0f0; atol = 2.0f-6))

    annulus = AnnulusDomain(1.0f0, 2.0f0)
    annulus_points = boundary_points(annulus, 12)
    radii = hypot.(
        [point[1] for point in annulus_points],
        [point[2] for point in annulus_points],
    )
    @test count(isapprox(1.0f0; atol = 2.0f-6), radii) == 4
    @test count(isapprox(2.0f0; atol = 2.0f-6), radii) == 8
    @test length(boundary_points(annulus, 7; component = :inner)) == 7

    targets = TargetDisks(-1.0f0, 0.25f0, 1.0f0, 0.25f0)
    @test point_location(targets, (0.0f0, 0.0f0)) == Interior
    @test point_location(targets, (0.0f0, 1.0f0)) == Exterior
    crossing = first_boundary_hit(targets, (0.0f0, 0.0f0), (0.0f0, 2.0f0))
    @test !isnothing(crossing)
    @test crossing.component_index == 2
    @test crossing.point == Point(0.0f0, 0.75f0)
    @test point_location(targets, crossing.point) == OnBoundary

    ball = BallDomain((0.0, 0.0), 1.0)
    @test ball isa AbstractDomain
    @test ball((0.0, 0.0)) == Interior
    @test ball((1.0, 0.0)) == OnBoundary
end

@testset "custom domain protocol" begin
    interval = IntervalDomain(1.0f0)
    @test interval isa AbstractDomain
    @test boundary_components(interval) == (:low, :high)
    @test point_location(interval, (0.0f0,)) == Interior
    @test point_location(interval, (1.0f0,)) == OnBoundary
    @test point_location(interval, (2.0f0,)) == Exterior
    @test boundary_component(interval, (-1.0f0,)) == 1
    @test boundary_component(interval, (1.0f0,)) == 2
    @test boundary_measure(interval, :low) == 1.0

    crossing = first_boundary_hit(interval, (0.0f0,), (2.0f0,))
    @test crossing.component_index == 2
    @test crossing.point == Point(1.0f0)
    @test_throws ArgumentError boundary_points(interval, 4)
end
