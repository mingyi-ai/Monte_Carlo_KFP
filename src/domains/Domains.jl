"""
Domain geometry and the `AbstractDomain` protocol.

Every concrete shape owns its interior classification, boundary components,
sampling, and segment-crossing behavior. 
"""
module Domains

using LinearAlgebra: dot, norm
using ..MonteCarloKFP: Point

include("protocol.jl")
include("segments.jl")
include("sphere.jl")
include("shapes.jl")

export AbstractDomain,
    SphereSide,
    InsideSphere,
    OutsideSphere,
    PointLocation,
    Interior,
    OnBoundary,
    Exterior,
    BoundaryHit,
    point_location,
    boundary_components,
    boundary_component,
    boundary_measure,
    boundary_point,
    boundary_points,
    boundary_level,
    domain_dimension,
    domain_components,
    domain_scale,
    domain_level,
    component_level,
    domain_hit,
    component_measure,
    component_point,
    first_boundary_hit,
    SquareDomain,
    BoxDomain,
    LShapeDomain,
    BallDomain,
    ShellDomain,
    AnnulusDomain,
    TargetDisks

end
