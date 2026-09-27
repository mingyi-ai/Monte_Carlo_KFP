module MonteCarloKFP

include("random.jl")
include("dynamics.jl")
include("domains/Domains.jl")
using .Domains
include("simulation.jl")

# Simulation workflow
export SimulationConfig,
    SimulationResult,
    simulate,
    exit_rate,
    exit_probability,
    conditional_exit_probability

# Dynamics and the extension protocol
export Point,
    phase_momentum,
    phase_position,
    AbstractDynamics,
    advance,
    noise_dimension,
    KineticBrownian,
    UnderdampedLangevin

# Domains
export AbstractDomain,
    InsideSphere,
    OutsideSphere,
    SquareDomain,
    BoxDomain,
    LShapeDomain,
    BallDomain,
    ShellDomain,
    AnnulusDomain,
    TargetDisks

# Classification and boundary queries
export Interior,
    OnBoundary,
    Exterior,
    BoundaryHit,
    point_location,
    boundary_components,
    boundary_component,
    boundary_measure,
    boundary_point,
    boundary_points,
    first_boundary_hit

# Unexported but public extension points: `SphereSide` and `PointLocation`
# (the enum types behind `InsideSphere`/`OutsideSphere` and
# `Interior`/`OnBoundary`/`Exterior`, which stay exported), the domain protocol
# primitives (`domain_dimension`, `domain_components`, `domain_scale`,
# `domain_level`, `component_level`, `domain_hit`, `component_measure`,
# `component_point`), plus `boundary_level`.

end
