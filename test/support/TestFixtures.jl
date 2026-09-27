module TestFixtures

using MonteCarloKFP

export FreeFall, IntervalDomain, OverdampedLangevin

# A minimal custom dynamics: subtype AbstractDynamics, prescribe the noise shape,
# add one `advance` method. FreeFall is deterministic, so it consumes no normals.
struct FreeFall{T<:AbstractFloat} <: AbstractDynamics
    gravity::T
end

MonteCarloKFP.noise_dimension(::FreeFall, ::Point) = 0

MonteCarloKFP.advance(
    dynamics::FreeFall{T},
    state::Point{2,T},
    dt::T,
    ::NTuple{0,T},
) where {T} = Point(state[1] - dynamics.gravity * dt, state[2] + state[1] * dt)

# A d-dimensional position-space dynamics: the state is R^d and every component
# is diffusive. Used to exercise generic-dimension domains and results.
struct OverdampedLangevin{T<:AbstractFloat} <: AbstractDynamics
    noise::T
end

MonteCarloKFP.noise_dimension(::OverdampedLangevin, ::Point{N,T}) where {N,T} =
    N

MonteCarloKFP.advance(
    dynamics::OverdampedLangevin{T},
    state::Point{N,T},
    dt::T,
    normal_increment::NTuple{N,T},
) where {N,T} = Point(
    ntuple(
        index ->
            state[index] + dynamics.noise * sqrt(dt) * normal_increment[index],
        Val(N),
    ),
)

# A custom AbstractDomain with no patches: the one-dimensional interval (-h, h),
# implementing the domain primitives directly. Shows the protocol extension point.
struct IntervalDomain{T<:AbstractFloat} <: AbstractDomain
    halfwidth::T
end

MonteCarloKFP.domain_dimension(::IntervalDomain) = 1
MonteCarloKFP.domain_components(::IntervalDomain) = (:low, :high)
MonteCarloKFP.domain_scale(domain::IntervalDomain) = domain.halfwidth

MonteCarloKFP.component_level(
    domain::IntervalDomain,
    point::Point{1},
    index::Int,
) = index == 1 ? -point[1] - domain.halfwidth : point[1] - domain.halfwidth

MonteCarloKFP.domain_level(domain::IntervalDomain, point::Point{1}) = max(
    MonteCarloKFP.component_level(domain, point, 1),
    MonteCarloKFP.component_level(domain, point, 2),
)

MonteCarloKFP.component_measure(::IntervalDomain, ::Int) = 1.0

MonteCarloKFP.component_point(domain::IntervalDomain, index::Int, ::Real) =
    Point(index == 1 ? -domain.halfwidth : domain.halfwidth)

function MonteCarloKFP.domain_hit(
    domain::IntervalDomain,
    from::Point{1},
    to::Point{1},
)
    delta = to[1] - from[1]
    iszero(delta) && return nothing

    earliest_hit = nothing
    for (component_index, boundary_value) in
        ((1, -domain.halfwidth), (2, domain.halfwidth))
        segment_fraction = (boundary_value - from[1]) / delta
        zero(segment_fraction) <= segment_fraction <= one(segment_fraction) ||
            continue
        if isnothing(earliest_hit) ||
           segment_fraction < earliest_hit.segment_fraction
            earliest_hit = MonteCarloKFP.BoundaryHit(
                Point(boundary_value),
                component_index,
                segment_fraction,
            )
        end
    end
    return earliest_hit
end

end
