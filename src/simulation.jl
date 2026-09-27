using ProgressMeter: Progress, next!

"""Configuration shared by all dynamics and domains."""
struct SimulationConfig{T<:AbstractFloat}
    trajectories::Int
    dt::T
    max_steps::Int
    seed::UInt64
    function SimulationConfig(
        trajectories::Int,
        dt::T,
        max_steps::Int,
        seed::UInt64,
    ) where {T<:AbstractFloat}
        trajectories > 0 ||
            throw(ArgumentError("trajectories must be positive"))
        dt > zero(T) || throw(ArgumentError("dt must be positive"))
        max_steps > 0 || throw(ArgumentError("max_steps must be positive"))
        new{T}(trajectories, dt, max_steps, seed)
    end
end

function SimulationConfig(;
    trajectories::Integer,
    dt::AbstractFloat,
    max_steps::Integer,
    seed::Integer = 0,
)
    seed >= 0 || throw(ArgumentError("seed must be non-negative"))
    return SimulationConfig(Int(trajectories), dt, Int(max_steps), UInt64(seed))
end

"""
Per-trajectory results for an `N`-dimensional simulation.

For trajectory `i`:

- `final_points[i]` is the last Euler–Maruyama integration point;
- `exit_points[i]` is the first interpolated boundary intersection;
- `exit_component_indices[i]` is the one-based index into
  `boundary_components(domain)`, or `0` if the trajectory was censored;
- `step_counts[i]` is the number of integration steps performed.

A censored trajectory did not exit before the finite time horizon. Its exit point
has `NaN` components and its exit component index is zero.
"""
struct SimulationResult{N,T<:AbstractFloat}
    final_points::Vector{Point{N,T}}
    exit_points::Vector{Point{N,T}}
    exit_component_indices::Vector{Int}
    step_counts::Vector{Int}
end

Base.length(result::SimulationResult)::Int =
    length(result.exit_component_indices)

"""Simulate one trajectory and write its outcome into slot `trajectory_index`."""
function _simulate_trajectory!(
    result::SimulationResult{N,T},
    trajectory_index::Int,
    dynamics::AbstractDynamics,
    domain::AbstractDomain,
    initial_point::Point{N,T},
    config::SimulationConfig{T},
)::Nothing where {N,T}
    random_stream = path_rng(config.seed, trajectory_index)
    current_point = initial_point

    for step_index = 1:config.max_steps
        noise_channel_count = noise_dimension(dynamics, current_point)
        standard_normal_increment =
            ntuple(_ -> normal!(random_stream, T), noise_channel_count)
        next_point = advance(
            dynamics,
            current_point,
            config.dt,
            standard_normal_increment,
        )
        boundary_hit = first_boundary_hit(domain, current_point, next_point)

        if !isnothing(boundary_hit)
            result.final_points[trajectory_index] = next_point
            result.exit_points[trajectory_index] = boundary_hit.point
            result.exit_component_indices[trajectory_index] =
                boundary_hit.component_index
            result.step_counts[trajectory_index] = step_index
            return nothing
        end
        current_point = next_point
    end

    result.final_points[trajectory_index] = current_point
    return nothing
end

"""
    simulate(dynamics, domain, start, config; threaded=true, progress=false) -> SimulationResult{N,T}

Simulate independent trajectories using Euler–Maruyama integration. `start` is
an `N`-tuple in the domain's phase space, so the result state dimension is part
of the call signature. Random streams are derived from
`(config.seed, trajectory_index)`, making a run reproducible and independent of
thread scheduling. CPU execution is the stable baseline backend.

Set `progress=true` to display a thread-safe progress bar, or supply a
preconfigured `ProgressMeter.Progress` to control its appearance and output.
Progress reporting is disabled by default so callers can manage progress around
multiple simulations without nested bars.
"""
function simulate(
    dynamics::AbstractDynamics,
    domain::AbstractDomain,
    start::NTuple{N,<:Real},
    config::SimulationConfig{T};
    threaded::Bool = true,
    progress::Union{Bool,Progress} = false,
)::SimulationResult{N,T} where {N,T}
    expected_dimension = domain_dimension(domain)
    N == expected_dimension || throw(
        ArgumentError(
            "start has $N components but the domain is $(expected_dimension)-dimensional",
        ),
    )
    initial_point = Point(map(value -> T(value), start))
    point_location(domain, initial_point) == Interior ||
        throw(ArgumentError("start must be in the interior of the domain"))

    trajectory_count = config.trajectories
    result = SimulationResult(
        Vector{Point{N,T}}(undef, trajectory_count),
        fill(Point(ntuple(_ -> T(NaN), Val(N))), trajectory_count),
        zeros(Int, trajectory_count),
        fill(config.max_steps, trajectory_count),
    )
    progress_meter = if progress isa Progress
        progress
    elseif progress
        Progress(trajectory_count; desc = "Simulating trajectories:")
    else
        nothing
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads :dynamic for trajectory_index = 1:trajectory_count
            _simulate_trajectory!(
                result,
                trajectory_index,
                dynamics,
                domain,
                initial_point,
                config,
            )
            isnothing(progress_meter) || next!(progress_meter)
        end
    else
        for trajectory_index = 1:trajectory_count
            _simulate_trajectory!(
                result,
                trajectory_index,
                dynamics,
                domain,
                initial_point,
                config,
            )
            isnothing(progress_meter) || next!(progress_meter)
        end
    end

    return result
end

"""A vector start is rejected: the state dimension belongs to the type."""
function simulate(
    ::AbstractDynamics,
    ::AbstractDomain,
    start::AbstractVector,
    ::SimulationConfig,
)
    throw(
        ArgumentError(
            "start must be a tuple, not $(typeof(start)); the state dimension is part of the type",
        ),
    )
end

"""Fraction of trajectories that reached any absorbing boundary before the horizon."""
function exit_rate(result::SimulationResult)::Float64
    return count(!iszero, result.exit_component_indices) / length(result)
end

"""
Finite-horizon probability of reaching `component_index`. Censored trajectories
count as failures, so this is a lower bound on an eventual hitting probability.
"""
function exit_probability(
    result::SimulationResult,
    component_index::Integer,
)::Float64
    component_index > 0 ||
        throw(ArgumentError("component_index must be positive"))
    return count(==(component_index), result.exit_component_indices) /
           length(result)
end

"""Probability of `component_index`, conditional on having exited before the horizon."""
function conditional_exit_probability(
    result::SimulationResult,
    component_index::Integer,
)::Float64
    component_index > 0 ||
        throw(ArgumentError("component_index must be positive"))
    exited_count = count(!iszero, result.exit_component_indices)
    iszero(exited_count) && return NaN
    return count(==(component_index), result.exit_component_indices) /
           exited_count
end
