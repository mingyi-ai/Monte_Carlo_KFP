"""
A point in `N`-dimensional phase space, stored as an immutable tuple.

`Point` behaves as an `AbstractVector`: `point[index]` and `length(point)` work,
and generic array code accepts it. `Point(v, x)` is the two-dimensional kinetic
case. The point itself carries no position/momentum interpretation; see
[`phase_position`](@ref) and [`phase_momentum`](@ref).
"""
struct Point{N,T<:AbstractFloat} <: AbstractVector{T}
    data::NTuple{N,T}

    function Point{N,T}(data::NTuple{N,T}) where {N,T<:AbstractFloat}
        N > 0 ||
            throw(ArgumentError("a point must have at least one component"))
        return new{N,T}(data)
    end
end

Base.size(::Point{N}) where {N} = (N,)
Base.IndexStyle(::Type{<:Point}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(state::Point, index::Int) =
    state.data[index]
Base.show(io::IO, state::Point) = print(io, "Point", state.data)

function Point(values::NTuple{N,<:Real}) where {N}
    floats = map(float, values)
    T = promote_type(typeof.(floats)...)
    return Point{N,T}(convert(NTuple{N,T}, floats))
end

Point() = throw(ArgumentError("a point must have at least one component"))
Point(values::Real...) = Point(Tuple(values))
Point(values::AbstractVector{<:Real}) = Point(Tuple(values))

# Preserve the fixed-size, immutable representation for ordinary vector
# arithmetic instead of falling back to Base's allocated `Vector` results for
# generic `AbstractVector`s.
@inline Base.:+(left::Point{N}, right::Point{N}) where {N} =
    Point(ntuple(index -> left[index] + right[index], Val(N)))

@inline Base.:-(left::Point{N}, right::Point{N}) where {N} =
    Point(ntuple(index -> left[index] - right[index], Val(N)))

@inline Base.:*(factor::Real, point::Point{N}) where {N} =
    Point(ntuple(index -> factor * point[index], Val(N)))

@inline Base.:*(point::Point, factor::Real) = factor * point

@inline Base.:/(point::Point{N}, divisor::Real) where {N} =
    Point(ntuple(index -> point[index] / divisor, Val(N)))

"""
    phase_momentum(state) -> Point
    phase_position(state) -> Point

The two halves of an even-dimensional kinetic phase-space state.

Kinetic models use the `(v, x)` layout: the leading half is the momentum (the
diffusive block) and the trailing half is the position, each of length `n` for a
state of dimension `2n`. Requiring that even dimension and splitting the halves
is the only place the kinetic interpretation enters; `Point` itself is a plain
vector. (`position` is already a `Base` binding, hence the `phase_` prefix.)
"""
@inline function phase_blocks(state::Point{N,T}) where {N,T}
    iseven(N) ||
        throw(ArgumentError("a kinetic state must have even dimension; got $N"))
    half = N ÷ 2
    return (
        Point(ntuple(index -> state[index], Val(N ÷ 2))),
        Point(ntuple(index -> state[half+index], Val(N ÷ 2))),
    )
end

phase_momentum(state::Point) = phase_blocks(state)[1]
phase_position(state::Point) = phase_blocks(state)[2]

"""
Abstract supertype for phase-space dynamics.

A value `dynamics <: AbstractDynamics` participates in a simulation by
providing two methods:

    noise_dimension(dynamics, state) -> Int
    advance(dynamics, state, dt, normal_increment) -> Point

`AbstractDynamics` is a nominal tag for dispatch and documentation, not a
checked trait: as in the rest of Julia, the interface is the generic functions
and conformance is a method, not a declaration.
"""
abstract type AbstractDynamics end

"""
    advance(dynamics, state, dt, normal_increment) -> Point

One Euler–Maruyama step of `dynamics`. For a phase-space state `y`, drift `a`,
diffusion `b`, and driving Wiener processes `W₁, …, W_p`, the step is

    advance(dynamics, y, dt, normal_increment)
        = y + a(y) * dt + b(y) * sqrt(dt) * normal_increment,

where `normal_increment = (z₁, …, z_p)` is a tuple of independent standard
normals with

    zₖ = (Wₖ(t + dt) - Wₖ(t)) / sqrt(dt),    zₖ ~ Normal(0, 1).

The count `p` is prescribed by [`noise_dimension`](@ref); the tuple length is
part of the dynamics protocol, so `advance` can validate it. Standardizing the
increment to unit variance keeps `dt` the only place the time-discretization
scale lives; a deterministic dynamics simply ignores the argument.

The increment is passed in rather than drawn here. That keeps `advance` a pure
map of `(state, dt, normal_increment)`, so all randomness stays in the caller's
per-trajectory stream and a run is reproducible without the dynamics knowing
anything about RNGs.
"""
function advance end

"""
    noise_dimension(dynamics, state) -> Int

Number of independent standard-normal increments one step of `dynamics` at
`state` consumes, i.e. the length of the tuple passed to [`advance`](@ref). The
dynamics prescribes the noise shape here, so the caller knows how many normals
to draw.

The return value must be inferable as a constant — derive it from the state
dimension or a type parameter, not from a runtime field — so the increment tuple
stays type-stable.

By convention a diagonal-noise dynamics drives the leading
`noise_dimension` components of the state — the diffusive block first. The
kinetic models follow this directly: their leading half is the momentum.
"""
function noise_dimension end

"""
    KineticBrownian(; noise=1.0f0)

Kinetic Brownian dynamics on an even-dimensional `(v, x)` state:

`dV = noise dW`, `dX = V dt`.

The state dimension is `2n` with `n = length(v) = length(x)`. The leading `n`
components are the momentum and receive the `n` Wiener increments.
"""
struct KineticBrownian{T<:AbstractFloat} <: AbstractDynamics
    noise::T
end

KineticBrownian(; noise = 1.0f0) = KineticBrownian(noise)

function noise_dimension(::KineticBrownian, ::Point{N,T}) where {N,T}
    iseven(N) || throw(
        ArgumentError(
            "KineticBrownian requires an even-dimensional state; got $N",
        ),
    )
    return N ÷ 2
end

@inline function advance(
    dynamics::KineticBrownian{T},
    state::Point{N,T},
    dt::T,
    normal_increment::NTuple{P,T},
) where {N,T,P}
    half = N ÷ 2
    P == half || throw(
        ArgumentError(
            "KineticBrownian expects $half normal increments; got $P",
        ),
    )
    noise = dynamics.noise
    next_momentum = ntuple(
        index -> state[index] + noise * sqrt(dt) * normal_increment[index],
        Val(N ÷ 2),
    )
    next_position =
        ntuple(index -> state[half+index] + state[index] * dt, Val(N ÷ 2))
    return Point{N,T}((next_momentum..., next_position...))
end

"""
    UnderdampedLangevin(gradient; friction, noise=1.0f0)

Underdamped Langevin dynamics on an even-dimensional `(v, x)` state:

`dV = (-friction*V - gradient(X))dt + noise dW`, `dX = V dt`.

`gradient(position_state)` must return an indexable of length `length(x)`, such
as an `NTuple` or a `Point`. The leading `n` components are the momentum and
receive the `n` Wiener increments.
"""
struct UnderdampedLangevin{T<:AbstractFloat,F} <: AbstractDynamics
    gradient::F
    friction::T
    noise::T
end

function UnderdampedLangevin(gradient; friction, noise = 1.0f0)
    T = promote_type(typeof(float(friction)), typeof(float(noise)))
    return UnderdampedLangevin{T,typeof(gradient)}(
        gradient,
        convert(T, friction),
        convert(T, noise),
    )
end

function noise_dimension(::UnderdampedLangevin, ::Point{N,T}) where {N,T}
    iseven(N) || throw(
        ArgumentError(
            "UnderdampedLangevin requires an even-dimensional state; got $N",
        ),
    )
    return N ÷ 2
end

@inline function advance(
    dynamics::UnderdampedLangevin{T},
    state::Point{N,T},
    dt::T,
    normal_increment::NTuple{P,T},
) where {N,T,P}
    half = N ÷ 2
    P == half || throw(
        ArgumentError(
            "UnderdampedLangevin expects $half normal increments; got $P",
        ),
    )
    velocity, position_state = phase_blocks(state)
    gradient = dynamics.gradient(position_state)
    friction, noise = dynamics.friction, dynamics.noise
    next_velocity = ntuple(
        index ->
            velocity[index] +
            (-friction * velocity[index] - gradient[index]) * dt +
            noise * sqrt(dt) * normal_increment[index],
        Val(N ÷ 2),
    )
    next_position = ntuple(
        index -> position_state[index] + velocity[index] * dt,
        Val(N ÷ 2),
    )
    return Point{N,T}((next_velocity..., next_position...))
end
