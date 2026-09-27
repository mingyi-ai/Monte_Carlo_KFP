"""
Concrete analytic domains.

Each shape implements the domain protocol directly. Segment and sphere helpers
share only local calculations; they do not define interior semantics or expose a
public patch-composition API.
"""

# --- square -----------------------------------------------------------------

"""
    SquareDomain(halfwidth=1.0f0)

Two-dimensional open square `(-h, h)²` with boundary components in the order
`:bottom`, `:right`, `:top`, `:left`.
"""
struct SquareDomain{T<:AbstractFloat} <: AbstractDomain
    halfwidth::T

    function SquareDomain{T}(halfwidth::T) where {T<:AbstractFloat}
        isfinite(halfwidth) && halfwidth > zero(T) ||
            throw(ArgumentError("halfwidth must be positive and finite"))
        return new{T}(halfwidth)
    end
end

SquareDomain(halfwidth::Real = 1.0f0) = begin
    h = float(halfwidth)
    SquareDomain{typeof(h)}(h)
end

function square_segments(domain::SquareDomain{T}) where {T}
    h = domain.halfwidth
    return (
        (Point(-h, -h), Point(h, -h)),
        (Point(h, -h), Point(h, h)),
        (Point(h, h), Point(-h, h)),
        (Point(-h, h), Point(-h, -h)),
    )
end

domain_dimension(::SquareDomain) = 2
domain_components(::SquareDomain) = (:bottom, :right, :top, :left)
domain_scale(domain::SquareDomain) = domain.halfwidth

domain_level(domain::SquareDomain, point::Point) =
    max(abs(point[1]) - domain.halfwidth, abs(point[2]) - domain.halfwidth)

function component_level(domain::SquareDomain, point::Point, index::Int)
    h = domain.halfwidth
    index == 1 && return -point[2] - h
    index == 2 && return point[1] - h
    index == 3 && return point[2] - h
    return -point[1] - h
end

domain_hit(domain::SquareDomain, from::Point{2}, to::Point{2}) =
    earliest_segment_hit(square_segments(domain), from, to)

component_measure(domain::SquareDomain, index::Int) =
    segment_measure(square_segments(domain)[index])
component_point(domain::SquareDomain, index::Int, parameter) =
    segment_point(square_segments(domain)[index], parameter)

# --- axis-aligned box -------------------------------------------------------

"""
    BoxDomain(lower, upper)

Axis-aligned `N`-dimensional box `lower[i] < point[i] < upper[i]`. Boundary
components are ordered `:lower_1`, `:upper_1`, ..., `:lower_N`, `:upper_N`.
"""
struct BoxDomain{T<:AbstractFloat,N} <: AbstractDomain
    lower::Point{N,T}
    upper::Point{N,T}

    function BoxDomain{T,N}(lower::Point{N,T}, upper::Point{N,T}) where {T,N}
        N >= 2 ||
            throw(ArgumentError("a box domain needs at least two dimensions"))
        all(index -> lower[index] < upper[index], 1:N) || throw(
            ArgumentError("box lower bounds must be smaller than upper bounds"),
        )
        all(isfinite, lower) && all(isfinite, upper) ||
            throw(ArgumentError("box bounds must be finite"))
        return new{T,N}(lower, upper)
    end
end

function BoxDomain(lower::NTuple{N,<:Real}, upper::NTuple{N,<:Real}) where {N}
    T = promote_type(float(eltype(lower)), float(eltype(upper)))
    lo = Point(ntuple(index -> T(lower[index]), Val(N)))
    hi = Point(ntuple(index -> T(upper[index]), Val(N)))
    return BoxDomain{T,N}(lo, hi)
end

domain_dimension(::BoxDomain{T,N}) where {T,N} = N
domain_components(::BoxDomain{T,N}) where {T,N} = ntuple(Val(2N)) do index
    axis = (index + 1) ÷ 2
    Symbol(iseven(index) ? "upper_" : "lower_", axis)
end

domain_scale(domain::BoxDomain{T,N}) where {T,N} = maximum(
    ntuple(
        index -> max(abs(domain.lower[index]), abs(domain.upper[index])),
        Val(N),
    ),
)

function domain_level(domain::BoxDomain{T,N}, point::Point) where {T,N}
    return maximum(
        ntuple(
            axis -> max(
                domain.lower[axis] - point[axis],
                point[axis] - domain.upper[axis],
            ),
            Val(N),
        ),
    )
end

function component_level(domain::BoxDomain, point::Point, index::Int)
    axis = (index + 1) ÷ 2
    return iseven(index) ? point[axis] - domain.upper[axis] :
           domain.lower[axis] - point[axis]
end

function domain_hit(
    domain::BoxDomain{T,N},
    from::Point{N},
    to::Point{N},
) where {T,N}
    direction = to - from
    best = nothing
    for component = 1:2N
        axis = (component + 1) ÷ 2
        iszero(direction[axis]) && continue
        boundary = iseven(component) ? domain.upper[axis] : domain.lower[axis]
        fraction = (boundary - from[axis]) / direction[axis]
        zero(fraction) <= fraction <= one(fraction) || continue
        point = Point(
            ntuple(
                index -> muladd(fraction, direction[index], from[index]),
                Val(N),
            ),
        )
        all(
            index ->
                index == axis ||
                domain.lower[index] <= point[index] <= domain.upper[index],
            1:N,
        ) || continue
        if isnothing(best) ||
           fraction < best.segment_fraction ||
           (
               fraction == best.segment_fraction &&
               component < best.component_index
           )
            best = BoundaryHit(point, component, fraction)
        end
    end
    return best
end

function component_measure(domain::BoxDomain{T,N}, index::Int) where {T,N}
    axis = (index + 1) ÷ 2
    return prod(
        ntuple(
            position -> begin
                target = position < axis ? position : position + 1
                domain.upper[target] - domain.lower[target]
            end,
            Val(N - 1),
        ),
    )
end

function component_point(
    domain::BoxDomain{T,N},
    index::Int,
    parameter::Tuple,
) where {T,N}
    length(parameter) == N - 1 || throw(
        ArgumentError(
            "a box face in $N dimensions needs $(N - 1) normalized coordinates; " *
            "got $(length(parameter))",
        ),
    )
    axis = (index + 1) ÷ 2
    boundary = iseven(index) ? domain.upper[axis] : domain.lower[axis]
    return Point(
        ntuple(Val(N)) do target
            target == axis && return boundary
            position = target < axis ? target : target - 1
            fraction = mod(T(parameter[position]), one(T))
            domain.lower[target] +
            fraction * (domain.upper[target] - domain.lower[target])
        end,
    )
end

component_point(domain::BoxDomain{T,2}, index::Int, parameter::Real) where {T} =
    component_point(domain, index, (parameter,))

# --- concave L shape --------------------------------------------------------

"""
    LShapeDomain(halfwidth=1.0f0)

The open square `(-h, h)²` with its top-right quadrant removed. The resulting
concave L has components `:bottom`, `:right`, `:notch_bottom`, `:notch_left`,
`:top`, and `:left`, ordered counterclockwise around the active region.
"""
struct LShapeDomain{T<:AbstractFloat} <: AbstractDomain
    halfwidth::T

    function LShapeDomain{T}(halfwidth::T) where {T<:AbstractFloat}
        isfinite(halfwidth) && halfwidth > zero(T) ||
            throw(ArgumentError("halfwidth must be positive and finite"))
        return new{T}(halfwidth)
    end
end

LShapeDomain(halfwidth::Real = 1.0f0) = begin
    h = float(halfwidth)
    LShapeDomain{typeof(h)}(h)
end

function lshape_segments(domain::LShapeDomain{T}) where {T}
    h, z = domain.halfwidth, zero(T)
    return (
        (Point(-h, -h), Point(h, -h)),
        (Point(h, -h), Point(h, z)),
        (Point(h, z), Point(z, z)),
        (Point(z, z), Point(z, h)),
        (Point(z, h), Point(-h, h)),
        (Point(-h, h), Point(-h, -h)),
    )
end

@inline function lshape_inside(domain::LShapeDomain, point::Point)
    h = domain.halfwidth
    return -h < point[1] < h &&
           -h < point[2] < h &&
           (point[1] < zero(point[1]) || point[2] < zero(point[2]))
end

domain_dimension(::LShapeDomain) = 2
domain_components(::LShapeDomain) =
    (:bottom, :right, :notch_bottom, :notch_left, :top, :left)
domain_scale(domain::LShapeDomain) = domain.halfwidth

function domain_level(domain::LShapeDomain, point::Point)
    segments = lshape_segments(domain)
    distance = point_segment_distance(point, segments[1])
    for index = 2:length(segments)
        distance = min(distance, point_segment_distance(point, segments[index]))
    end
    iszero(distance) && return distance
    return lshape_inside(domain, point) ? -distance : distance
end

function component_level(domain::LShapeDomain, point::Point, index::Int)
    distance = point_segment_distance(point, lshape_segments(domain)[index])
    iszero(distance) && return distance
    return lshape_inside(domain, point) ? -distance : distance
end

domain_hit(domain::LShapeDomain, from::Point{2}, to::Point{2}) =
    earliest_segment_hit(lshape_segments(domain), from, to)

component_measure(domain::LShapeDomain, index::Int) =
    segment_measure(lshape_segments(domain)[index])
component_point(domain::LShapeDomain, index::Int, parameter) =
    segment_point(lshape_segments(domain)[index], parameter)

# --- sphere-based shapes --------------------------------------------------

"""
A ball: the region inside a sphere (`InsideSphere`) or outside it
(`OutsideSphere`). One boundary component, `:sphere`.
"""
struct BallDomain{T<:AbstractFloat,N} <: AbstractDomain
    center::Point{N,T}
    radius::T
    domain_side::SphereSide

    function BallDomain{T,N}(
        center::Point{N,T},
        radius::T,
        domain_side::SphereSide,
    ) where {T<:AbstractFloat,N}
        N >= 2 || throw(
            ArgumentError("a ball domain needs a dimension of at least 2"),
        )
        isfinite(radius) && radius > zero(T) ||
            throw(ArgumentError("a ball radius must be positive and finite"))
        return new{T,N}(center, radius, domain_side)
    end
end

function BallDomain(
    center,
    radius::Real,
    domain_side::SphereSide = InsideSphere,
)
    center_point = point_from(center)
    N = length(center_point)
    T = promote_type(typeof(center_point[1]), typeof(float(radius)))
    return BallDomain{T,N}(
        Point(ntuple(index -> T(center_point[index]), Val(N))),
        T(radius),
        domain_side,
    )
end

domain_dimension(::BallDomain{T,N}) where {T,N} = N
domain_components(::BallDomain) = (:sphere,)
domain_scale(domain::BallDomain{T}) where {T} =
    max(one(T), norm(domain.center) + domain.radius)
component_level(domain::BallDomain, point::Point, ::Int) =
    sphere_level(domain.center, domain.radius, domain.domain_side, point)
domain_level(domain::BallDomain, point::Point) =
    sphere_level(domain.center, domain.radius, domain.domain_side, point)
component_measure(domain::BallDomain{T,N}, ::Int) where {T,N} =
    sphere_surface_area(N, domain.radius)
component_point(domain::BallDomain, ::Int, parameter) =
    sphere_point(domain.center, domain.radius, parameter)

function domain_hit(domain::BallDomain, from::Point, to::Point)
    candidate = sphere_hit(domain.center, domain.radius, from, to)
    isnothing(candidate) && return nothing
    point, fraction = candidate
    return BoundaryHit(point, 1, fraction)
end

"""
A shell between two concentric spheres: outside the inner one and inside the
outer one. Boundary components `:inner` and `:outer`.
"""
struct ShellDomain{T<:AbstractFloat,N} <: AbstractDomain
    center::Point{N,T}
    inner::T
    outer::T

    function ShellDomain{T,N}(
        center::Point{N,T},
        inner::T,
        outer::T,
    ) where {T<:AbstractFloat,N}
        N >= 2 || throw(
            ArgumentError("a shell domain needs a dimension of at least 2"),
        )
        zero(T) < inner < outer || throw(
            ArgumentError("radii must satisfy 0 < inner_radius < outer_radius"),
        )
        return new{T,N}(center, inner, outer)
    end
end

function ShellDomain(center, inner_radius::Real, outer_radius::Real)
    center_point = point_from(center)
    N = length(center_point)
    inner, outer = promote(float(inner_radius), float(outer_radius))
    T = promote_type(typeof(center_point[1]), typeof(inner))
    return ShellDomain{T,N}(
        Point(ntuple(index -> T(center_point[index]), Val(N))),
        T(inner),
        T(outer),
    )
end

domain_dimension(::ShellDomain{T,N}) where {T,N} = N
domain_components(::ShellDomain) = (:inner, :outer)
domain_scale(domain::ShellDomain{T}) where {T} =
    norm(domain.center) + domain.outer

function component_level(domain::ShellDomain, point::Point, index::Int)
    distance = norm(point - domain.center)
    return index == 1 ? domain.inner - distance : distance - domain.outer
end

function domain_level(domain::ShellDomain, point::Point)
    distance = norm(point - domain.center)
    return max(domain.inner - distance, distance - domain.outer)
end

component_measure(domain::ShellDomain{T,N}, index::Int) where {T,N} =
    sphere_surface_area(N, index == 1 ? domain.inner : domain.outer)
component_point(domain::ShellDomain, index::Int, parameter) = sphere_point(
    domain.center,
    index == 1 ? domain.inner : domain.outer,
    parameter,
)

function domain_hit(domain::ShellDomain, from::Point, to::Point)
    return earliest_hit((
        sphere_hit(domain.center, domain.inner, from, to),
        sphere_hit(domain.center, domain.outer, from, to),
    ))
end

"""
The active region outside two disjoint absorbing spheres. Boundary components
`:left` and `:right`.
"""
struct TargetDisks{T<:AbstractFloat,N} <: AbstractDomain
    centers::NTuple{2,Point{N,T}}
    radii::NTuple{2,T}

    function TargetDisks{T,N}(
        centers::NTuple{2,Point{N,T}},
        radii::NTuple{2,T},
    ) where {T<:AbstractFloat,N}
        N >= 2 || throw(
            ArgumentError("a target domain needs a dimension of at least 2"),
        )
        all(radius -> isfinite(radius) && radius > zero(T), radii) ||
            throw(ArgumentError("target radii must be positive and finite"))
        return new{T,N}(centers, radii)
    end
end

function TargetDisks(centers::NTuple{2}, radii::NTuple{2,<:Real})
    first_center = point_from(centers[1])
    N = length(first_center)
    first_radius, second_radius = promote(float(radii[1]), float(radii[2]))
    T = promote_type(typeof(first_center[1]), typeof(first_radius))
    promoted = ntuple(
        index -> begin
            point = point_from(centers[index])
            length(point) == N || throw(
                ArgumentError("both target centers must have $N components"),
            )
            Point(ntuple(axis -> T(point[axis]), Val(N)))
        end,
        Val(2),
    )
    return TargetDisks{T,N}(promoted, (T(first_radius), T(second_radius)))
end

"""
    TargetDisks(left_x, left_radius, right_x, right_radius)

Two-dimensional convenience: disks centered on the position axis (`v = 0`).
"""
function TargetDisks(
    left_x::Real,
    left_radius::Real,
    right_x::Real,
    right_radius::Real,
)
    left, left_r, right, right_r = promote(
        float(left_x),
        float(left_radius),
        float(right_x),
        float(right_radius),
    )
    left_r > 0 || throw(ArgumentError("left_radius must be positive"))
    right_r > 0 || throw(ArgumentError("right_radius must be positive"))
    left < right || throw(ArgumentError("left_x must be less than right_x"))
    left + left_r < right - right_r ||
        throw(ArgumentError("target disks must be disjoint and non-tangent"))
    centers = (Point(zero(left), left), Point(zero(right), right))
    return TargetDisks(centers, (left_r, right_r))
end

domain_dimension(::TargetDisks{T,N}) where {T,N} = N
domain_components(::TargetDisks) = (:left, :right)
domain_scale(domain::TargetDisks{T}) where {T} = max(
    norm(domain.centers[1]) + domain.radii[1],
    norm(domain.centers[2]) + domain.radii[2],
)

component_level(domain::TargetDisks, point::Point, index::Int) = sphere_level(
    domain.centers[index],
    domain.radii[index],
    OutsideSphere,
    point,
)
domain_level(domain::TargetDisks, point::Point) =
    max(component_level(domain, point, 1), component_level(domain, point, 2))
component_measure(domain::TargetDisks{T,N}, index::Int) where {T,N} =
    sphere_surface_area(N, domain.radii[index])
component_point(domain::TargetDisks, index::Int, parameter) =
    sphere_point(domain.centers[index], domain.radii[index], parameter)

function domain_hit(domain::TargetDisks, from::Point, to::Point)
    return earliest_hit((
        sphere_hit(domain.centers[1], domain.radii[1], from, to),
        sphere_hit(domain.centers[2], domain.radii[2], from, to),
    ))
end

# --- convenience shapes ----------------------------------------------------

"""Two-dimensional annulus centered at the origin: a `ShellDomain`."""
function AnnulusDomain(inner_radius::Real, outer_radius::Real)
    inner, outer = promote(float(inner_radius), float(outer_radius))
    zero(inner) < inner < outer || throw(
        ArgumentError("radii must satisfy 0 < inner_radius < outer_radius"),
    )
    return ShellDomain(Point(zero(inner), zero(inner)), inner, outer)
end
