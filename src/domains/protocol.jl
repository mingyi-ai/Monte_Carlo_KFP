"""
Domain protocol
===============

A domain is a subset of `N`-dimensional phase space with a signed level that is
negative inside. The protocol is a set of generic functions implemented directly
by each concrete domain type. It deliberately does not prescribe patches, meshes,
constructive solid geometry, or another boundary representation.

The protocol primitives are:

    domain_dimension(domain)              -> Int
    domain_components(domain)             -> Tuple of Symbols
    domain_scale(domain)                  -> Real
    domain_level(domain, point)           -> Real
    component_level(domain, point, index) -> Real
    domain_hit(domain, from, to)          -> Union{Nothing, BoundaryHit}
    component_measure(domain, index)      -> Real
    component_point(domain, index, parameter) -> Point

Everything user-facing — `point_location`, `boundary_component`,
`first_boundary_hit`, `boundary_measure`, `boundary_point`, `boundary_points` —
is implemented once on top of these, so a new domain shape only describes its
own geometry. Julia has no trait bounds; as with `AbstractDynamics`, the
interface is the generic functions and conformance is a set of methods. The
provided square, box, sphere-based, and concave L-shaped domains all implement
these primitives directly. More complex mesh and arbitrary high-dimensional
geometry are outside the package's current scope.
"""
@enum PointLocation::UInt8 begin
    Interior
    OnBoundary
    Exterior
end

"""
The first boundary intersection along a path segment.

`component_index` is the one-based index into `boundary_components(domain)`, and
`segment_fraction` locates the hit along `from → to` in `[0, 1]`.
"""
struct BoundaryHit{N,T<:AbstractFloat}
    point::Point{N,T}
    component_index::Int
    segment_fraction::T
end

# --- domain protocol --------------------------------------------------------

"""
Abstract supertype for domains.

A domain is an `N`-dimensional region described by the primitives listed in the
module docstring. Concrete shapes implement those primitives directly.
"""
abstract type AbstractDomain end

"""Ambient dimension `N` of the domain."""
function domain_dimension end

"""Tuple of stable boundary-component names."""
function domain_components end

"""Characteristic length used to scale the default classification tolerance."""
function domain_scale end

"""Signed level of the whole domain; negative inside, zero on the boundary."""
function domain_level end

"""Signed level of boundary component `index`."""
function component_level end

"""First boundary crossing along `from → to` as a `BoundaryHit`, or `nothing`."""
function domain_hit end

"""`(N-1)`-dimensional measure of boundary component `index`."""
function component_measure end

"""Point on boundary component `index` from normalized `parameter` in `[0, 1)`."""
function component_point end

# --- generic domain operations ---------------------------------------------

point_from(point::Point) = point
function point_from(point)
    values = Tuple(point)
    isempty(values) &&
        throw(ArgumentError("a point must have at least one component"))
    return Point(values)
end

"""Coerce a point to the domain's dimension, rejecting other dimensions."""
function domain_point(domain::AbstractDomain, value)
    point = point_from(value)
    dimension = domain_dimension(domain)
    length(point) == dimension || throw(
        ArgumentError(
            "this domain is $(dimension)-dimensional; got a $(length(point))-dimensional point",
        ),
    )
    return point
end

"""Stable names of a domain's boundary components."""
boundary_components(domain::AbstractDomain) = domain_components(domain)

function component_index(domain::AbstractDomain, component::Integer)
    1 <= component <= length(domain_components(domain)) ||
        throw(ArgumentError("invalid boundary component: $component"))
    return Int(component)
end

function component_index(domain::AbstractDomain, component::Symbol)
    index = findfirst(==(component), domain_components(domain))
    isnothing(index) &&
        throw(ArgumentError("unknown boundary component: $component"))
    return index
end

function default_boundary_tolerance(domain::AbstractDomain, point::Point)
    scale = domain_scale(domain)
    T = promote_type(typeof(point[1]), typeof(scale))
    return T(16) * eps(T) * max(one(T), T(scale))
end

"""Classify a point as `Interior`, `OnBoundary`, or `Exterior`."""
function point_location(domain::AbstractDomain, point; atol = nothing)
    validated_point = domain_point(domain, point)
    tolerance =
        isnothing(atol) ? default_boundary_tolerance(domain, validated_point) :
        atol
    tolerance >= 0 || throw(ArgumentError("atol must be non-negative"))

    level = domain_level(domain, validated_point)
    level > tolerance && return Exterior
    return level >= -tolerance ? OnBoundary : Interior
end

"""Calling a domain is shorthand for `point_location(domain, point)`."""
(domain::AbstractDomain)(point; atol = nothing) =
    point_location(domain, point; atol = atol)

"""Signed level of one boundary component: negative inside, zero on it."""
function boundary_level(domain::AbstractDomain, point, component)
    validated_point = domain_point(domain, point)
    return component_level(
        domain,
        validated_point,
        component_index(domain, component),
    )
end

"""Return the boundary-component ID for a point, or zero if it is not on the boundary."""
function boundary_component(domain::AbstractDomain, point; atol = nothing)
    validated_point = domain_point(domain, point)
    tolerance =
        isnothing(atol) ? default_boundary_tolerance(domain, validated_point) :
        atol
    tolerance >= 0 || throw(ArgumentError("atol must be non-negative"))

    domain_level(domain, validated_point) > tolerance && return 0
    for index = 1:length(domain_components(domain))
        abs(component_level(domain, validated_point, index)) <= tolerance &&
            return index
    end
    return 0
end

"""Measure of one boundary component."""
function boundary_measure(domain::AbstractDomain, component)
    return component_measure(domain, component_index(domain, component))
end

"""Point at normalized coordinates on one boundary component."""
function boundary_point(domain::AbstractDomain, component, parameter)
    return component_point(
        domain,
        component_index(domain, component),
        parameter,
    )
end

"""
    boundary_points(domain, count; component=:all)

Generate deterministic equal-arclength samples on a two-dimensional domain. All
components are concatenated in protocol order unless one component ID or name is
selected. Equal-arclength sampling is a one-dimensional concept, so this is only
defined for `domain_dimension(domain) == 2`; higher-dimensional boundaries are
sampled through `boundary_point` with normalized coordinates.
"""
function boundary_points(
    domain::AbstractDomain,
    count::Integer;
    component = :all,
)
    count >= 0 || throw(ArgumentError("count must be non-negative"))
    domain_dimension(domain) == 2 || throw(
        ArgumentError(
            "equal-arclength sampling is defined for two-dimensional domains; got a " *
            "$(domain_dimension(domain))-dimensional domain",
        ),
    )
    point_type = typeof(boundary_point(domain, 1, 0.0))
    count == 0 && return point_type[]

    if component !== :all
        index = component_index(domain, component)
        return [boundary_point(domain, index, (i - 1) / count) for i = 1:count]
    end

    measures = [
        float(component_measure(domain, index)) for
        index = 1:length(domain_components(domain))
    ]
    starts = [zero(eltype(measures)); cumsum(measures[1:(end-1)])]
    spacing = sum(measures) / count
    points = point_type[]
    sizehint!(points, count)

    for i = 1:count
        distance = (i - 1) * spacing
        index = searchsortedlast(starts, distance)
        parameter = (distance - starts[index]) / measures[index]
        push!(points, boundary_point(domain, index, parameter))
    end
    return points
end

"""Return the first boundary intersection along `from → to`, or `nothing`."""
function first_boundary_hit(domain::AbstractDomain, from, to)
    return domain_hit(
        domain,
        domain_point(domain, from),
        domain_point(domain, to),
    )
end
