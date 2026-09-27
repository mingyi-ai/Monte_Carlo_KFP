"""Private helpers for spherical boundary components in arbitrary dimensions."""

"""Which side of a supporting sphere the active region occupies."""
@enum SphereSide::Int8 begin
    InsideSphere
    OutsideSphere
end

"""Surface area of an `(n-1)`-sphere: `2 π^{n/2} / Γ(n/2) · r^{n-1}`, via factorials."""
function sphere_surface_area(ambient_dimension::Integer, radius)
    n = Int(ambient_dimension)
    n >= 2 || throw(
        ArgumentError("a sphere needs an ambient dimension of at least 2"),
    )
    if iseven(n)
        half = n ÷ 2
        return 2 * π^half / factorial(half - 1) * radius^(n - 1)
    end
    half = (n - 1) ÷ 2
    return 2 * (4π)^half * factorial(half) / factorial(2 * half) *
           radius^(n - 1)
end

@inline function sphere_level(
    center::Point{N,T},
    radius::T,
    side::SphereSide,
    point::Point{N,T},
) where {N,T}
    radial_level = norm(point - center) - radius
    return side == InsideSphere ? radial_level : -radial_level
end

"""First intersection of `from → to` with the sphere, or `nothing`."""
function sphere_hit(
    center::Point{N,T},
    radius::T,
    from::Point{N,T},
    to::Point{N,T},
) where {N,T}
    direction = to - from
    offset = from - center
    quadratic = dot(direction, direction)
    iszero(quadratic) && return nothing
    linear = 2 * dot(offset, direction)
    constant = dot(offset, offset) - radius^2
    discriminant = linear^2 - 4 * quadratic * constant
    discriminant < zero(discriminant) && return nothing

    root = sqrt(discriminant)
    denominator = 2 * quadratic
    best = nothing
    for candidate in
        ((-linear - root) / denominator, (-linear + root) / denominator)
        zero(candidate) <= candidate <= one(candidate) || continue
        point = Point(
            ntuple(
                axis -> muladd(candidate, direction[axis], from[axis]),
                Val(N),
            ),
        )
        if isnothing(best) || candidate < best[2]
            best = (point, candidate)
        end
    end
    return best
end

"""Point on the sphere from normalized spherical coordinates."""
function sphere_point(
    center::Point{N,T},
    radius::T,
    parameter::Tuple,
) where {N,T}
    M = N - 1
    length(parameter) == M || throw(
        ArgumentError(
            "a sphere in $N dimensions needs $M normalized angles; got $(length(parameter))",
        ),
    )
    angles = ntuple(Val(M)) do index
        full_turn = index < M ? T(π) : T(2π)
        full_turn * mod(T(parameter[index]), one(T))
    end
    coordinates = ntuple(Val(N)) do axis
        prefix = one(T)
        for index = 1:(axis-1)
            prefix *= sin(angles[index])
        end
        if axis < N
            center[axis] + radius * prefix * cos(angles[axis])
        else
            center[axis] + radius * prefix
        end
    end
    return Point(coordinates)
end

sphere_point(center::Point{2,T}, radius::T, parameter::Real) where {T} =
    sphere_point(center, radius, (parameter,))

"""Earliest of a tuple of `(point, fraction)` candidates, with a component index."""
function earliest_hit(hits::Tuple)
    best = nothing
    for (index, candidate) in pairs(hits)
        isnothing(candidate) && continue
        point, fraction = candidate
        if isnothing(best) ||
           fraction < best[2] ||
           (fraction == best[2] && index < best[3])
            best = (point, fraction, index)
        end
    end
    isnothing(best) && return nothing
    return BoundaryHit(best[1], best[3], best[2])
end
