"""Private helpers for finite line-segment boundary components in two dimensions."""

@inline cross2(left::Point{2}, right::Point{2}) =
    left[1] * right[2] - left[2] * right[1]

"""First intersection of a path segment with a finite boundary segment."""
function segment_hit(
    boundary_from::Point{2},
    boundary_to::Point{2},
    from::Point{2},
    to::Point{2},
)
    direction = to - from
    boundary_direction = boundary_to - boundary_from
    denominator = cross2(direction, boundary_direction)
    iszero(denominator) && return nothing

    offset = boundary_from - from
    fraction = cross2(offset, boundary_direction) / denominator
    boundary_fraction = cross2(offset, direction) / denominator
    zero(fraction) <= fraction <= one(fraction) || return nothing
    zero(boundary_fraction) <= boundary_fraction <= one(boundary_fraction) ||
        return nothing

    point = Point(
        ntuple(axis -> muladd(fraction, direction[axis], from[axis]), Val(2)),
    )
    return point, fraction
end

"""Earliest hit among an ordered tuple of finite boundary segments."""
function earliest_segment_hit(segments::Tuple, from::Point{2}, to::Point{2})
    best = nothing
    for (component, segment) in pairs(segments)
        candidate = segment_hit(segment[1], segment[2], from, to)
        isnothing(candidate) && continue
        point, fraction = candidate
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

@inline segment_measure(segment::Tuple{Point{2,T},Point{2,T}}) where {T} =
    norm(segment[2] - segment[1])

function segment_point(
    segment::Tuple{Point{2,T},Point{2,T}},
    parameter,
) where {T}
    fraction = mod(T(parameter), one(T))
    return segment[1] + fraction * (segment[2] - segment[1])
end

function point_segment_distance(
    point::Point{2},
    segment::Tuple{Point{2,T},Point{2,T}},
) where {T}
    from, to = segment
    direction = to - from
    offset = point - from
    raw_fraction = dot(offset, direction) / dot(direction, direction)
    fraction = clamp(raw_fraction, zero(raw_fraction), one(raw_fraction))
    return norm(offset - fraction * direction)
end
