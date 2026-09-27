bump(value) = abs(value) < 1 ? exp(-1 / (1 - value^2)) : zero(value)

on_boundary_coordinate(value, boundary) =
    isapprox(value, boundary; atol = 2eps(typeof(value)), rtol = 0)

"""
Prescribed data on the regular/outgoing boundary in conventional `(x, v)`
coordinates. Return `nothing` on the interior and the transport-influx boundary.
"""
function prescribed_boundary_value(x, v)
    if on_boundary_coordinate(abs(v), one(v))
        return zero(promote_type(typeof(x), typeof(v)))
    elseif on_boundary_coordinate(x, -one(x)) && v < zero(v)
        return bump(2 * (v + 0.5f0))
    elseif on_boundary_coordinate(x, one(x)) && v > zero(v)
        return bump(2 * (v - 0.5f0))
    end
    return nothing
end

"""
Evaluate prescribed data at a discrete exit. Exact trajectories leave a side
only through its outgoing half. An Euler step can cross with a small wrong-sign
interpolated velocity near grazing, so project that velocity to the grazing
point rather than treating the influx boundary as prescribed data.
"""
function boundary_payoff(x, v)
    if on_boundary_coordinate(abs(v), one(v))
        return zero(promote_type(typeof(x), typeof(v)))
    elseif on_boundary_coordinate(x, -one(x))
        return bump(2 * (min(v, zero(v)) + 0.5f0))
    elseif on_boundary_coordinate(x, one(x))
        return bump(2 * (max(v, zero(v)) - 0.5f0))
    end
    error("exit point (x, v) = ($x, $v) is not on the square boundary")
end

"Move an influx or grazing boundary point just inside for its one-sided trace."
function simulation_start(x, v)
    inset = 32eps(typeof(x))
    if on_boundary_coordinate(x, -one(x)) && v >= zero(v)
        x += inset
    elseif on_boundary_coordinate(x, one(x)) && v <= zero(v)
        x -= inset
    end
    return (v, x) # package order
end
