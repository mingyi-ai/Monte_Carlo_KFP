# Harmonic measure for underdamped Langevin dynamics with two absorbing wells.
# Optional controls: MCKFP_TRAJECTORIES, MCKFP_DT, MCKFP_MAX_STEPS,
# MCKFP_SEED, MCKFP_ALPHA, MCKFP_RADIUS, and MCKFP_OUTPUT.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using MonteCarloKFP
using Plots

env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float32, get(ENV, name, string(default)))

potential(x) = (x^2 - one(x))^2 / 4 + oftype(x, 0.1) * (x + one(x))
potential_gradient(x) = x * (x^2 - one(x)) + oftype(x, 0.1)

function bisect_root(
    f,
    lower::T,
    upper::T;
    iterations = 64,
) where {T<:AbstractFloat}
    lower_value = f(lower)
    upper_value = f(upper)
    signbit(lower_value) == signbit(upper_value) &&
        throw(ArgumentError("root is not bracketed by ($lower, $upper)"))
    for _ = 1:iterations
        midpoint = (lower + upper) / 2
        midpoint_value = f(midpoint)
        iszero(midpoint_value) && return midpoint
        if signbit(midpoint_value) == signbit(lower_value)
            lower, lower_value = midpoint, midpoint_value
        else
            upper = midpoint
        end
    end
    return (lower + upper) / 2
end

left_minimum = bisect_root(potential_gradient, -1.5f0, -0.5f0)
saddle = bisect_root(potential_gradient, -0.5f0, 0.5f0)
right_minimum = bisect_root(potential_gradient, 0.5f0, 1.5f0)
radius = env_float("MCKFP_RADIUS", 0.35f0)
dt = env_float("MCKFP_DT", 0.008f0)
alpha = env_float("MCKFP_ALPHA", 0.995f0)
friction = (1.0f0 - alpha) / dt

gradient(position) = (potential_gradient(position[1]),)
dynamics = UnderdampedLangevin(gradient; friction = friction, noise = 1.0f0)
domain = TargetDisks(left_minimum, radius, right_minimum, radius)
config = SimulationConfig(
    trajectories = env_int("MCKFP_TRAJECTORIES", 20_000),
    dt = dt,
    max_steps = env_int("MCKFP_MAX_STEPS", 10_000),
    seed = env_int("MCKFP_SEED", 7),
)
result = simulate(dynamics, domain, (0.0f0, saddle), config; progress = true)

println(
    "critical points (left, saddle, right): ",
    (left_minimum, saddle, right_minimum),
)
println("exit rate: ", exit_rate(result))
println(
    "left-target probability among exits: ",
    conditional_exit_probability(result, 1),
)
println(
    "right-target probability among exits: ",
    conditional_exit_probability(result, 2),
)

function target_angles(result, component_index, center_x)
    return [
        mod(atan(point[1], point[2] - center_x), 2π) for (point, index) in
        zip(result.exit_points, result.exit_component_indices) if
        index == component_index
    ]
end

left_angles = target_angles(result, 1, left_minimum)
right_angles = target_angles(result, 2, right_minimum)
circle_angle = range(0, 2π; length = 361)

phase_plot = plot(
    left_minimum .+ radius .* cos.(circle_angle),
    radius .* sin.(circle_angle);
    color = :grey,
    label = "left target",
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "v",
    title = "Absorbing-target exits",
)
plot!(
    phase_plot,
    right_minimum .+ radius .* cos.(circle_angle),
    radius .* sin.(circle_angle);
    color = :grey,
    label = "right target",
)
for (component, label, color) in
    ((1, "left exits", :blue), (2, "right exits", :red))
    indices = findall(==(component), result.exit_component_indices)
    scatter!(
        phase_plot,
        [result.exit_points[index][2] for index in indices],
        [result.exit_points[index][1] for index in indices];
        label = label,
        color = color,
        markersize = 1.5,
        markerstrokewidth = 0,
        alpha = 0.3,
    )
end
scatter!(
    phase_plot,
    [saddle],
    [0.0f0];
    color = :black,
    label = "start",
    markersize = 5,
)

function harmonic_measure_plot(angles, target_name, color)
    panel = plot(;
        xlabel = "boundary angle (radians)",
        ylabel = "probability per bin",
        title = "$target_name harmonic measure",
        xlims = (0, 2π),
        legend = false,
    )
    if !isempty(angles)
        histogram!(
            panel,
            angles;
            bins = range(0, 2π; length = 49),
            normalize = :probability,
            color = color,
            alpha = 0.75,
        )
    end
    return panel
end

left_angle_plot = harmonic_measure_plot(left_angles, "Left target", :blue)
right_angle_plot = harmonic_measure_plot(right_angles, "Right target", :red)

position_range = range(-1.6f0, 1.6f0; length = 300)
potential_plot = plot(
    position_range,
    potential.(position_range);
    xlabel = "x",
    ylabel = "U(x)",
    title = "Double-well potential",
    label = false,
    color = :black,
)
scatter!(
    potential_plot,
    [left_minimum, saddle, right_minimum],
    potential.([left_minimum, saddle, right_minimum]);
    label = "critical points",
    color = [:blue, :black, :red],
)

figure = plot(
    phase_plot,
    potential_plot,
    left_angle_plot,
    right_angle_plot;
    layout = (2, 2),
    size = (1100, 800),
)
default_output =
    joinpath(@__DIR__, "artifacts", "double_well_harmonic_measure.png")
output = get(ENV, "MCKFP_OUTPUT", default_output)
mkpath(dirname(output))
savefig(figure, output)
println("saved plot: ", abspath(output))
