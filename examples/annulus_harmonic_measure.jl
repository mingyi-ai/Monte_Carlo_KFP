# Harmonic measure for kinetic Brownian motion in an annulus.
# Optional controls: MCKFP_TRAJECTORIES, MCKFP_DT, MCKFP_MAX_STEPS,
# MCKFP_SEED, and MCKFP_OUTPUT.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using MonteCarloKFP
using Plots

env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float32, get(ENV, name, string(default)))

const INNER_RADIUS = 2.0f0
const OUTER_RADIUS = 4.0f0
const START = (-1.0f0, 3.0f0) # package order: (v, x)

config = SimulationConfig(
    trajectories = env_int("MCKFP_TRAJECTORIES", 10_000),
    dt = env_float("MCKFP_DT", 0.002f0),
    max_steps = env_int("MCKFP_MAX_STEPS", 10_000),
    seed = env_int("MCKFP_SEED", 42),
)
domain = AnnulusDomain(INNER_RADIUS, OUTER_RADIUS)
result = simulate(
    KineticBrownian(noise = 1.0f0),
    domain,
    START,
    config;
    progress = true,
)

inner_probability = conditional_exit_probability(result, 1)
outer_probability = conditional_exit_probability(result, 2)
println("trajectories: ", length(result))
println("exit rate: ", exit_rate(result))
println("inner-boundary probability among exits: ", inner_probability)
println("outer-boundary probability among exits: ", outer_probability)

function exit_angles(result, component_index)
    return [
        mod(atan(point[1], point[2]), 2π) for (point, index) in
        zip(result.exit_points, result.exit_component_indices) if
        index == component_index
    ]
end

inner_angles = exit_angles(result, 1)
outer_angles = exit_angles(result, 2)

circle_angle = range(0, 2π; length = 361)
phase_plot = plot(
    OUTER_RADIUS .* cos.(circle_angle),
    OUTER_RADIUS .* sin.(circle_angle);
    color = :grey,
    label = "outer boundary",
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "v",
    title = "Exit locations",
)
plot!(
    phase_plot,
    INNER_RADIUS .* cos.(circle_angle),
    INNER_RADIUS .* sin.(circle_angle);
    color = :grey,
    label = "inner boundary",
)
for (component, label, color) in
    ((1, "inner exits", :blue), (2, "outer exits", :red))
    indices = findall(==(component), result.exit_component_indices)
    scatter!(
        phase_plot,
        [result.exit_points[index][2] for index in indices],
        [result.exit_points[index][1] for index in indices];
        markersize = 1.5,
        markerstrokewidth = 0,
        alpha = 0.35,
        label = label,
        color = color,
    )
end
scatter!(
    phase_plot,
    [START[2]],
    [START[1]];
    color = :black,
    label = "start",
    markersize = 5,
)

function harmonic_measure_plot(angles, boundary_name, color)
    panel = plot(;
        xlabel = "boundary angle (radians)",
        ylabel = "probability per bin",
        title = "$boundary_name harmonic measure",
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

inner_angle_plot = harmonic_measure_plot(inner_angles, "Inner boundary", :blue)
outer_angle_plot = harmonic_measure_plot(outer_angles, "Outer boundary", :red)
figure = plot(
    phase_plot,
    inner_angle_plot,
    outer_angle_plot;
    layout = (1, 3),
    size = (1450, 450),
)
default_output = joinpath(@__DIR__, "artifacts", "annulus_harmonic_measure.png")
output = get(ENV, "MCKFP_OUTPUT", default_output)
mkpath(dirname(output))
savefig(figure, output)
println("saved plot: ", abspath(output))
