# Monte Carlo solution of a square Dirichlet problem near its singular set.
# Optional controls: MCKFP_TRAJECTORIES, MCKFP_DT, MCKFP_MAX_STEPS,
# MCKFP_SEED, MCKFP_GRID_X, MCKFP_GRID_V, and MCKFP_OUTPUT.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using MonteCarloKFP
using Plots
using ProgressMeter: Progress, next!

env_int(name, default) = parse(Int, get(ENV, name, string(default)))
env_float(name, default) = parse(Float32, get(ENV, name, string(default)))

include("square_boundary.jl")

function observed_boundary_mean(result)
    exited = findall(!iszero, result.exit_component_indices)
    isempty(exited) && return NaN
    total = sum(exited) do index
        point = result.exit_points[index] # package order: (v, x)
        boundary_payoff(point[2], point[1])
    end
    return total / length(exited)
end

function main()
    trajectory_count = env_int("MCKFP_TRAJECTORIES", 25_000)
    dt = env_float("MCKFP_DT", 0.001f0)
    max_steps = env_int("MCKFP_MAX_STEPS", 6_000)
    base_seed = env_int("MCKFP_SEED", 1_000)
    x_count = env_int("MCKFP_GRID_X", 17)
    v_count = env_int("MCKFP_GRID_V", 21)
    x_count >= 3 || error("MCKFP_GRID_X must be at least 3")
    v_count >= 3 || error("MCKFP_GRID_V must be at least 3")

    x_values = collect(range(-1.0f0, 1.0f0; length = x_count))
    v_values = collect(range(-1.0f0, 1.0f0; length = v_count))
    solution = Matrix{Float64}(undef, v_count, x_count)
    exit_rates = Matrix{Float64}(undef, v_count, x_count)
    domain = SquareDomain(1.0f0)
    dynamics = KineticBrownian(noise = 1.0f0)

    run_index = 0
    grid_progress = Progress(x_count * v_count; desc = "Grid points:", dt = 0)
    for (v_index, v) in pairs(v_values), (x_index, x) in pairs(x_values)
        prescribed_value = prescribed_boundary_value(x, v)
        if !isnothing(prescribed_value)
            solution[v_index, x_index] = prescribed_value
            exit_rates[v_index, x_index] = 1.0
        else
            run_index += 1
            config = SimulationConfig(
                trajectories = trajectory_count,
                dt = dt,
                max_steps = max_steps,
                seed = base_seed + run_index,
            )
            trajectory_progress =
                Progress(trajectory_count; desc = "  Trajectories:", offset = 1)
            result = simulate(
                dynamics,
                domain,
                simulation_start(x, v),
                config;
                progress = trajectory_progress,
            )
            solution[v_index, x_index] = observed_boundary_mean(result)
            exit_rates[v_index, x_index] = exit_rate(result)
        end
        next!(grid_progress)
    end

    finite_values = filter(isfinite, vec(solution))
    interior_exit_rate =
        sum(exit_rates[2:(end-1), 2:(end-1)]) / ((v_count - 2) * (x_count - 2))
    println("simulated grid points: ", run_index)
    println("mean interior exit rate: ", interior_exit_rate)
    println(
        "finite solution range: ",
        isempty(finite_values) ? "none" : extrema(finite_values),
    )

    heatmap_plot = heatmap(
        x_values,
        v_values,
        solution;
        xlabel = "x",
        ylabel = "v",
        title = "Estimated Dirichlet solution",
        color = :viridis,
    )
    surface_plot = surface(
        x_values,
        v_values,
        solution;
        xlabel = "x",
        ylabel = "v",
        zlabel = "u(x, v)",
        title = "Square-boundary experiment",
        color = :viridis,
        camera = (35, 35),
    )
    exit_plot = heatmap(
        x_values,
        v_values,
        exit_rates;
        xlabel = "x",
        ylabel = "v",
        title = "Finite-horizon exit rate",
        color = :magma,
        clims = (0, 1),
    )

    figure = plot(
        heatmap_plot,
        surface_plot,
        exit_plot;
        layout = (1, 3),
        size = (1450, 450),
    )
    default_output = joinpath(@__DIR__, "artifacts", "square_dirichlet.png")
    output = get(ENV, "MCKFP_OUTPUT", default_output)
    mkpath(dirname(output))
    savefig(figure, output)
    println("saved plot: ", abspath(output))
end

main()
