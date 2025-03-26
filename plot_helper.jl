module PlotHelper

using Plots, StatsBase, LinearAlgebra, Colors

export histogram_annulus, plot_boundary_value, plot3d_square

function histogram_annulus(start_point::Vector{Float32}, exit_points::Array{Float32, 2}, inner_radius::Float32, outer_radius::Float32; epsilon=1e-2, num_bins=360, scale_factor=1.0)
    # Helper function to compute angles
    function angles_from_points(points)
        return [mod(atan(p[2], p[1]), 2π) for p in points]
    end

    bar_color = RGB(153/255, 0, 0)  # Strong Red
    grey_color = RGB(230/255, 230/255, 230/255)  # Soft Grey

    # Extract inner and outer exits
    outer_exits = [exit_points[:, i] for i in 1:size(exit_points, 2) if abs(norm(exit_points[:, i]) - outer_radius) < epsilon]
    inner_exits = [exit_points[:, i] for i in 1:size(exit_points, 2) if abs(norm(exit_points[:, i]) - inner_radius) < epsilon]

    # Convert to angles
    angles_outer = angles_from_points(outer_exits)
    angles_inner = angles_from_points(inner_exits)

    # Define x-axis ticks for angle
    ticks_angles = [0, π/2, π, 3π/2, 2π]
    ticks_labels = ["0", "π/2", "π", "3π/2", "2π"]

    plt_hist_annular_outer = histogram(angles_outer, bins=120,
            xlabel="Angle (radians)", ylabel="Frequency", 
            title="Angular Histogram for Outer Boundary",
            label="Frequency", legend=:topright,
            xticks=(ticks_angles, ticks_labels))
    plt_hist_annular_inner = histogram(angles_inner, bins=120,
            xlabel="Angle (radians)", ylabel="Frequency", 
            title="Angular Histogram for Inner Boundary",
            label="Frequency", legend=:topleft,
            xticks=(ticks_angles, ticks_labels))

    # Compute histograms for circular plot
    bin_edges = collect(range(0, stop=2π, length=num_bins+1))
    hist_outer = fit(Histogram, angles_outer, bin_edges)
    hist_inner = fit(Histogram, angles_inner, bin_edges)

    frequencies_outer = hist_outer.weights / maximum(hist_outer.weights)
    frequencies_inner = hist_inner.weights / maximum(hist_inner.weights)

    bar_thicknesses_outer = scale_factor .* frequencies_outer
    bar_thicknesses_inner = scale_factor .* frequencies_inner

    # Prepare circular plot
    θ = range(0, stop=2π, length=200)
    outer_circle_x = outer_radius .* cos.(θ)
    outer_circle_y = outer_radius .* sin.(θ)
    inner_circle_x = inner_radius .* cos.(θ)
    inner_circle_y = inner_radius .* sin.(θ)

    plt_hist_annular_circular = plot(outer_circle_x, outer_circle_y, lw=1.0, c=grey_color, aspect_ratio=:equal, label=false)
    plot!(plt_hist_annular_circular, inner_circle_x, inner_circle_y, lw=1.0, c=grey_color, label=false)

    # Function to add histogram bars in circular form
    function add_histogram_segments!(p, bin_edges, frequencies, base_radius, thicknesses, invert_direction=false)
        for i in 1:length(frequencies)
            r_base = base_radius
            r_mod = thicknesses[i]

            r_outer = r_base + r_mod
            r_inner = r_base

            if invert_direction
                r_outer = r_base
                r_inner = r_base - r_mod
            end

            θ_start, θ_end = bin_edges[i], bin_edges[i+1]
            npts = 10
            outer_arc = [(r_outer * cos(θ), r_outer * sin(θ)) for θ in range(θ_start, θ_end, length=npts)]
            inner_arc = [(r_inner * cos(θ), r_inner * sin(θ)) for θ in range(θ_end, θ_start, length=npts)]

            poly_points = vcat(outer_arc, inner_arc)
            xs, ys = first.(poly_points), last.(poly_points)

            max_brightness, min_brightness = 0.55, 0.25
            brightness = (max_brightness - min_brightness) * (1 - frequencies[i] / maximum(frequencies)) + min_brightness
            col = RGB(brightness*(153/255), 0, 0)

            plot!(p, xs, ys, seriestype=:shape, fillcolor=col, linecolor=:transparent, label=false)
        end
    end

    # Add histogram bars
    add_histogram_segments!(plt_hist_annular_circular, bin_edges, frequencies_outer, outer_radius, bar_thicknesses_outer, false)
    add_histogram_segments!(plt_hist_annular_circular, bin_edges, frequencies_inner, inner_radius, bar_thicknesses_inner, true)
    # Add start point
    scatter!(plt_hist_annular_circular, [start_point[1]], [start_point[2]], markersize=3, markercolor=:black, label="Start Point", marker=:dot)

    # Initialize 3D plot
    plt_hist_annular_3d = plot3d(outer_circle_x, outer_circle_y, zeros(length(outer_circle_x)), lw=0.5, c=:black, label=false, size=(800, 600))
    xlims!(-outer_radius-1.0, outer_radius+1.0)
    ylims!(-outer_radius-1.0, outer_radius+1.0)
    plot3d!(plt_hist_annular_3d, inner_circle_x, inner_circle_y, zeros(length(inner_circle_x)), lw=0.5, c=:black, label=false)

    # Function to add vertical bars in 3D
    function add_histogram_segments_3d!(p, bin_edges, frequencies, base_radius, heights)
        for i in 1:length(frequencies)
            r_base = base_radius
            z_height = heights[i]  # Bar height

            # Define angular span
            θ_start = bin_edges[i]
            θ_end = bin_edges[i+1]

            # Define base rectangle (bottom points)
            x_base = [r_base * cos(θ_start), r_base * cos(θ_end)]
            y_base = [r_base * sin(θ_start), r_base * sin(θ_end)]
            z_base = [0, 0]  # Always at z=0

            # Define top rectangle (same x-y but lifted to z_height)
            x_top = x_base
            y_top = y_base
            z_top = [z_height, z_height]

            # Compute grayscale intensity based on frequency
            max_brightness = 0.55
            min_brightness = 0.25
            brightness = (max_brightness - min_brightness) * (1 - frequencies[i] / maximum(frequencies)) + min_brightness
            col = RGB(brightness * (153/255), 0, 0)

            # Plot 4 vertical faces as rectangles
            for j in 1:2
                plot3d!(plt_hist_annular_3d, [x_base[j], x_top[j]], [y_base[j], y_top[j]], [z_base[j], z_top[j]],
                        lw=2, color=col, label=false)
            end
        end
    end

    # Add outer histogram bars
    add_histogram_segments_3d!(plt_hist_annular_3d, bin_edges, frequencies_outer, outer_radius, bar_thicknesses_outer)

    # Add inner histogram bars
    add_histogram_segments_3d!(plt_hist_annular_3d, bin_edges, frequencies_inner, inner_radius, bar_thicknesses_inner)

    # Add start point in 3D
    scatter3d!(plt_hist_annular_3d, [start_point[1]], [start_point[2]], [0], markersize=3, markercolor=:black, label=false, marker=:dot)

    display(plt_hist_annular_outer)
    display(plt_hist_annular_inner)
    display(plt_hist_annular_circular)
    display(plt_hist_annular_3d)
    return plt_hist_annular_outer, plt_hist_annular_inner, plt_hist_annular_circular, plt_hist_annular_3d
end

function plot_boundary_value(boundary_value)
    # Generate grid for the boundary of [-1,1] x [-1,1]
    v_values = range(-1.0f0, 1.0f0, length=100)  # For left and right edges
    x_values = range(-1.0f0, 1.0f0, length=100)  # For bottom and top edges

    # Prepare plot
    plot3d(xlabel="x", ylabel="v", zlabel="Boundary Value", title="3D Boundary Value on Square Boundary", legend=false, xlim = (-1,1), ylim=(-1,1), zlim = (0, 2))

    # Left and right edges (varying v, fixed x)
    plot3d!([-1.0f0 for v in v_values], v_values, [boundary_value([-1.0f0, v]) for v in v_values], lw=2, color=:blue)
    plot3d!([1.0f0 for v in v_values], v_values, [boundary_value([1.0f0, v]) for v in v_values], lw=2, color=:red)

    # Bottom and top edges (varying x, fixed v)
    plot3d!(x_values, [-1.0f0 for x in x_values], [boundary_value([x, -1.0f0]) for x in x_values], lw=2, color=:black)  # Bottom edge
    plot3d!(x_values, [1.0f0 for x in x_values], [boundary_value([x, 1.0f0]) for x in x_values], lw=2, color=:black)   # Top edge
end

function plot3d_square(x_values, v_values, grid_points, results, exit_rates)
    # Convert exit_rates dictionary to a matrix using broadcasting
    exit_rates_matrix = get.(Ref(exit_rates), grid_points, NaN)

    # Plot the heatmap
    plt_square_heatmap =  heatmap(x_values, v_values, exit_rates_matrix, xlabel="X", ylabel="V", color=:viridis, title="Exit Rates Heatmap", colorbar=true)

    results_matrix = get.(Ref(results), grid_points, NaN)

    
    # Find the center and range of each axis
    x_center, v_center, z_center = (mean(x_values), mean(v_values), mean(results_matrix))
    max_range = maximum([maximum(x_values) - minimum(x_values), maximum(v_values) - minimum(v_values), maximum(results_matrix) - minimum(results_matrix)]) / 2

    # Set equal axis limits around the center
    xlims = (x_center - max_range, x_center + max_range)
    ylims = (v_center - max_range, v_center + max_range)
    zlims = (z_center - max_range, z_center + max_range)

    plt_square_3d = plot3d(x_values,  # x-coordinates for surface
        v_values,  # y-coordinates for surface
        results_matrix,  # Transposed result matrix for surface
        st = :surface,  # Surface plot style
        xlabel = "x",
        ylabel = "v",
        zlabel = "u(x,v)",
        color = :viridis,  # Use a color gradient
        legend = false,  # Disable the legend
        xlims=xlims, ylims=ylims, zlims=zlims
        )
    
    
    display(plt_square_heatmap)
    display(plt_square_3d)

    return plt_square_heatmap, plt_square_3d
end

end  # module