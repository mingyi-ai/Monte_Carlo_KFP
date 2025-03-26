# MCKernel.jl
module MCKernel

# Custom package for GPU friendly random number generator
using Metal
include("RandomMetal.jl")
using .RandomMetal

export MCkernel_square, MCkernel_annulus, MCkernel_doublewell

"""
    simulate_particle_movement(start_points, end_points, exit_flags, time_step, num_steps, points_before_end, seeds)

Simulates particle movement, excutes exit check in a square domain.

# Arguments
- `start_points::MtlArray{Float32, 2}`: Initial positions and velocities of particles.
- `end_points::MtlArray{Float32, 2}`: Final positions and velocities of particles.
- `exit_flags::MtlArray{Bool, 1}`: Flags indicating if particles have exited the domain.
- `time_step::Float32`: Time step for the simulation.
- `num_steps::Int32`: Number of steps in the simulation.
- `points_before_end::MtlArray{Float32, 2}`: Positions and velocities of particles before they exit.
- `seeds::MtlArray{UInt32, 2}`: Random seeds for generating random numbers.

# Returns
- `Nothing`: This function modifies the input arrays in place.
"""
function MCkernel_square(start_points, end_points, exit_flags, time_step, num_steps, points_before_end, seeds)
    i = thread_position_in_grid_1d() # Thread ID

    x = start_points[1, i]
    v = start_points[2, i]

    x_prev, v_prev = 0.0f0, 0.0f0  # Initialize the point before exit

    seed1 = seeds[1, i]
    seed2 = seeds[2, i]

    for step in 1:num_steps
        # Boolean masks for exit conditions
        mask_x = (x < -1.0f0 || x > 1.0f0) ? 1 : 0
        mask_v = (v < -1.0f0 || v > 1.0f0) ? 1 : 0
        mask_exit = mask_x | mask_v  # Combine masks (exit if either condition is true)
        continue_mask = 1 - mask_exit  # 1 = active, 0 = exited

        # Generate two uniform distributed random numbers
        seed1 = xorshift32(seed1)
        seed2 = xorshift32(seed2)
        random_number1 = xorshift32_float(seed1)
        random_number2 = xorshift32_float(seed2)

        # Generate a normal distributed noise
        noise = box_muller(random_number1, random_number2)

        # Perturb the seeds to avoid deterministic patterns
        seed1 += UInt32(i)
        seed2 += UInt32(i)

        # Update position and velocity and store previous state if not exit
        x_prev, v_prev = continue_mask * x + mask_exit * x_prev, continue_mask * v + mask_exit * v_prev 
        x += continue_mask * (v * time_step)
        v += continue_mask * (sqrt(time_step) * noise)
    end
    # Exit check
    mask_x = (x < -1.0f0 || x > 1.0f0) ? 1 : 0
    mask_v = (v < -1.0f0 || v > 1.0f0) ? 1 : 0
    mask_exit = mask_x | mask_v  # Combine masks (exit if either condition is true)

    # Store results
    end_points[1, i] = x
    end_points[2, i] = v
    points_before_end[1, i] = x_prev
    points_before_end[2, i] = v_prev
    exit_flags[i] = mask_exit

    return
end


"""
    simulate_particle_movement_annulus(start_points, end_points, exit_flags, time_step, num_steps, seeds; inner_radius2=2.0f0^2, outer_radius2=4.0f0^2)

Simulates particle movement, executes exit check in an annulus domain.

# Arguments
- `start_points::MtlArray{Float32, 2}`: Initial positions and velocities of particles.
- `end_points::MtlArray{Float32, 2}`: Final positions and velocities of particles.
- `exit_flags::MtlArray{Bool, 2}`: Flags indicating if particles have exited the inner or outer domain.
- `time_step::Float32`: Time step for the simulation.
- `num_steps::Int32`: Number of steps in the simulation.
- `seeds::MtlArray{UInt32, 2}`: Random seeds for generating random numbers.
- `inner_radius2::Float32`: Square of the inner radius of the annulus.
- `outer_radius2::Float32`: Square of the outer radius of the annulus.

# Returns
- `Nothing`: This function modifies the input arrays in place.
"""
function MCkernel_annulus(start_points, end_points, exit_flags, time_step, num_steps, seeds, inner_radius2, outer_radius2)
    i = thread_position_in_grid_1d() # Thread ID

    seed1 = seeds[1, i]
    seed2 = seeds[2, i]

    x = start_points[1, i]
    v = start_points[2, i]

    # Fixed step iteration
    for step in 1:num_steps
        radius2 = x^2 + v^2

        # Mask to avoid branching
        mask_inner = radius2 <= inner_radius2
        mask_outer = radius2 >= outer_radius2
        continue_mask = (!mask_inner & !mask_outer)

        # Generate two uniform distributed random numbers
        seed1 = xorshift32(seed1)
        seed2 = xorshift32(seed2)
        random_number1 = xorshift32_float(seed1)
        random_number2 = xorshift32_float(seed2)

        # Generate a normal distributed noise
        noise = box_muller(random_number1, random_number2)

        # Perturb the seeds to avoid deterministic patterns
        seed1 += UInt32(i)
        seed2 += UInt32(i)

        # Update position and velocity (applies only if continue_mask is true)
        x += continue_mask * (v * time_step)
        v += continue_mask * (sqrt(time_step) * noise)
    end
    # Exit check
    radius2 = x^2 + v^2
    mask_inner = radius2 <= inner_radius2
    mask_outer = radius2 >= outer_radius2

    # Store results
    exit_flags[1, i] = mask_inner
    exit_flags[2, i] = mask_outer
    
    end_points[1, i] = x
    end_points[2, i] = v

    return
end

function MCkernel_doublewell(start_points, end_points, exit_flags, time_step, num_steps, points_before_end, seeds, alpha, grad, left_center_x, left_radius2, right_center_x, right_radius2)
    i = thread_position_in_grid_1d() # Thread ID

    seed1 = seeds[1, i]
    seed2 = seeds[2, i]

    x = start_points[1, i]
    v = start_points[2, i]

    x_prev, v_prev = 0.0f0, 0.0f0  # Initialize the point before exit

    # Fixed step iteration
    for step in 1:num_steps
        left_dist2 = (x - left_center_x)^2 + v^2
        right_dist2 = (x - right_center_x)^2 + v^2

        # Mask to avoid branching
        mask_left = left_dist2 <= left_radius2
        mask_right = right_dist2 <= right_radius2
        mask_exit = mask_left | mask_right
        continue_mask = !mask_exit

        # Generate two uniform distributed random numbers
        seed1 = xorshift32(seed1)
        seed2 = xorshift32(seed2)
        random_number1 = xorshift32_float(seed1)
        random_number2 = xorshift32_float(seed2)

        # Generate a normal distributed noise
        noise = box_muller(random_number1, random_number2)

        # Perturb the seeds to avoid deterministic patterns
        seed1 += UInt32(i)
        seed2 += UInt32(i)

        # Update position and velocity (applies only if continue_mask is true)
        x_prev, v_prev = continue_mask * x + mask_exit * x_prev, continue_mask * v + mask_exit * v_prev 
        grad_x = grad(x)
        x += continue_mask * (v * time_step)
        v += continue_mask * (-(1.0f0 - alpha) * v - grad_x * time_step +  sqrt(time_step) * noise)
    end
    # Exit check
    left_dist2 = (x - left_center_x)^2 + v^2
    right_dist2 = (x - right_center_x)^2 + v^2

    mask_left = left_dist2 <= left_radius2
    mask_right = right_dist2 <= right_radius2

    # Store results
    exit_flags[1, i] = mask_left
    exit_flags[2, i] = mask_right
    
    end_points[1, i] = x
    end_points[2, i] = v

    points_before_end[1, i] = x_prev
    points_before_end[2, i] = v_prev    

    return
end

end  # module MCKernel