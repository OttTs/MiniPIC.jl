"""
    simulate(x::AbstractVector{T}, v::AbstractVector{T}; settings...) where {T}

Run a one-dimensional electro-static particle-in-cell simulation with periodic boundary conditions.
The initial conditions are given by the particle's position `x` and velocity `v`.
A uniform background charge is calculated for quasi-neutrality.

The electric field is calculated with a second order finite difference scheme.
The particle positions `x` and velocities `v` are updated using the Leapfrog method.

# Settings
- `parameters::Parameters`: Simulation parameters
- `write_interval::Integer=1`: Number of iterations between two outputs
- `show_progress::Bool=true`: Flag whether to show the progress bar

Returns a DataFrame with the solution at the requested times.
"""
function simulate!(x::AbstractVector{T}, v::AbstractVector{T};
    parameters::Parameters,
    write_interval::Integer=1,
    write_particles::Bool=false,
    file_name::String,
    show_progress::Bool=true) where {T}

    progress = show_progress ? Progress(parameters.num_steps) : nothing

    # Prepare simulation (non-dimensionalization)
    Δx = parameters.mesh_length / parameters.num_cells
    x /= Δx
    v *= parameters.time_step / Δx
    α = parameters.time_step^2 * parameters.weighting * parameters.charge^2 /
        (Δx * parameters.mass * parameters.permittivity)
    E = zero(similar(x, parameters.num_cells))
    fieldsolver!(x, E, α, parameters.do_MC_PIC)

    df = initialize_output(file_name, parameters)
    write_state!(x, v, E, 0; df, parameters, file_name, write_particles=true)

    iter = 0
    while iter < parameters.num_steps
        # Perform until next output or end of simulation
        Nₜ = min(write_interval, parameters.num_steps - iter)
        simulate_steps!(x, v, E, α; num_steps = Nₜ, parameters.do_MC_PIC, progress)
        iter += Nₜ

        write_state!(x, v, E, iter;
            df,
            parameters,
            file_name,
            write_particles=(write_particles || iter==parameters.num_steps))
    end

    CSV.write(file_name * ".csv", df)

    show_progress && finish!(progress)
end

function simulate_steps!(x, v, E, α; num_steps, do_MC_PIC, progress)
    Nₓ = length(E)

    # Leapfrog, half step for velocity (backwards)
    @batch for i in eachindex(v)
        v[i] -= 0.5 * interpolate(x[i], E, do_MC_PIC)
    end

    # Main loop, full steps
    for _ in 1:num_steps
        @batch for i in eachindex(v)
            v[i] += interpolate(x[i], E, do_MC_PIC)
            x[i] = mod(x[i] + v[i], Nₓ)
        end
        fieldsolver!(x, E, α, do_MC_PIC)

        isnothing(progress) || next!(progress)
    end

    # Leapfrog, half step for velocity (forwards, to sync with position)
    @batch for i in eachindex(v)
        v[i] += 0.5 * interpolate(x[i], E, do_MC_PIC)
    end
end
