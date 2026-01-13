"""
    simulate!(x, v; parameters, file_name, write_interval=1, write_particles=false, show_progress=true)

Run a 1D electrostatic particle-in-cell simulation with periodic boundary conditions.

This function evolves the particle positions `x` and velocities `v` in-place using the
leapfrog integration scheme. The electric field is computed self-consistently by solving
the Poisson equation with a second-order finite difference method.

A uniform neutralizing background charge is automatically included for quasi-neutrality.

# Arguments
- `x::AbstractVector`: Initial particle positions (modified in-place, non-dimensionalized during simulation)
- `v::AbstractVector`: Initial particle velocities (modified in-place, non-dimensionalized during simulation)

# Keyword Arguments
- `parameters::Parameters`: Simulation parameters (required)
- `file_name::String`: Base name for output files (required). Creates `file_name.csv` and `file_name.jld2`
- `write_interval::Integer=1`: Number of time steps between diagnostic outputs
- `write_particles::Bool=false`: If `true`, write particle data at each output interval;
  particle data is always written at the first and last time step
- `show_progress::Bool=true`: If `true`, display a progress bar

# Output Files
- `file_name.csv`: Time series of kinetic and potential energy
- `file_name.jld2`: HDF5-compatible file containing simulation parameters and particle data

# Example
```julia
using MiniPIC

# Initialize particles
Nₚ = 10000
x = rand(Nₚ) * 2π  # Uniform in [0, 2π]
v = randn(Nₚ)       # Maxwellian velocity distribution

# Set up parameters
p = Parameters(
    mesh_length = 2π,
    num_cells = 64,
    time_step = 0.1,
    num_steps = 100
)

# Run simulation
simulate!(x, v; parameters=p, file_name="my_simulation")
```

See also: [`Parameters`](@ref), [`read_particles`](@ref), [`read_parameters`](@ref)
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

    # Leapfrog, half step for first push
    #fieldsolver!(x, E, α, do_MC_PIC) # Already done before calling this function
    @batch for i in eachindex(v)
        v[i] += 0.5 * interpolate(x[i], E, do_MC_PIC)
        x[i] = mod(x[i] + v[i], Nₓ)
    end

    # Main loop, full steps
    for _ in 2:num_steps
        fieldsolver!(x, E, α, do_MC_PIC)
        @batch for i in eachindex(v)
            v[i] += interpolate(x[i], E, do_MC_PIC)
            x[i] = mod(x[i] + v[i], Nₓ)
        end

        isnothing(progress) || next!(progress)
    end

    # Leapfrog, half step for velocity (forwards, to sync with position)
    fieldsolver!(x, E, α, do_MC_PIC)
    @batch for i in eachindex(v)
        v[i] += 0.5 * interpolate(x[i], E, do_MC_PIC)
    end
end
