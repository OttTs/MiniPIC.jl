"""
    simulate(x::AbstractVector{T}, v::AbstractVector{T}; settings...) where {T}

Run a one-dimensional electro-static particle-in-cell simulation with periodic boundary conditions.
The initial conditions are given by the particle's position `x` and velocity `v`.
A uniform background charge is calculated for quasi-neutrality.

The electric field is calculated with a second order finite difference scheme.
The particle positions `x` and velocities `v` are updated using the Leapfrog method.

# Settings
- `mesh_length::Number=1`: Length of the computational domain
- `num_cells::Integer`: Number of cells used to divide the domain
- `time_step::Number`: Simulation time step
- `num_steps::Integer`: Total number of time steps in the simulation
- `charge::Number=-ELEMENTARY_CHARGE`: Species charge
- `mass::Number=ELECTRON_MASS`: Species mass
- `weighting::Number=1`: Particle weighting, number of physical particles per simulation particle
- `permittivity::Number=VACUUM_PERMITTIVITY`: Absolute permittivity used to solve the Poisson equation
- `write_interval::Integer=num_steps`: Number of iterations between two outputs
- `write_particles::Bool=true`: Flag whether to write out the particle state
- `write_field::Bool=true`: Flag whether to write out the electric field
- `show_progress::Bool=true`: Flag whether to show the progress bar
- `do_MC_PIC::Bool=true`: If true, use (standard) MC-PIC scheme, else use EC-PIC scheme

Returns a DataFrame with the solution at the requested times.
"""
function simulate(x::AbstractVector{T}, v::AbstractVector{T};
        mesh_length::Number=1,
        num_cells::Integer,
        time_step::Number,
        num_steps::Integer,
        charge::Number=-ELEMENTARY_CHARGE,
        mass::Number=ELECTRON_MASS,
        weighting::Number=1,
        permittivity::Number=VACUUM_PERMITTIVITY,
        write_interval::Integer=num_steps,
        write_particles::Bool=true,
        write_field::Bool=false,
        show_progress::Bool=true,
        do_MC_PIC::Bool=true) where {T}

    progress = show_progress ? Progress(num_steps) : nothing

    # Initialize output
    output = (time=T[], Ekin=T[], Epot=T[])
    if write_particles; output=(output..., x=typeof(x)[], v=typeof(v)[]); end
    if write_field; output=(output..., E=typeof(x)[]); end
    output = DataFrame(output)

    # Prepare simulation (non-dimensionalization)
    Δx = mesh_length / num_cells
    x /= Δx
    v *= time_step / Δx
    α = time_step^2 * weighting * charge^2 / (Δx * mass * permittivity)
    E = zero(similar(x, num_cells))

    iter = 0
    while iter < num_steps
        # Perform until next output or end of simulation
        Nₜ = min(write_interval, num_steps - iter)
        simulate_steps!(x, v, E, α; num_steps = Nₜ, do_MC_PIC, progress)
        iter += Nₜ

        # Write output
        xₒ = Δx * x
        vₒ = (Δx / time_step) * v
        Eₒ = (mass / charge) * (Δx / time_step^2) * E
        data = (
            time_step * iter,
            kinetic_energy(vₒ, weighting * mass),
            potential_energy(Eₒ, Δx, permittivity))
        if write_particles; data=(data..., xₒ, vₒ); end
        if write_field; data=(data..., Eₒ); end
        push!(output, data)
    end

    show_progress && finish!(progress)

    return output
end

function simulate_steps!(x, v, E, α; num_steps, do_MC_PIC, progress)
    Nₓ = length(E)

    # Leapfrog, half step for velocity (backwards)
    @. v -= 0.5 * interpolate(x, $Ref(E), do_MC_PIC)

    # Main loop, full steps
    for _ in 1:num_steps
        @. v += interpolate(x, $Ref(E), do_MC_PIC)
        @. x = mod(x + v, Nₓ)
        fieldsolver!(x, E, α, do_MC_PIC)

        isnothing(progress) || next!(progress)
    end

    # Leapfrog, half step for velocity (forwards, to sync with position)
    @. v += 0.5 * interpolate(x, $Ref(E), do_MC_PIC)
end

kinetic_energy(v, m) = 0.5 * m * sum(vᵢ^2 for vᵢ in v)
potential_energy(E, Δx, ε) = 0.5 * ε * sum(Eᵢ^2 for Eᵢ in E) * Δx
