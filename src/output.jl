function initialize_output(file_name::String, parameters::Parameters)
    jldopen(file_name * ".jld2", "w") do file
        file["parameters"] = parameters
    end

    return DataFrame(time=Float64[], Ekin=Float64[], Epot=Float64[])
end


function write_state!(x, v, E, iter; df, parameters, file_name, write_particles)
    Δx = parameters.mesh_length / parameters.num_cells

    x = Δx * x
    v = (Δx / parameters.time_step) * v

    push!(df, (time=iter * parameters.time_step,
                Ekin=kinetic_energy(v, parameters.mass),
                Epot=potential_energy(E, Δx, parameters.permittivity)))

    if write_particles
        jldopen(file_name * ".jld2", "a+") do file
            file["x_$(iter)"] = x
            file["v_$(iter)"] = v
        end
    end
end


kinetic_energy(v, m) = 0.5 * m * sum(vᵢ^2 for vᵢ in v)
potential_energy(E, Δx, ε) = 0.5 * ε * sum(Eᵢ^2 for Eᵢ in E) * Δx


"""
    read_particles(file_name::String, iter::Integer) -> (x, v)

Read particle positions and velocities from a simulation output file.

# Arguments
- `file_name::String`: Base name of the output file (without `.jld2` extension)
- `iter::Integer`: Time step index to read (0 for initial conditions)

# Returns
- `x::Vector{Float64}`: Particle positions at the specified time step
- `v::Vector{Float64}`: Particle velocities at the specified time step

# Example
```julia
# Read initial conditions
x0, v0 = read_particles("my_simulation", 0)

# Read final state (assuming 100 time steps)
x_final, v_final = read_particles("my_simulation", 100)
```

See also: [`simulate!`](@ref), [`read_parameters`](@ref)
"""
function read_particles(file_name::String, iter::Integer)
    x, v = jldopen(file_name * ".jld2", "r") do file
        file["x_$(iter)"], file["v_$(iter)"]
    end
    return x, v
end


"""
    read_parameters(file_name::String) -> Parameters

Read simulation parameters from an output file.

# Arguments
- `file_name::String`: Base name of the output file (without `.jld2` extension)

# Returns
- `parameters::Parameters`: The simulation parameters used in the run

# Example
```julia
# Read parameters from a previous simulation
p = read_parameters("my_simulation")
println("Time step: ", p.time_step)
println("Number of cells: ", p.num_cells)
```

See also: [`simulate!`](@ref), [`read_particles`](@ref), [`Parameters`](@ref)
"""
function read_parameters(file_name::String)
    parameters = jldopen(file_name * ".jld2", "r") do file
        file["parameters"]
    end
    return parameters
end
