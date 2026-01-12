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


function read_particles(file_name::String, iter::Integer)
    x, v = jldopen(file_name * ".jld2", "r") do file
        file["x_$(iter)"], file["v_$(iter)"]
    end
    return x, v
end


function read_parameters(file_name::String)
    parameters = jldopen(file_name * ".jld2", "r") do file
        file["parameters"]
    end
    return parameters
end
