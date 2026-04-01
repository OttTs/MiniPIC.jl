module MiniPIC

using Polyester: @batch
using DataFrames: DataFrame
using CSV, JLD2
using ProgressMeter
using Statistics: mean, var

"""
    BOLTZMANN_CONSTANT

Boltzmann constant in SI units (J/K).

``k_B = 1.380649 \\times 10^{-23}`` J/K
"""
const BOLTZMANN_CONSTANT = 1.380649e-23

"""
    VACUUM_PERMITTIVITY

Vacuum permittivity (electric constant) in SI units (F/m).

``\\varepsilon_0 = 8.8541878128 \\times 10^{-12}`` F/m
"""
const VACUUM_PERMITTIVITY = 8.8541878128E-12

"""
    ELEMENTARY_CHARGE

Elementary charge in SI units (C).

``e = 1.60217663 \\times 10^{-19}`` C
"""
const ELEMENTARY_CHARGE = 1.60217663E-19

"""
    ELECTRON_MASS

Electron mass in SI units (kg).

``m_e = 9.1093837 \\times 10^{-31}`` kg
"""
const ELECTRON_MASS = 9.1093837E-31

"""
    Parameters(; mesh_length=1, num_cells, time_step, num_steps,
                 charge=-ELEMENTARY_CHARGE, mass=ELECTRON_MASS,
                 weighting=1, permittivity=VACUUM_PERMITTIVITY,
                 do_MC_PIC=true)

A struct to hold simulation parameters.

# Arguments
- `mesh_length::Number=1`: Length of the computational domain
- `num_cells::Integer`: Number of cells used to divide the domain
- `time_step::Number`: Simulation time step
- `num_steps::Integer`: Total number of time steps in the simulation
- `charge::Number=-ELEMENTARY_CHARGE`: Species charge
- `mass::Number=ELECTRON_MASS`: Species mass
- `weighting::Number=1`: Particle weighting, number of physical particles per simulation particle
- `permittivity::Number=VACUUM_PERMITTIVITY`: Absolute permittivity used to solve the Poisson equation
- `do_MC_PIC::Bool=true`: If true, use (standard) MC-PIC scheme, else use EC-PIC scheme
"""
struct Parameters
    mesh_length::Float64
    num_cells::Int
    time_step::Float64
    num_steps::Int
    charge::Float64
    mass::Float64
    weighting::Float64
    permittivity::Float64
    do_MC_PIC::Bool
    function Parameters(;
        mesh_length::Number=1,
        num_cells::Integer,
        time_step::Number,
        num_steps::Integer,
        charge::Number=-ELEMENTARY_CHARGE,
        mass::Number=ELECTRON_MASS,
        weighting::Number=1,
        permittivity::Number=VACUUM_PERMITTIVITY,
        do_MC_PIC::Bool=true)
        return new(
            Float64(mesh_length), Int(num_cells), Float64(time_step), Int(num_steps),
            Float64(charge), Float64(mass), Float64(weighting), Float64(permittivity),
            do_MC_PIC)
    end
end

include("fieldsolver.jl")
include("output.jl")
include("simulation.jl")

export BOLTZMANN_CONSTANT, VACUUM_PERMITTIVITY
export ELEMENTARY_CHARGE, ELECTRON_MASS
export Parameters, simulate!
export read_particles, read_parameters

end
