# MiniPIC.jl

[![Build Status](https://github.com/OttTs/MiniPIC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/OttTs/MiniPIC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://OttTs.github.io/MiniPIC.jl/dev)

A minimal 1D electrostatic Particle-in-Cell (PIC) simulation package for Julia.

## Features

- **Vlasov-Poisson solver** for collisionless plasma simulations
- **Leapfrog integration** for symplectic time-stepping
- **Linear shape functions** for charge deposition and field interpolation
- **Periodic boundary conditions**
- **Two PIC schemes**: MC-PIC (momentum-conserving) and EC-PIC (energy-conserving)
- **Multi-threaded** particle operations using [Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/OttTs/MiniPIC.jl")
```

## Quick Start

```julia
using MiniPIC

# Initialize particles
Nₚ = 10000
L = 2π
x = rand(Nₚ) * L      # Uniform positions
v = randn(Nₚ)         # Maxwellian velocities

# Set up simulation parameters
p = Parameters(
    mesh_length = L,
    num_cells = 64,
    time_step = 0.1,
    num_steps = 200,
    charge = -1.0,
    mass = 1.0,
    weighting = L / Nₚ,
    permittivity = 1.0
)

# Run simulation
simulate!(x, v; parameters=p, file_name="simulation")
```

## Output

The simulation produces two output files:
- `simulation.csv` — Time series of kinetic and potential energy
- `simulation.jld2` — Particle data and parameters (HDF5-compatible)

Read results with:
```julia
using CSV, DataFrames
df = CSV.read("simulation.csv", DataFrame)

# Read particle data at specific simulation step
x_final, v_final = read_particles("simulation", 200)
```

## Documentation

For more details on the physical model and numerical methods, see the [documentation](https://OttTs.github.io/MiniPIC.jl/dev).

## License

MIT License