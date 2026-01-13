# MiniPIC.jl

A minimal 1D electrostatic Particle-in-Cell (PIC) simulation package for Julia.

## Overview

MiniPIC.jl provides a simple yet efficient implementation of the particle-in-cell method for simulating collisionless plasmas. It solves the Vlasov-Poisson system using:

- **Leapfrog integration** for particle trajectories
- **Second-order finite differences** for the Poisson equation
- **Linear shape functions** for charge deposition and field interpolation
- **Periodic boundary conditions**

The package supports two PIC schemes:
- **MC-PIC** (Momentum-Conserving): Standard PIC with centered electric field nodes
- **EC-PIC** (Energy-Conserving): PIC with staggered electric field nodes

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/OttTs/MiniPIC.jl")
```

## Quick Start

```julia
using MiniPIC

# Set up a simple Landau damping simulation
Nₚ = 10000        # Number of particles
L = 2π            # Domain length

# Initialize particles: uniform positions with sinusoidal perturbation
x = rand(Nₚ) * L
v = randn(Nₚ)     # Maxwellian velocity distribution

# Define simulation parameters
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
simulate!(x, v; parameters=p, file_name="landau_damping")

# Read results
using CSV, DataFrames
df = CSV.read("landau_damping.csv", DataFrame)
```

## Contents

```@contents
Pages = ["physical-model.md", "numerics.md", "api.md"]
Depth = 2
```