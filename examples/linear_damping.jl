# f(x, v, 0) = 1 / (√(2π) * L) * (1 + α * cos(k * x)) * exp(-v^2 / 2)
# linear landau damping: α = 0.05, nonlinear: α = 0.5

using MiniPIC
using GLMakie
using Roots

Nₚ = 500000
α = 0.05
k = 0.5
L = 2π / k

x = [find_zero(x -> y - (x + α / k * sin(k * x)) / L, L*y) for y in range(0, 1, length=Nₚ)]
v = randn(Nₚ)

output = simulate(x, v;
    mesh_length=L,
    num_cells=1000,
    time_step=0.02,
    num_steps=1000,
    charge=-1,
    mass=1,
    weighting=L/Nₚ,
    permittivity=1,
    write_interval=1,
    write_particles=false,
    do_MC_PIC=false
)


fig = Figure()
ax = Axis(fig[1, 1], yscale=log10)
lines!(ax, output.time, output.Epot)

fig