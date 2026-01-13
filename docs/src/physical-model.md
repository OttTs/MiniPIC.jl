# Physical Model

## Vlasov-Poisson System

MiniPIC.jl simulates collisionless plasmas by solving the Vlasov-Poisson system in one spatial dimension.

### Vlasov Equation

The evolution of the distribution function ``f(x, v, t)`` of a charged species is governed by the Vlasov equation:

```math
\frac{\partial f}{\partial t} + v \frac{\partial f}{\partial x} + \frac{q}{m} E \frac{\partial f}{\partial v} = 0
```

where:
- ``q`` is the particle charge
- ``m`` is the particle mass
- ``E`` is the electric field

### Poisson Equation

The electric field is determined self-consistently through the electrostatic Poisson equation:

```math
\frac{\partial^2 \varphi}{\partial x^2} = -\frac{\rho}{\varepsilon}
```

```math
E = -\frac{\partial \varphi}{\partial x}
```

where ``\varphi`` is the electric potential and ``\varepsilon`` is the permittivity.

### Charge Density

The charge density is obtained by integrating the distribution function over velocity space:

```math
\rho = q \int f \, \mathrm{d}v - Q
```

where ``Q`` is a uniform neutralizing background charge:

```math
Q = \frac{q}{L} \iint f \, \mathrm{d}x \, \mathrm{d}v
```

This background ensures quasi-neutrality of the plasma.

## Boundary Conditions

The simulation uses **periodic boundary conditions** in both position and field quantities. Particles leaving the domain on one side re-enter from the opposite side.