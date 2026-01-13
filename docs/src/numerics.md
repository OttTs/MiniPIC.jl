# Numerical Methods

## Particle Discretization

The distribution function is approximated by ``N_p`` simulation particles:

```math
f(x,v,t) \approx \sum_{p=1}^{N_p} \omega \cdot S(x - x_p) \cdot \delta(v - v_p)
```

Each simulation particle has:
- Position ``x_p``
- Velocity ``v_p``
- Weighting ``\omega = N/N_p`` (number of physical particles per simulation particle)

The shape function ``S`` is introduced together with the discretization of the Poisson equation below.

The trajectory of each particle follows Newton's equations:

```math
\begin{aligned}
\frac{\mathrm{d}x_p}{\mathrm{d}t} &= v_p \\
\frac{\mathrm{d}v_p}{\mathrm{d}t} &= \frac{q}{m} E(x_p)
\end{aligned}
```

## Time Integration

The **Leapfrog method** is used to integrate the particle trajectories:

```math
\begin{aligned}
x_p^{n+1} &= x_p^n + \Delta t \, v_p^{n+1/2} \\
v_p^{n+1/2} &= v_p^{n-1/2} + \frac{q}{m} E^n(x_p^n)
\end{aligned}
```

## Poisson Solver

The computational domain is discretized by ``N_x`` equally spaced grid points with step size ``\Delta x``.

The Laplacian is approximated by a central difference quotient:

```math
\frac{\varphi_{i+1} - 2\varphi_{i} + \varphi_{i-1}}{\Delta x^2} = -\frac{1}{\varepsilon} \rho_i
```

The electric field is computed on a staggered grid:

```math
E_{i+1/2} = -\frac{\varphi_{i+1} - \varphi_{i}}{\Delta x}
```

### Charge Deposition

The charge density on the grid is obtained by depositing particles using the shape function:

```math
\rho_i = q \cdot \omega \sum_{p=1}^{N_p} S(x_i - x_p) - \frac{q \cdot \omega \cdot N_p}{L}
```

The **linear shape function** (first-order b-spline) is:

```math
S(r) = \begin{cases}
(\Delta x + r) / \Delta x^2 & -\Delta x < r < 0 \\
(\Delta x - r) / \Delta x^2 & 0 \leq r < \Delta x \\
0 & \text{otherwise}
\end{cases}
```

### Gauge Condition

The Poisson equation with periodic boundary conditions is singular (the potential is defined up to a constant). We fix the gauge by setting ``\varphi_1 = 0``.

### Field Interpolation

The electric field at particle positions is obtained by linear interpolation:

```math
E(x_p) = \sum_{i=1}^{N_x} E(x_{i+1/2}) \, S(x_{i+1/2} - x_p)
```

## Non-dimensionalization

To improve numerical efficiency, the following non-dimensional variables are introduced:

```math
\bar{x} = \frac{x}{\Delta x}, \quad
\bar{v} = \frac{\Delta t \, v}{\Delta x}, \quad
\bar{\varphi} = -\frac{\varepsilon \, \varphi}{\omega \, q \, \Delta x}, \quad
\bar{E} = \frac{q}{m} \frac{\Delta t^2}{\Delta x} E
```

### Particle Equations

In non-dimensional form, the particle equations become:

```math
\begin{aligned}
\bar{x}_p^{n+1} &= \bar{x}_p^n + \bar{v}_p^{n+1/2} \\
\bar{v}_p^{n+1/2} &= \bar{v}_p^{n-1/2} + \bar{E}^n(\bar{x}_p^n)
\end{aligned}
```

This makes particle updates as efficient as possible (no multiplications required).

### Field Equations

The finite difference equations reduce to:

```math
\begin{aligned}
\bar{\varphi}_{i+1} - 2\bar{\varphi}_{i} + \bar{\varphi}_{i-1} &= \sum_{p=1}^{N_p} \bar{S}(\bar{x}_i - \bar{x}_p) - \frac{N_p}{N_x} \\
\bar{E}_{i+1/2} &= \alpha (\bar{\varphi}_{i+1} - \bar{\varphi}_{i})
\end{aligned}
```

where the single remaining parameter is:

```math
\alpha = \frac{\omega \, q^2 \, \Delta t^2}{m \, \varepsilon \, \Delta x}
```

This parameter characterizes the strength of the electric field on the particles.

### Shape Function

The non-dimensional shape function is:

```math
\bar{S}(r) = \begin{cases}
1 + r & -1 < r < 0 \\
1 - r & 0 \leq r < 1 \\
0 & \text{otherwise}
\end{cases}
```