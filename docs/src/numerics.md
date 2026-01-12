## Numerics
We first approximate the density function by $N_p$ simulation particles
$$
f(x,v,t)\approx\sum_p^{N_p}{\omega\cdot S(x - x_p)\cdot\delta(v - v_p)}.
$$

Next to its position $x_p$ and velocity $v_p$, each simulation particle has a weighting $\omega=N/N_p$, representing the number of physical particles per simulation particle. In order to keep things simple, each particle has the same weighting. We introduce the shape function $S$ together with the discretization of the Poisson equation below.

With this, the trajectory of each particle is obtained by
$$
\begin{aligned}
    \frac{\mathrm{d}x_p}{\mathrm{d}t} &= v_p\\
    \frac{\mathrm{d}v_p}{\mathrm{d}t} &= \frac{q}{m}E(x_p)
\end{aligned}
$$

### Time integration
We use the Leapfrog method to integrate the particle trajectories, leading to
$$
\begin{aligned}
    x_p^{n+1} &= x_p^n + \Delta t\,v_p^{n+1/2}\\
    v_p^{n+1/2} &= v_p^{n-1/2} + \frac{q}{m}E^n(x_p^n)
\end{aligned}
$$

### Poisson solver
First, the computational domain is discretized by $N_x$ equally spaced grid points with step size $\Delta x$.
The second order derivative in the Poisson equation is then approximated by a central difference quotient, leading to
$$
\frac{\varphi_{i+1} - 2\varphi_{i} + \varphi_{i-1}}{\Delta x^2} = -\frac{1}{\varepsilon}\rho_i.
$$

The electric field is afterwards obtained on a staggered grid also using a central difference quotient for the first derivative
$$
E_{i+1/2} = - \frac{\varphi_{i+1} - \varphi_{i}}{\Delta x}
$$

For the charge density, we can insert the approximation for $f$ into our definition of the charge density. This leads to
$$
    \rho_i = q\cdot\omega\sum_p^{N_p}S(x_i-x_p) - \frac{q\cdot\omega\cdot N_p}{L}.
$$
The shape function is given by
$$
S(r) = \begin{cases}
(\Delta x + r) / \Delta x^2 & -\Delta x < r < 0 \\
(\Delta x - r) / \Delta x^2 & 0 < r < \Delta x
\end{cases}.
$$


There is one remaining issue: The resulting system is singular. This is because the reference point for the zero electric potential can be chosen arbitrarily. We can simply set $\varphi_1=0$.

The electric field at the position of a particle is given by linear interpolation
$$
E(x_p) = \sum_{i}^{N_x} E(x_{i+1/2})S(x_{i+1/2} - x_p)
$$

## Nondimensionalization
We can reduce the number of parameters by introducing the following nondimensionalized variables:
$$
    \bar{x}=\frac{x}{\Delta x},\quad
    \bar{v}=\frac{\Delta t\,v}{\Delta x},\quad
    \bar{\varphi}=-\frac{\varepsilon\,\varphi}{\omega\,q\Delta x},\quad
    \bar{E} = \frac{q}{m}\frac{\Delta t^2}{\Delta x}E
$$

With these new variables, we have for the particle paths
$$
\begin{aligned}
    \bar{x}_p^{n+1} &= \bar{x}_p^n + \bar{v}_p^{n+1/2},\\
    \bar{v}_p^{n+1/2} &= \bar{v}_p^{n-1/2} + \bar{E}^n(\bar{x}_p^n),
\end{aligned}
$$
which makes the calculations per particle as fast as possible.

The finite differences formulation reduces to
$$
\begin{aligned}
\bar{\varphi}_{i+1} - 2\bar{\varphi}_{i} + \bar{\varphi}_{i-1} &= \sum_p^{N_p}\bar{S}(\bar{x}-\bar{x}_p) - \frac{N_p}{N_x}\\
\bar{E}_{i+1/2} &= \alpha(\bar{\varphi}_{i+1} - \bar{\varphi}_{i}),
\end{aligned}
$$
where we introduced the only remaining parameter
$$
\alpha = \frac{\omega\,q^2\Delta t^2}{m\,\varepsilon\Delta x}
$$
that describes the strength of the electric field on the particles.

with the transformed shape function
$$
\bar{S}(r) = \begin{cases}
1+r & -1 < r < 0 \\
1-r & 0 < r < 1
\end{cases}.
$$