## Physical model

The evolution of the density function ``f(x, v, t)`` of a charged species is given by the Vlasov equation

$
\frac{\partial f}{\partial t} + v\cdot\frac{\partial f}{\partial t} + \frac{q}{m}E\cdot\frac{\partial f}{\partial v} = 0,
$

with the species charge ``q`` and mass ``m``.

The electric field is given by the electrostatic Poisson equation

$
\varphi_{xx} = - \frac{\rho}{\varepsilon}\\
E = - \varphi_x
$

where $\varphi$ is the electric potential and $\varepsilon$ is the absolute permittivity.

The charge density is given by

$
\rho = q \int f\,\mathrm{d}v - Q
$

with $Q=\frac{q}{L}\int\int f\,\mathrm{d}x\,\mathrm{d}v$ being a uniform background charge in order to ensure quasi-neutrality.