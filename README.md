


## Algorithm
(need some text here)

### Define simulation parameters
During this step, we define the simulation parameters, including:
- $\Delta t $ - timestep length ($\mathrm{s}$), which should be selected carefully as it significantly affects the stability of the simulation. To determine the appropriate value, the CFL condition must be taken into account.
- $N$ - the number of timesteps.

### Define fluids properties
Some physical properties of the fluid are:
- Density ($\rho$), which is the amount of mass per unit volume of the substance, measured in $\mathrm{kg/m^3}$.
- Dynamic viscosity ($\mu$), measured in $\mathrm{Pa\cdot s}$, which is a measure of the fluid's resistance to flow.
- Thermal conductivity ($k$), measured in $\mathrm{W/(m\cdot K)}$, which is a measure of the material's ability to conduct heat.
- Specific heat capacity ($C_p$), measured in $\mathrm{J/(kg\cdot K)}$, which is the amount of heat required to raise the temperature of one unit of mass of the fluid by one degree Kelvin.

### Define mesh
(need some text here)

### Time loop (Simulation)

#### Solve for tentitive velocity

#### Solve pressure Poisson equation

### Update velocities

### Calculate heat transfer


## Brief explanation
### Domain 
### Boundary Conditions (BC)
The purpose of boundary conditions is to determine unknown function values on domain's boundary.
- Dirichlet boundary conditions, which specify the value of the function at the boundary. For example, $T(0.0,y) = 200.0$.
- Neumann boundary conditions, which specify the normal derivative of the function at the boundary. For example, $\frac{\partial T}{\partial x}(L,y) = 0.0$.

(about combining BCs and creating list of it)

### Initial Conditions (IC)
The initial conditions describe the field values at the start of the simulation, i.e., at $t=0$. To set the initial state for a field in FEniCS, we can use the following syntax:

```Python
# Set the initial circular conditions with a center at (0.5, 0.5) and a magnitude of 60
t_prev.interpolate(
    Expression(
        '60 * (0.5*(1 - tanh(8*(sqrt(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2)) - 0.1))))',
     degree=1)
)
```
Here, `t_prev` is assigned an initial condition using `Expression()` with an expression evaluated for each point in the domain. The example above defines circular initial conditions with a center at $(0.5,0.5)$ and a magnitude of 60.

## Derivation of weak forms
Complete derrivation see in []
