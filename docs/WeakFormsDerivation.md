In order to derive weak form of any equatiion we need to multiply both right hands side and left hands side to corresponding test function.

<!-- ## Momentum equation 
 -->

## Heat transfer
First we write heat-transfer equation in **strong form**.
```math
\frac{\partial T}{\partial t} 
= 
  \frac{k}{\rho C_p} \nabla^2\, T 
- (\mathbf{u}\cdot\nabla)\, T 
```


Rewrite left hand side as follows:
```math
\frac{\partial T}{\partial t} = \frac{1}{\Delta t} \left(T^{(n+1)}-T^{(n)}\right)
```

Here, $(n)$ and $(n+1)$ denote the temperature field at the current and next time step, respectively.

Multiplying by test function $w$:
```math
\frac{1}{\Delta t} \left(T^{(n+1)}-T^{(n)}\right) \cdot w
= \frac{k}{\rho C_p} \nabla^2  T  \cdot w
- (\mathbf{u}\cdot\nabla)\, T \cdot w
```

Then integrate over domain $\Omega$:
```math
\frac{1}{\Delta t} \int\limits_\Omega \left(T^{(n+1)}-T^{(n)}\right) \cdot w
= \frac{k}{\rho C_p} \int\limits_\Omega  \nabla^2  T  \cdot w
- \int\limits_\Omega (\mathbf{u}\cdot\nabla)\, T \cdot w
```

### Integration by parts and 
Let's take a closure look at right hands side's first term. We can notice that it's a part of integration by parts equation. Let's write its obviously:
```math
\int\limits_\Omega \nabla \cdot \left(\nabla\,T\cdot\,w\right)
= \int\limits_\Omega  \nabla^2\,T\cdot\,w
+ \int\limits_\Omega  (\nabla\,T) \cdot (\nabla\,w)
```

Rearrange for $\int\limits_\Omega  \nabla^2\,T\cdot\,w$:
```math
\int\limits_\Omega  \nabla^2\,T\cdot\,w 
= \int\limits_\Omega \nabla \cdot\left(\nabla\,T\cdot\,w\right)
- \int\limits_\Omega  (\nabla\,T) \cdot (\nabla\,w)
```

### Application of Gauss divergence theorem
Using the divergence theorem, we can rewrite equation as follows:
```math
\int\limits_\Omega  \nabla^2\,T\cdot\,w 
= 
  \int\limits_{\partial\Omega} \left(\nabla\,T\cdot\,w\right) 
  \cdot \mathbf{n}
- \int\limits_\Omega  (\nabla\,T) \cdot (\nabla\,w)
```

where $\partial\Omega$ denotes the boundary of the domain $\Omega$ and $\mathbf{n}$ is the outward unit normal vector on the boundary. Assuming that the temperature is known on the boundary, which means that the test function $w$ equals zero on the boundary, the term becomes zero, and we are left with:
```math
\int\limits_\Omega  \nabla^2\,T\cdot\,w 
=
- \int\limits_\Omega (\nabla\,T) \cdot (\nabla\,w)
```

### Weak form for heat transfer equation
Substituting back into the equation, we obtain the weak form of the heat-transfer equation as:
```math
\frac{1}{\Delta t} \int\limits_\Omega \left(T^{(n+1)}-T^{(n)}\right) \cdot w
=
- \frac{k}{\rho C_p} \int\limits_\Omega (\nabla\,T) \cdot (\nabla\,w) 
- \int\limits_\Omega (\mathbf{u}\cdot\nabla)\, T \cdot w
```
where $n$ is the current time step and $n+1$ is the next time step, and $w$ is the test function.


<!-- ### Side notes
Heat-transfer equation is typical scalar [convection-diffusion equation](https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation) that writes the following way:
$$
\frac{\partial c}{\partial t} = \nabla\cdot\left(D\nabla\,c\right) - \nabla\cdot(\mathbf{v}c) + R
$$
where: -->


