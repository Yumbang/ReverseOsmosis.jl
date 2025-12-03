# PDE Formulation of Reverse Osmosis Membrane Transport

This document reformulates the iterative steady-state model (Model 1) as a system of partial differential equations (PDEs) in space and time.

## Current Approach (Iterative Steady-State)

The current `element_filtration` function solves for steady-state at each time step:
- Iteratively solves for concentration at each spatial element
- Assumes quasi-steady-state within each time step
- Fouling accumulates between time steps

## PDE Formulation

### Coordinate System

- **x**: Axial position along membrane [m], x ∈ [0, L]
- **t**: Time [s]
- **y**: Cross-membrane direction (not explicitly modeled, averaged)

### State Variables

- **Q(x,t)**: Volumetric flow rate in feed channel [m³/h]
- **C(x,t)**: Salt concentration in feed channel [kg/m³]
- **P(x,t)**: Pressure in feed channel [Pa]
- **R_m(x,t)**: Membrane resistance [Pa·s/m]

### Governing Equations

#### 1. Water Mass Balance
Permeate withdrawal reduces feed flow along the membrane:

```
∂Q/∂x = -v_w(x,t) · W
```

where:
- W = membrane width [m]
- v_w = transmembrane water flux [m/s]

**In dimensional form:**
```
∂Q/∂x [m³/h/m] = -v_w [m/s] · W [m] · 3600 [s/h]
```

#### 2. Salt Mass Balance
Salt is concentrated as water permeates:

```
∂(Q·C)/∂x = -C_p(x,t) · v_w(x,t) · W
```

Expanding:
```
Q·∂C/∂x + C·∂Q/∂x = -C_p · v_w · W
```

Substituting equation (1):
```
Q·∂C/∂x = C_p · v_w · W - C · v_w · W
Q·∂C/∂x = (C_p - C) · v_w · W
Q·∂C/∂x = -C · reject · v_w · W
```

where:
- C_p = (1 - reject)·C [kg/m³] (permeate concentration)
- reject = salt rejection coefficient [-]

**Simplified:**
```
∂C/∂x = -C · reject · v_w · W / Q
```

#### 3. Momentum Equation (Pressure Drop)
Hagen-Poiseuille flow through spacer channel:

```
∂P/∂x = -12·K·μ·u / H²
```

where:
- K = spacer resistance coefficient [-]
- μ = dynamic viscosity [Pa·s]
- u = axial velocity in channel [m/s]
- H = channel height [m]

**Relating velocity to flow rate:**
```
u = Q / (3600 · W · H)  [m/s]
```

**Therefore:**
```
∂P/∂x = -12·K·μ·Q / (3600 · W · H³)
```

#### 4. Constitutive Relations

**Water flux (resistance model):**
```
v_w(x,t) = max(0, [P(x,t) - π(C(x,t),T)] / R_m(x,t))
```

**Osmotic pressure (van't Hoff):**
```
π(C,T) = 2·C·R_gas·(T + 273.15) / M_avg
```
where R_gas = 8314.5 J/(kmol·K), M_avg = 58.44 g/mol for NaCl

**Viscosity (Vogel equation):**
```
μ(T) = 2.414×10⁻⁵ × 10^(247.8/(T+133.15))  [Pa·s]
```

#### 5. Fouling Dynamics
Membrane resistance increases with cumulative flux:

```
∂R_m/∂t = k_fp · v_w(x,t)
```

where k_fp = fouling potential coefficient [Pa·s/m²]

### Complete PDE System

**Spatial evolution (along membrane):**
```
∂Q/∂x = -v_w(P,C,R_m,T) · W

∂C/∂x = -C · reject · v_w(P,C,R_m,T) · W / Q

∂P/∂x = -12·K·μ(T)·Q / (3600 · W · H³)
```

**Temporal evolution (fouling):**
```
∂R_m/∂t = k_fp · v_w(P,C,R_m,T)
```

**Constitutive relations:**
```
v_w = max(0, [P - 2·C·R_gas·(T+273.15)/M_avg] / R_m)
μ = 2.414×10⁻⁵ × 10^(247.8/(T+133.15))
```

### Boundary Conditions

**Inlet (x = 0):**
```
Q(0,t) = Q_feed(t)
C(0,t) = C_feed(t)
P(0,t) = P_feed(t)
```

**Outlet (x = L):**
Natural boundary conditions (determined by integration)

**Initial conditions (t = 0):**
```
R_m(x,0) = R_m_pristine (uniform)
```

### Dimensionless Form

Define characteristic scales:
- L₀ = membrane length [m]
- Q₀ = inlet flow rate [m³/h]
- P₀ = inlet pressure [Pa]
- C₀ = inlet concentration [kg/m³]
- R₀ = pristine resistance [Pa·s/m]

**Dimensionless variables:**
```
x* = x/L₀
Q* = Q/Q₀
C* = C/C₀
P* = P/P₀
R* = R_m/R₀
t* = t·k_fp·v₀/R₀  (characteristic fouling time)
```

where v₀ = P₀/R₀ is characteristic velocity.

**Dimensionless parameters:**
```
Da = W·L₀·v₀/Q₀                    (Damköhler number - permeation/flow ratio)
Re = K·μ·Q₀·L₀/(3600·W·H³·P₀)     (Reynolds-like number - pressure drop)
Π = 2·C₀·R_gas·T/(M_avg·P₀)        (Osmotic pressure ratio)
```

**Dimensionless PDEs:**
```
∂Q*/∂x* = -Da · v_w*

∂C*/∂x* = -Da · reject · C* · v_w* / Q*

∂P*/∂x* = -Re · Q*

∂R*/∂t* = v_w*
```

where `v_w* = max(0, P* - Π·C*) / R*`

## Numerical Solution Methods

### 1. Method of Lines (Current Approach)

**Spatial discretization:**
- Divide [0, L] into N segments: x_i = i·Δx, i = 0..N
- Approximate derivatives: ∂/∂x ≈ (f(x_{i+1}) - f(x_i))/Δx

**Time discretization:**
- Explicit Euler: R^{n+1} = R^n + Δt·k_fp·v_w^n
- This is what the current code does

**Pros:**
- Simple to implement
- Decouples spatial and temporal integration
- Natural for steady-state at each time step

**Cons:**
- Requires iteration at each element for steady-state
- Time step limited by fouling dynamics

### 2. Finite Difference Method

**Spatial grid:** x_i = i·Δx
**Temporal grid:** t_n = n·Δt

**Discretization:**
```
(Q_{i+1} - Q_i)/Δx = -v_w(P_i,C_i,R_i) · W

(C_{i+1} - C_i)/Δx = -C_i · reject · v_w(P_i,C_i,R_i) · W / Q_i

(P_{i+1} - P_i)/Δx = -12·K·μ·Q_i / (3600 · W · H³)

(R_i^{n+1} - R_i^n)/Δt = k_fp · v_w(P_i^n,C_i^n,R_i^n)
```

**Algorithm:**
1. March in space (x-direction) at fixed time
2. March in time with fouling update

### 3. Finite Element Method

**Weak formulation** using test functions φ:

```
∫[∂Q/∂x + v_w·W]·φ dx = 0
∫[∂(QC)/∂x + C_p·v_w·W]·φ dx = 0
∫[∂P/∂x + 12Kμu/H²]·φ dx = 0
```

**Pros:**
- Better for complex geometries
- Higher-order accuracy

**Cons:**
- More complex implementation
- Overkill for 1D problems

### 4. ODE-based Approach (DifferentialEquations.jl)

Treat spatial position as "time" and solve ODEs:

```julia
using DifferentialEquations

function RO_spatial!(du, u, p, x)
    Q, C, P = u
    W, H, reject, K, μ, R_m, T, M_avg = p

    # Calculate water flux
    π = 2 * C * 8314.5 * (T + 273.15) / 58.44
    v_w = max(0.0, (P - π) / R_m)

    # Spatial derivatives
    du[1] = -v_w * W * 3600  # dQ/dx
    du[2] = -C * reject * v_w * W / Q * 3600  # dC/dx
    du[3] = -12 * K * μ * Q / (3600 * W * H^3)  # dP/dx
end

# Solve from inlet (x=0) to outlet (x=L)
u0 = [Q_feed, C_feed, P_feed]
prob = ODEProblem(RO_spatial!, u0, (0.0, L), params)
sol = solve(prob, Tsit5())
```

**Pros:**
- Leverages sophisticated ODE solvers
- Adaptive stepping
- Automatic stiffness detection

**Cons:**
- Requires DifferentialEquations.jl dependency
- Less control over spatial discretization

## Implementation Example

Here's a complete PDE-based implementation:

```julia
using DifferentialEquations

"""
PDE-based single membrane element simulation.
Solves spatial ODEs for steady-state, then updates fouling.
"""
function element_filtration_pde(
    feed::Water,
    element::MembraneElement;
    dt::Union{Float64, Nothing}=nothing
)
    # Extract parameters
    W = element.width
    H = element.height
    L = element.dx
    K = element.spacer_resistance
    reject = element.salt_rejection
    R_m = element.R_m
    k_fp = element.k_fp
    T = feed.T

    # Viscosity
    μ = 2.414e-5 * 10^(247.8 / (T + 273.15 - 140))

    # Define spatial ODE system
    function spatial_ode!(du, u, p, x)
        Q, C, P = u
        R_m, T = p

        # Osmotic pressure
        π = 2 / 58.44 * C * 8.3145e3 * (T + 273.15)

        # Water flux
        v_w = max(0.0, (P - π) / R_m)

        # Spatial derivatives
        du[1] = -v_w * W * 3600              # dQ/dx [m³/h/m]
        du[2] = -C * reject * v_w * W * 3600 / Q  # dC/dx [kg/m³/m]
        du[3] = -12 * K * μ * Q / (3600 * W * H^3)  # dP/dx [Pa/m]

        return nothing
    end

    # Initial conditions (inlet)
    u0 = [feed.Q, feed.C, feed.P]

    # Parameters
    params = (R_m, T)

    # Solve spatial ODE from x=0 to x=L
    prob = ODEProblem(spatial_ode!, u0, (0.0, L), params)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)

    # Extract outlet conditions
    Q_out, C_out, P_out = sol.u[end]

    # Calculate average water flux for fouling
    v_w_avg = -(feed.Q - Q_out) / (W * L * 3600)

    # Permeate properties
    Q_p = feed.Q - Q_out
    C_p = (1 - reject) * feed.C
    P_p = 1e5  # Atmospheric

    # Fouling
    if isnothing(dt)
        ΔR_m = k_fp * v_w_avg
    else
        ΔR_m = k_fp * v_w_avg * dt
    end

    permeate = Water(Q_p, T, C_p, P_p)
    brine = Water(Q_out, T, C_out, P_out)

    return (brine, permeate, ΔR_m)
end
```

## Advantages of PDE Formulation

1. **Physical clarity**: Explicit conservation laws
2. **Flexibility**: Easy to modify physics (add dispersion, reaction, etc.)
3. **Advanced numerics**: Leverage sophisticated PDE solvers
4. **Adaptive refinement**: Automatic mesh adaptation where needed
5. **Uncertainty quantification**: Easier for sensitivity analysis
6. **Multi-physics coupling**: Natural framework for adding energy equations, etc.

## Disadvantages

1. **Computational cost**: ODE solvers may be slower than direct iteration
2. **Dependencies**: Requires DifferentialEquations.jl
3. **Learning curve**: More complex for simple cases
4. **Debugging**: Less transparent than explicit iteration

## Extensions

### 1. Add Concentration Polarization Layer

Include boundary layer mass transfer:

```
∂C_membrane/∂x = C_bulk · exp(v_w / k_mt)
```

where k_mt is from Sherwood correlation.

### 2. Energy Equation

Add temperature variation:

```
∂T/∂x = heat_transfer_terms
```

### 3. Unsteady Feed Conditions

Solve full space-time PDE:

```
∂Q/∂t + u·∂Q/∂x = -v_w·W
```

This requires 2D (x,t) finite difference or moving mesh methods.

### 4. Multi-Stage Systems

Couple multiple membrane elements as boundary conditions:

```
Q_inlet^{module_i+1} = Q_outlet^{module_i}
```

## Comparison with Current Code

| Aspect | Current (Iterative) | PDE-based |
|--------|---------------------|-----------|
| Speed | Fast (optimized) | Moderate (depends on solver) |
| Accuracy | High (tight tolerance) | Configurable |
| Flexibility | Limited | High |
| Complexity | Low | Moderate |
| Dependencies | None | DifferentialEquations.jl |
| Physical insight | Implicit | Explicit |

## Recommendations

1. **Keep current approach** for production (fast, well-tested)
2. **Use PDE formulation** for:
   - Research and model development
   - Verification of current code
   - Extensions (multi-physics, optimization, etc.)
   - Publications (clearer presentation)
3. **Hybrid approach**: PDE for complex cases, iteration for simple cases

## Further Reading

1. **Membrane transport**: "Membrane Technology and Applications" by Baker
2. **Numerical methods**: "Numerical Recipes" for ODE/PDE methods
3. **DifferentialEquations.jl**: Official documentation
4. **Dimensionless analysis**: Buckingham π theorem

## References

- Film theory: Brian, P.L.T. (1965). "Concentration polarization in reverse osmosis desalination with variable flux and incomplete salt rejection." *Ind. Eng. Chem. Fundam.*
- RO modeling: Geraldes, V., et al. (2001). "Flow and mass transfer modelling of nanofiltration." *Journal of Membrane Science*
