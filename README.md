# ReverseOsmosis.jl

A Julia package for simulating reverse osmosis (RO) membrane filtration processes with detailed physics modeling and fouling dynamics.

[![Julia Version](https://img.shields.io/badge/julia-v1.10+-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

ReverseOsmosis.jl provides a comprehensive framework for simulating single-pass reverse osmosis systems with:

- **Physics-based modeling**: Implementation of solution-diffusion theory, concentration polarization, and membrane fouling
- **Spatial discretization**: Axial profiles of concentration, pressure, and fouling along membrane modules
- **Flexible architecture**: Support for different membrane models and fouling mechanisms
- **Performance optimized**: ~100x faster than naive implementations through careful avoidance of Unitful.jl in hot loops
- **Type-safe interfaces**: Unitful.jl integration at API boundaries for dimensional correctness

## Features

### Water Stream Modeling
- `Water`: Basic water streams with flow, temperature, concentration, and pressure
- `Water2`: Extended streams with foulant particle tracking
- Mass-weighted mixing and pressurization operations

### Membrane Hierarchy
- **MembraneElement**: Discretized membrane segments for spatial resolution
- **MembraneModule**: Collections of elements representing spiral-wound modules
- **PressureVessel**: Series arrangements of modules (typical: 7 modules/vessel)

### Physical Models

#### Model 1: Resistance Model (MembraneElement)
Simple and fast model suitable for preliminary design:
- Water flux: `v_w = (P - π) / R_m`
- Fixed salt rejection coefficient
- Fouling: `ΔR_m = k_fp × v_w`

#### Model 2: Solution-Diffusion with Cake Layer (MembraneElement2)
Advanced model for detailed fouling studies:
- Water flux: `v_w = A × (P - π)` with temperature correction
- Salt flux: `v_s = B × C_membrane` (solution-diffusion)
- Cake-enhanced concentration polarization (CEMT model)
- Critical flux-based cake deposition
- Carman-Kozeny equation for cake resistance

### Process Configurations
- **Single-Pass RO**: Feed → Pump → Membrane → Permeate + Brine
- **Semi-Batch RO**: Moved to [SemiBatchReverseOsmosis.jl](https://github.com/Yumbang/SemiBatchReverseOsmosis.git)

## Installation

```julia
# In Julia REPL, press ] to enter Pkg mode
pkg> add https://github.com/Yumbang/ReverseOsmosis.jl.git

# Or using Pkg
using Pkg
Pkg.add(url="https://github.com/Yumbang/ReverseOsmosis.jl.git")
```

## Quick Start

### Basic Single-Pass RO Simulation

```julia
using ReverseOsmosis
using Unitful

# Create a 7-module pressure vessel (typical industrial configuration)
# Dow Filmtec SW30HR-380 specifications
vessel = pristine_vessel(
    7,          # 7 modules in series
    50,         # 50 axial segments per module
    1.016,      # 1.016 m length (40 inches)
    7e-4,       # 0.7 mm channel height
    37.0,       # 37 m² membrane area
    6.8e10,     # Membrane resistance [Pa·s/m]
    0.67e9,     # Fouling potential [Pa·s/m²]
    16.0,       # Spacer resistance coefficient
    0.995       # 99.5% salt rejection
)

# Create process with 80% pump efficiency
process = SinglePassRO(0.8, vessel)

# Define seawater feed
feed = Water(
    10.0,   # 10 m³/h flow rate
    25.0,   # 25°C temperature
    35.0,   # 35 kg/m³ TDS (typical seawater)
    1e5     # Atmospheric pressure [Pa]
)

# Simulate 1 hour at 55 bar operating pressure
brine, permeate, _, _, power = process_singlepass_RO!(
    feed,
    55e5u"Pa";          # 55 bar
    process=process,
    dt=3600.0u"s",      # 1 hour time step
    fouling=true        # Enable fouling
)

# Results
println("Permeate production: $(permeate.Q) m³/h")
println("Permeate TDS: $(permeate.C) kg/m³")
println("Recovery rate: $((permeate.Q / feed.Q) * 100)%")
println("Power consumption: $(power)")
```

### Long-Term Fouling Simulation

```julia
using ReverseOsmosis, Unitful, DataFrames

vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
process = SinglePassRO(0.8, vessel)

results = DataFrame(
    hour = Int[],
    permeate_flow = Float64[],
    recovery = Float64[],
    power_kW = Float64[]
)

# Simulate 168 hours (1 week)
for hour in 1:168
    feed = Water(10.0, 25.0, 35.0, 1e5)

    brine, permeate, _, _, power = process_singlepass_RO!(
        feed, 55e5u"Pa";
        process=process,
        dt=3600.0u"s",
        fouling=true
    )

    push!(results, (
        hour,
        permeate.Q,
        (permeate.Q / feed.Q) * 100,
        ustrip(u"kW", power)
    ))
end

# Analyze performance decline
println("Initial permeate: $(results.permeate_flow[1]) m³/h")
println("Final permeate: $(results.permeate_flow[end]) m³/h")
println("Flux decline: $((1 - results.permeate_flow[end]/results.permeate_flow[1]) * 100)%")
```

### Advanced: Cake Layer Fouling

```julia
using ReverseOsmosis

# Use MembraneElement2 with solution-diffusion model
vessel = pristine_vessel2(
    7, 50, 1.016, 7e-4, 37.0,
    1.5e15,     # Membrane resistance [m⁻¹] (different units!)
    3.5e-8,     # Salt permeability [kg/m²·s]
    0.8,        # Fouling propensity
    16.0        # Spacer resistance
)

process = SinglePassRO(0.8, vessel)

# Feed with foulant particles
feed = Water2(
    10.0,   # Flow rate
    25.0,   # Temperature
    35.0,   # Salt concentration
    0.1,    # Foulant concentration [kg/m³]
    1e5     # Pressure
)

# Simulate with cake layer tracking
brine, permeate, _, _, power = process_singlepass_RO!(
    feed, 55e5u"Pa";
    process=process,
    dt=3600.0u"s",
    fouling=true
)

# Check cake accumulation
profile = profile_vessel(process.pressure_vessel)
println("Max cake thickness: $(maximum(profile.cake_thickness))")
println("Max cake resistance: $(maximum(profile.R_c))")
```

## Scientific Background

### Osmotic Pressure
Van't Hoff equation assuming complete dissociation:
```
π = (2 × C × R × T) / M_avg
```
where C is salt concentration, R is gas constant, T is temperature, M_avg is molecular weight.

### Concentration Polarization
Film theory with cake-enhanced mass transfer:
```
C_membrane = C_bulk × exp(v_w / k_CEMT)
```
where k_CEMT accounts for both bulk mass transfer and diffusion through cake layer.

### Fouling Models

**Simple Resistance Model:**
```
dR_m/dt = k_fp × v_w
```
Membrane resistance increases proportionally to water flux.

**Cake Layer Model:**
```
dδ/dt = m_fc / (ρ_p × (1 - ε_p))
R_c = 180 × (1 - ε_p) / (ρ_p × d_p² × ε_p³) × m_cake
```
where δ is cake thickness, m_fc is foulant deposition flux (from critical flux model), and R_c follows Carman-Kozeny equation.

## Performance

The codebase is optimized for performance:
- Element filtration functions avoid Unitful.jl in hot loops (~100x speedup)
- Type-stable implementations with explicit annotations
- Iterative solvers with convergence tolerance (1e-6 or 1e-4)
- Efficient memory usage through in-place mutations (`foul!`)

Typical performance on modern hardware:
- Single element: ~10 μs
- Full vessel (7 modules × 50 segments): ~3 ms
- 168-hour simulation (fouling enabled): ~0.5 s

## Documentation

- **[src/README.md](src/README.md)**: Detailed code architecture and API reference
- **[examples/README.md](examples/README.md)**: Usage examples and tutorials
- **Inline documentation**: All functions have comprehensive docstrings

## Project Structure

```
ReverseOsmosis.jl/
├── src/
│   ├── ReverseOsmosis.jl      # Main module
│   ├── fluid.jl                # Water stream types and operations
│   ├── membrane.jl             # Membrane hierarchy and constructors
│   ├── filtration.jl           # Core physics and filtration algorithms
│   ├── ro_units.jl             # Custom units (optional)
│   └── Process/
│       ├── ROProcess.jl        # Process container module
│       └── SinglePass.jl       # Single-pass RO implementation
├── examples/
│   ├── README.md               # Example usage guide
│   ├── test_semibatch.jl       # ⚠️ Outdated (moved to separate repo)
│   └── test_ode.jl             # ⚠️ Outdated (moved to separate repo)
├── Project.toml                # Package dependencies
└── README.md                   # This file
```

## Dependencies

- **Julia**: ≥1.10
- **DataFrames.jl**: 1.7+
- **Unitful.jl**: 1.22+
- **Printf**: Standard library

## Related Projects

- **[SemiBatchReverseOsmosis.jl](https://github.com/Yumbang/SemiBatchReverseOsmosis.git)**: Semi-batch RO with brine recirculation, ODE-based dynamics

## Scientific References

Key papers implemented in this package:

1. **Critical Flux Model**: Chong, T.H., et al. (2008). "Combined effect of membrane fouling and concentration polarization on flux behavior in crossflow nanofiltration." *Journal of Membrane Science*, 314(1-2), 89-95. [DOI:10.1016/j.memsci.2008.01.030](https://doi.org/10.1016/j.memsci.2008.01.030)

2. **Cake-Enhanced Mass Transfer**: Jiang, S., et al. (2021). "A comprehensive CFD-based model for the prediction of fouling development in a spiral-wound reverse osmosis membrane." *Desalination*, 510, 115289. [DOI:10.1016/j.desal.2021.115289](https://doi.org/10.1016/j.desal.2021.115289)

3. **Solution-Diffusion Theory**: See standard membrane transport references (Geankoplis, Baker & Wijmans, etc.)

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

Areas for contribution:
- Additional process configurations (multi-stage, ERD systems)
- New fouling models
- Validation against experimental data
- Test suite development
- Performance benchmarks
- Additional examples

## Citation

If you use this package in your research, please cite:

```bibtex
@software{reverseosmosisjl,
  author = {Yumbang},
  title = {ReverseOsmosis.jl: A Julia Package for Reverse Osmosis Simulation},
  year = {2024},
  url = {https://github.com/Yumbang/ReverseOsmosis.jl}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Yumbang ([@Yumbang](https://github.com/Yumbang))

## Acknowledgments

This package was developed for research in membrane desalination and water treatment processes.

---

**Note**: Semi-batch RO functionality has been moved to a separate repository. For closed-circuit operation with brine recirculation, see [SemiBatchReverseOsmosis.jl](https://github.com/Yumbang/SemiBatchReverseOsmosis.git).
