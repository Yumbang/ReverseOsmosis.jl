# ReverseOsmosis.jl Source Code

This directory contains the core implementation of the ReverseOsmosis.jl package for simulating reverse osmosis membrane processes.

## Architecture Overview

The codebase is organized in a hierarchical structure mirroring the physical system:

```
Water Streams (fluid.jl)
    ↓
Membrane Elements (membrane.jl)
    ↓
Filtration Physics (filtration.jl)
    ↓
Process Models (Process/)
```

## File Descriptions

### Core Data Structures

#### `fluid.jl`
Defines water stream representations and operations.

**Key Types:**
- `Water`: Basic water stream (Q, T, C, P, m_avg)
- `Water2`: Extended stream with foulant concentration tracking

**Key Functions:**
- `mix()`: Mass-weighted stream mixing
- `pressurize()`: Apply pressure to streams
- `profile_water()`: Export to DataFrames for analysis

#### `membrane.jl`
Defines membrane hierarchy and constructors.

**Hierarchical Structure:**
```
MembraneElement/Element2  (single discretized segment)
    ↓
MembraneModule  (array of elements)
    ↓
PressureVessel  (array of modules)
```

**Key Types:**
- `MembraneElement`: Simple resistance model with salt rejection
- `MembraneElement2`: Advanced solution-diffusion model with cake layer
- `MembraneModule`: Collection of elements in series
- `PressureVessel`: Collection of modules in series

**Key Functions:**
- `pristine_membrane()`, `pristine_vessel()`: Constructors for pristine membranes
- `pristine_membrane2()`, `pristine_vessel2()`: Constructors for Element2 type
- `profile_membrane()`, `profile_vessel()`: DataFrame export with units
- `foul!()`: Apply fouling (mutates membrane properties)

### Physics & Simulation

#### `filtration.jl`
Core filtration algorithms and physics models.

**Module Structure:**
- `Filtration`: Submodule with physics implementations

**Key Functions:**
- `osmo_p()`: Osmotic pressure (van't Hoff equation)
- `pump()`: Pump modeling with energy consumption
- `element_filtration()`: Single element simulation (MembraneElement)
- `element_filtration2()`: Advanced simulation (MembraneElement2)
- `module_filtration()`: Module-level cascade
- `vessel_filtration()`: Vessel-level cascade

**Models Implemented:**
1. **Resistance Model** (element_filtration):
   - v_w = (P - π) / R_m
   - Fixed salt rejection coefficient
   - Simple fouling: ΔR_m = k_fp × v_w

2. **Solution-Diffusion Model** (element_filtration2):
   - Water flux: v_w = A × (P - π)
   - Salt flux: v_s = B × C_membrane
   - Cake-enhanced concentration polarization (CEMT)
   - Critical flux fouling model
   - Carman-Kozeny equation for cake resistance

**Performance Note:**
All filtration functions avoid Unitful.jl in hot loops for ~100x speedup. Units handled at API boundaries only.

#### `ro_units.jl`
Custom units for RO processes (currently commented out).

- `lmh`: Liters per m² per hour (flux)
- `lmhbar`: Liters per m² per hour per bar (permeability)

### Process Models

#### `Process/ROProcess.jl`
Container module for different process configurations.

**Current Processes:**
- SinglePass: Single-pass RO system

#### `Process/SinglePass.jl`
Single-pass RO process implementation.

**Key Type:**
- `SinglePassRO`: Stateful process structure

**Key Function:**
- `process_singlepass_RO!()`: Simulate one time step
  - Pressurizes feed
  - Runs vessel_filtration
  - Applies fouling (if enabled)
  - Returns brine, permeate, power consumption

**Process Flow:**
```
Unpressurized Feed
    ↓
High-Pressure Pump (with efficiency)
    ↓
Pressure Vessel (multiple modules in series)
    ↓
Final Brine + Mixed Permeate
```

## Design Patterns

### Type Dispatch
The codebase uses Julia's multiple dispatch extensively:
- `Water` + `MembraneElement` → `element_filtration`
- `Water2` + `MembraneElement2` → `element_filtration2`

### Mutable State
- Membrane structures are mutable to track fouling over time
- `foul!()` functions mutate membrane properties in-place
- Process structures (`SinglePassRO`) maintain state across time steps

### Performance Optimization
- Hot loops (filtration algorithms) avoid Unitful.jl
- Type-stable implementations with explicit type annotations
- Iterative solvers with convergence tolerance (1e-6 or 1e-4)

### Units Handling
- API boundaries use Unitful.jl for type safety
- Internal calculations use raw Float64 (SI units)
- Profile functions attach units for DataFrame export

## Typical Usage Flow

```julia
using ReverseOsmosis

# 1. Create membrane
vessel = pristine_vessel(
    7,      # 7 modules
    50,     # 50 segments per module
    1.016,  # 1.016 m length
    7e-4,   # 0.7 mm channel height
    37.0,   # 37 m² area
    6.8e10, 0.67e9, 16.0, 0.995  # membrane properties
)

# 2. Create process
process = SinglePassRO(0.8, vessel)  # 80% pump efficiency

# 3. Define feed
feed = Water(10.0, 25.0, 35.0, 1e5)  # 10 m³/h, 25°C, 35 kg/m³ TDS, 1 bar

# 4. Simulate
brine, permeate, _, _, power = process_singlepass_RO!(
    feed,
    55e5u"Pa";          # 55 bar pressure
    process=process,
    dt=3600.0u"s",      # 1 hour time step
    fouling=true
)

# 5. Analyze results
println("Permeate flow: $(permeate.Q) m³/h")
println("Permeate TDS: $(permeate.C) kg/m³")
println("Power: $(power)")
```

## Extension Points

To add new features:

1. **New fouling models**: Modify `element_filtration` or create `element_filtration3`
2. **New process configurations**: Create new module in `Process/`
3. **New membrane types**: Add `MembraneElement3` with new properties
4. **Custom units**: Uncomment and extend `ro_units.jl`

## Testing

Currently no formal test suite. Examples in `../examples/` serve as integration tests.

## References

Key scientific papers implemented:
- Chong et al. (2008), J. Membrane Sci. 314:89-95 - Critical flux model
- Jiang et al. (2021), Desalination 510:115289 - Cake-enhanced mass transfer
