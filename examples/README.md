# Examples Directory

**⚠️ WARNING: Examples in this directory are currently outdated and non-functional.**

The example files in this directory reference semi-batch RO functionality that has been moved to a separate repository: [SemiBatchReverseOsmosis.jl](https://github.com/Yumbang/SemiBatchReverseOsmosis.git)

## Current Status

### `test_semibatch.jl` ❌
- **Status**: Non-functional
- **Issue**: References undefined types `CirculationPipe`, `SemiBatchRO`, and function `process_semibatch_RO!`
- **Reason**: Semi-batch code moved to separate repository

### `test_ode.jl` ❌
- **Status**: Non-functional
- **Issue**: Implements ODE-based semi-batch RO simulation
- **Dependencies**: Requires `DifferentialEquations.jl`, `DiffEqCallbacks.jl` (not in current dependencies)
- **Reason**: Semi-batch code moved to separate repository

## Working with Single-Pass RO

For functional examples using the current codebase, here's a basic single-pass RO simulation:

```julia
using ReverseOsmosis
using Unitful

# Create a pristine 7-module pressure vessel
# Each module: 1.016 m length, 37 m² area, 50 axial segments
vessel = pristine_vessel(
    7,          # Number of modules
    50,         # Segments per module
    1.016,      # Module length [m]
    7e-4,       # Channel height [m]
    37.0,       # Membrane width [m] (rolled area ≈ 37 m²)
    6.8e10,     # Membrane resistance [Pa·s/m]
    0.67e9,     # Fouling potential [Pa·s/m²]
    16.0,       # Spacer resistance [-]
    0.995       # Salt rejection [-]
)

# Create single-pass RO process with 80% pump efficiency
process = SinglePassRO(0.8, vessel)

# Define seawater feed at atmospheric pressure
feed = Water(
    10.0,   # Flow rate [m³/h]
    25.0,   # Temperature [°C]
    35.0,   # TDS concentration [kg/m³] (~35 g/L seawater)
    1e5     # Pressure [Pa] (atmospheric)
)

# Simulate one hour of operation at 55 bar
brine, permeate, _, _, power = process_singlepass_RO!(
    feed,
    55e5u"Pa";              # Operating pressure [Pa]
    process=process,
    dt=3600.0u"s",          # Time step (1 hour) [s]
    fouling=true,           # Apply fouling
    profile_process=false   # Don't return detailed profiles
)

# Display results
println("=== Single-Pass RO Results ===")
println("Feed: ", feed)
println("Brine: ", brine)
println("Permeate: ", permeate)
println("Power consumption: ", power)
println("Recovery: $(permeate.Q / feed.Q * 100)%")
println("Salt rejection: $((1 - permeate.C / feed.C) * 100)%")
```

## Long-Term Fouling Simulation

To simulate fouling over extended operation:

```julia
using ReverseOsmosis
using Unitful
using DataFrames
using Plots

# Initialize system
vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
process = SinglePassRO(0.8, vessel)

# Storage for results
results = []

# Simulate 100 hours
for hour in 1:100
    feed = Water(10.0, 25.0, 35.0, 1e5)

    brine, permeate, _, _, power = process_singlepass_RO!(
        feed, 55e5u"Pa";
        process=process,
        dt=3600.0u"s",
        fouling=true
    )

    # Store results
    push!(results, (
        hour = hour,
        permeate_flow = permeate.Q,
        permeate_TDS = permeate.C,
        power = ustrip(u"kW", power),
        recovery = permeate.Q / feed.Q * 100
    ))
end

# Convert to DataFrame and plot
df = DataFrame(results)

plot(df.hour, df.permeate_flow,
     xlabel="Time [hours]",
     ylabel="Permeate Flow [m³/h]",
     title="Membrane Fouling Effect",
     label="Permeate Production")
```

## Advanced: Spatial Profiling

To analyze spatial variations along the membrane:

```julia
# Enable detailed profiling
brine, permeate, brines_profile, permeates_profile, power = process_singlepass_RO!(
    feed, 55e5u"Pa";
    process=process,
    dt=3600.0u"s",
    profile_process=true  # ← Enable profiling
)

# brines_profile is a DataFrame with one row per membrane element
# Columns: Q, T, C, P (all with units)

using Plots

# Plot concentration polarization along vessel
plot(1:nrow(brines_profile), brines_profile.C,
     xlabel="Element Number",
     ylabel="Brine Concentration",
     title="Concentration Profile Along Vessel",
     label="Brine TDS")
```

## Using Advanced Model (Water2 / MembraneElement2)

For cake layer fouling simulations:

```julia
using ReverseOsmosis

# Create vessel with Element2 type
vessel = pristine_vessel2(
    7,          # Number of modules
    50,         # Segments per module
    1.016,      # Module length [m]
    7e-4,       # Channel height [m]
    37.0,       # Membrane width [m]
    1.5e15,     # Membrane resistance [m⁻¹] (different units!)
    3.5e-8,     # Salt permeability [kg/m²·s]
    0.8,        # Fouling propensity [-]
    16.0        # Spacer resistance [-]
)

process = SinglePassRO(0.8, vessel)

# Feed with foulant concentration
feed = Water2(
    10.0,   # Flow rate [m³/h]
    25.0,   # Temperature [°C]
    35.0,   # Salt concentration [kg/m³]
    0.1,    # Foulant concentration [kg/m³]
    1e5     # Pressure [Pa]
)

# Simulate (note: uses element_filtration2 internally)
brine, permeate, _, _, power = process_singlepass_RO!(
    feed, 55e5u"Pa";
    process=process,
    dt=3600.0u"s",
    fouling=true
)

# Check cake layer accumulation
vessel_profile = profile_vessel(process.pressure_vessel)
println("Maximum cake thickness: ", maximum(vessel_profile.cake_thickness))
println("Maximum cake resistance: ", maximum(vessel_profile.R_c))
```

## For Semi-Batch RO

If you need semi-batch RO functionality (closed-circuit with brine recirculation), please refer to the separate repository:

**[SemiBatchReverseOsmosis.jl](https://github.com/Yumbang/SemiBatchReverseOsmosis.git)**

That repository contains:
- Circulation pipe modeling
- Semi-batch process structures
- CC (closed-circuit) and purge modes
- ODE-based dynamic simulations
- Recovery optimization examples

## Contributing Examples

If you create new working examples, please consider contributing them back to the repository via pull request!
