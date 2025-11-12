"""
    Filtration

Core module implementing reverse osmosis filtration algorithms.

This module contains the fundamental physics-based models for simulating
membrane filtration processes, including:
- Osmotic pressure calculations (van't Hoff equation)
- Pump energy modeling
- Element-level filtration with concentration polarization
- Module and vessel-level filtration cascades
- Fouling dynamics (resistance increase and cake layer formation)

# Performance Note
Filtration functions intentionally avoid Unitful.jl in hot loops for
~100x performance improvement. Units are handled at API boundaries only.

# Exported Functions
- `osmo_p`: Osmotic pressure calculation
- `pump`: Pump modeling with energy consumption
- `element_filtration`: Single element filtration (MembraneElement)
- `element_filtration2`: Advanced element filtration (MembraneElement2)
- `module_filtration`: Module-level filtration cascade
- `vessel_filtration`: Vessel-level filtration cascade
"""
module Filtration

export
    # Filtration
    osmo_p, pump,
    element_filtration, element_filtration2, module_filtration, vessel_filtration

using Unitful
using ..ReverseOsmosis: Water, Water2, pressurize, profile_water, mix
using ..ReverseOsmosis: MembraneElement, MembraneElement2, MembraneModule, PressureVessel, profile_membrane, pristine_membrane

# ===================================================================
# Pump Modeling
# ===================================================================

"""
    pump(feed::Union{Water, Water2}, pressure_to_apply::Float64; efficiency::Float64=0.8)

Model high-pressure pump with energy consumption calculation.

# Arguments
- `feed::Union{Water, Water2}`: Feed water stream
- `pressure_to_apply::Float64`: Pressure increase [Pa]
- `efficiency::Float64=0.8`: Pump isentropic efficiency (0-1)

# Returns
- Tuple of (pressurized_water, power_consumption)
  - `pressurized_water`: Water stream with increased pressure
  - `power_consumption`: Electrical power consumption [W] with Unitful units

# Formula
Power = (Volumetric flow × Pressure increase) / efficiency

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 1e5)
pressurized, power = pump(feed, 50e5; efficiency=0.8)
# Returns water at 51 bar and power consumption in watts
```
"""
function pump(feed::Union{Water, Water2}, pressure_to_apply::Float64; efficiency::Float64=0.8)
    @assert (0.0 ≤ efficiency ≤ 1.0) "Pump efficiency must be between 0 and 1."
    pressurized_water = pressurize(feed, pressure_to_apply)
    power_consumption = (feed.Q / 3600) * pressure_to_apply / efficiency * u"W"
    return (pressurized_water, power_consumption)
end

# ===================================================================
# Osmotic Pressure Calculations
# ===================================================================

"""
    osmo_p(concentration::Float64, temperature::Float64)::Float64
    osmo_p(concentration::Float64, temperature::Float64, m_avg::Float64)::Float64
    osmo_p(water::Union{Water, Water2})::Float64

Calculate osmotic pressure using van't Hoff equation.

# Arguments
- `concentration::Float64`: Salt concentration (TDS) [kg/m³]
- `temperature::Float64`: Temperature [°C]
- `m_avg::Float64`: Average ion molecular weight [g/mol] (default: 58.44 for NaCl)
- `water::Union{Water, Water2}`: Water stream (extracts C, T, m_avg)

# Returns
- `Float64`: Osmotic pressure [Pa]

# Formula
π = (i × C × R × T) / M_avg

Where:
- i = 2 (van't Hoff factor for 1:1 salt)
- C = concentration [kg/m³]
- R = 8.3145 kJ/(kmol·K) (gas constant)
- T = temperature [K]
- M_avg = average molecular weight [g/mol]

# Note
Assumes complete dissociation. For seawater or mixed electrolytes,
experimental corrections may be needed at high concentrations.

# Examples
```julia
π = osmo_p(35.0, 25.0)          # Seawater at 25°C: ~27 bar
π = osmo_p(35.0, 25.0, 58.44)   # Explicitly specify NaCl
feed = Water(10.0, 25.0, 35.0, 1e5)
π = osmo_p(feed)                # From Water struct
```
"""
function osmo_p(concentration::Float64, temperature::Float64)::Float64
    return 2 / 58.44 * concentration * 8.3145e3 * (temperature + 273.15) # Assuming all of the ions are Na⁺ ions
end

function osmo_p(concentration::Float64, temperature::Float64, m_avg::Float64)::Float64
    return 2 / m_avg * concentration * 8.3145e3 * (temperature + 273.15)
end

function osmo_p(water::Union{Water, Water2})::Float64
    osmo_p(water.C, water.T, water.m_avg)
end


# ===================================================================
# Fouling Helper Functions
# ===================================================================

"""
    cake_mass(v, u, cf)

Calculate cake layer mass deposition flux on membrane surface.
Uses critical flux concept with empirical correlations.

# Arguments
- `v`: Transmembrane water flux [m/s]
- `u`: Tangential brine velocity [m/s]
- `cf`: Foulant concentration in bulk [kg/m³]

# Returns
- Cake mass deposition flux [kg/m²·s]

# Model Description
Critical flux model where deposition only occurs above threshold:
- v_crit = 15 L/m²/h × (u/0.1 m/s)^0.4
- Deposition flux = cf × max(0, v - v_crit)

Higher crossflow velocity increases critical flux (shear removes foulants).

# Reference
Chong et al. (2008), J. Membrane Sci. 314:89-95
https://doi.org/10.1016/j.memsci.2008.01.030
"""
function cake_mass(v, u, cf)
    v_crit = 15.0 / 3600.0 / 1e3 * (u / 0.1)^0.4
    mf_c = cf * max(0.0, v - v_crit)
    return mf_c
end

"""
    k_CEMT(u, H, W, μ, ρ, ϵ, δ)

Calculate Cake-Enhanced Mass Transfer (CEMT) coefficient.
Accounts for how cake layer affects salt transport to membrane surface.

# Arguments
- `u`: Brine velocity [m/s]
- `H`: Channel height [m]
- `W`: Channel width [m]
- `μ`: Water dynamic viscosity [Pa·s]
- `ρ`: Water density [kg/m³]
- `ϵ`: Cake layer porosity [-]
- `δ`: Cake layer thickness [m]

# Returns
- Cake-enhanced mass transfer coefficient [m/s]

# Model Description
Combines:
1. Bulk mass transfer (Sherwood correlation)
2. Diffusion through porous cake layer
3. Tortuosity effects in cake

Key equations:
- Sh = 0.065 × Re^0.875 × Sc^0.25 (Schock-Miquel correlation)
- τ_c = 1 - 2ln(ϵ) (cake tortuosity)
- 1/k_CEMT = δ(1/D_c - 1/D_b) + 1/k (series resistances)

# Reference
Jiang et al. (2021), Desalination 510:115289
https://doi.org/10.1016/j.desal.2021.115289
"""
function k_CEMT(u, H, W, μ, ρ, ϵ, δ)
    D_b = 1.51e-9                       # Hydraulic dispersion coefficient [m²/s]
    τ_c = 1 - 2log(ϵ)                   # Tortuosity of the cake layer [-]
    D_c = D_b * ϵ / τ_c                 # Diffusivity of solute within cake layer [m²/s]
    ν   = μ / ρ                         # Kinematic viscosity of water [m²/s]
    dh  = 2 * W * H / (W + H)           # Hydraulic diameter of the spacer channel [m]
    Re  = (4 * dh * u) / ν              # Reynolds number [-]
    Sc  = ν / D_b                       # Schmidt number  [-]
    Sh  = 0.065 * Re^0.875 * Sc^0.25    # Sherwood number [-]
    k   = Sh * D_b / dh                 # Mass transfer coefficient [m/s]

    k_CEMT = 1 / (δ * (1 / D_c - 1 / D_b) + 1 / k) # Cake-enhanced mass transfer coefficient [m/s]

    return k_CEMT
end

# ===================================================================
# Element Filtration Algorithms (Core Simulation Functions)
# ===================================================================
# Performance critical: Intentionally avoids Unitful.jl for ~100x speedup
# All inputs/outputs are in SI base units (Pa, m, kg, s, °C)
# ===================================================================

"""
    element_filtration(feed::Water, element::MembraneElement; dt::Union{Float64, Nothing}=nothing)

Simulate single-pass filtration through one membrane element.
Uses resistance model with fixed salt rejection coefficient.

# Arguments
- `feed::Water`: Inlet water stream
- `element::MembraneElement`: Membrane element to simulate
- `dt::Union{Float64, Nothing}=nothing`: Time step [s] (if provided, returns absolute fouling change)

# Returns
- Tuple of `(brine, permeate, ΔR_m)`
  - `brine::Water`: Concentrate stream exiting element
  - `permeate::Water`: Permeate stream produced
  - `ΔR_m::Float64`: Resistance increase [Pa·s/m] or rate [Pa·s/m·s] depending on dt

# Model Description
Iteratively solves for brine concentration considering:
1. **Water flux**: v_w = (P - π) / R_m (resistance model)
2. **Salt transport**: Based on rejection coefficient
3. **Pressure drop**: Hagen-Poiseuille flow through spacer
4. **Concentration polarization**: Iterative mass balance
5. **Fouling**: ΔR_m = k_fp × v_w (× dt if provided)

# Convergence
- Tolerance: 1e-6 relative error
- Max iterations: 100
- Fails with assertion if non-convergent

# Physical Constraints
- Asserts positive brine velocity (no backflow)
- Handles zero permeate flux (when P < π)

# Note
Temperature-dependent viscosity: Vogel equation
μ = 2.414×10⁻⁵ × 10^(247.8/(T+133.15))

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 55e5)
element = MembraneElement(7e-4, 37.0/50, 1.016/50, 16.0, 6.8e10, 0.67e9, 0.995)
brine, permeate, dR = element_filtration(feed, element; dt=3600.0)
```
"""
function element_filtration(
    feed::Water, element::MembraneElement; dt::Union{Float64, Nothing}=nothing
)::Tuple{Water, Water, Float64}
    C_feed = feed.C
    Q_feed = feed.Q
    P_feed = feed.P
    T_feed = feed.T
    M_feed = feed.m_avg
    # Membrane cross-section flux in m/s
    U_feed = Q_feed / element.width / element.height / 3600

    K      = element.spacer_resistance      # Spacer resistance coefficient
    k_fp   = element.k_fp                   # Fouling potential coefficient
    reject = element.salt_rejection         # Membrane salt rejection rate
    μ      = 2.414e-5 * 10^(247.8 / (T_feed + 273.15 - 140)) # Water viscosity coefficient
    
    local err         ::Float64 = 1.0
    local idx_iter    ::Int64   = 0
    local c_cal       ::Float64 = C_feed
    local v_w_guess   ::Float64 = 0.0
    local u_guess     ::Float64 = 0.0
    local osmo_p_guess::Float64 = 0.0

    # Main loop
    while err > 1e-6
        c_guess = c_cal
        osmo_p_guess = osmo_p(c_guess, T_feed, M_feed)
        v_w_guess    = max(P_feed - osmo_p_guess, 0.0) / element.R_m
        u_guess      = U_feed - v_w_guess * element.dx / element.height
        c_cal        = ((C_feed * U_feed * element.height) - ((1-reject)*C_feed) * (v_w_guess * element.dx)) / (u_guess * element.height)
        err          = abs((c_cal - c_guess) / c_cal)
        idx_iter     += 1
        if idx_iter ≥ 100
            @assert false "\nFailed to converge. Iteration exceeded 100 times. Feed water: $(feed) osmo_p_guess: $(osmo_p_guess) v_w_guess: $(v_w_guess)"
        end
    end
    
    # @assert (osmo_p_guess - P_feed < 0.0) "Transmembrane pressure is negative (TMP: $(P_feed - osmo_p_guess)). Forward osmosis isn't currently supported."
    # @assert (v_w_guess > 0.0) "Permeate membrane flux is negative (v_w_guess: $(v_w_guess)). Forward osmosis isn't currently supported.\nFeed: $(feed)"
    @assert (u_guess   > 0.0) "Brine flux is negative (u_guess: $(u_guess)). Backward flow isn't currently supported."

    Q_p = v_w_guess * element.width * element.dx * 3600 # Permeate flowrate [m³/h]
    Q_b = Q_feed - Q_p                                  # Brine flowrate    [m³/h]
    C_p = (1-reject)*C_feed
    C_b = c_cal
    P_p = 1e-10
    P_b = P_feed - ((12 * K * μ * u_guess * element.dx) / element.height^2)

    # If dt isn't given, return ΔR_m which is d(R_m)/dt.
    # If dt is given, return ΔR_m which is actual increase of membrane resistance.
    if isnothing(dt)
        ΔR_m = k_fp * v_w_guess
    else
        ΔR_m = k_fp * v_w_guess * dt
    end

    permeate = Water(Q_p, T_feed, C_p, P_p)
    brine    = Water(Q_b, T_feed, C_b, P_b)

    return (brine, permeate, ΔR_m)
end

"""
    element_filtration2(feed::Water2, element::MembraneElement2; dt::Union{Float64, Nothing}=nothing)

Advanced element filtration with solution-diffusion model and cake layer dynamics.
Tracks foulant concentration and cake layer formation.

# Arguments
- `feed::Water2`: Inlet water stream with foulant concentration
- `element::MembraneElement2`: Membrane element with cake layer tracking
- `dt::Union{Float64, Nothing}=nothing`: Time step [s] (if provided, returns absolute changes)

# Returns
- Tuple of `(brine, permeate, ΔCake_thickness, ΔR_c)`
  - `brine::Water2`: Concentrate stream exiting element
  - `permeate::Water2`: Permeate stream produced
  - `ΔCake_thickness::Float64`: Cake thickness increase [m] or rate [m/s]
  - `ΔR_c::Float64`: Cake resistance increase [m⁻¹] or rate [m⁻¹/s]

# Model Description
**Solution-Diffusion Model:**
- Water flux: v_w = A × (P - π) where A = 1/(R_m + R_c)/μ
- Salt flux: v_s = B × C_membrane

**Concentration Polarization:**
- Film theory with cake enhancement
- C_membrane = C_bulk × exp(v_w / k_CEMT)
- k_CEMT accounts for mass transfer through cake layer

**Fouling Dynamics:**
- Critical flux model for cake deposition
- Carman-Kozeny equation for cake resistance:
  R_c = 180(1-ε)/ρ_p/d_p²/ε³ × m_cake
- Foulants completely rejected (Cf_permeate = 0)

**Temperature Correction:**
- TCF_A = exp(4140(1/T - 1/20°C)) for water permeability
- TCF_B = exp(-4968(1/T - 1/20°C)) for salt permeability

# Physical Parameters (Hard-coded)
- Cake porosity: ε_p = 0.36
- Particle density: ρ_p = 1400 kg/m³
- Particle diameter: d_p = 20 nm
- Critical flux: 15 L/m²/h @ 0.1 m/s

# Convergence
- Tolerance: 1e-4 relative error (looser than element_filtration)
- Max iterations: 100
- Solves for membrane surface concentration

# Physical Constraints
- Asserts positive brine velocity
- Prevents negative fluxes
- Limits foulant deposition to available mass

# Note
Cake layer reduces effective channel height, changing local hydraulics.

# Examples
```julia
feed = Water2(10.0, 25.0, 35.0, 0.1, 55e5)  # 0.1 kg/m³ foulant
element = MembraneElement2(7e-4, 37.0/50, 1.016/50, 16.0, 1.5e15, 0.0, 3.5e-8, 0.0, 0.8)
brine, permeate, dThick, dRc = element_filtration2(feed, element; dt=3600.0)
```
"""
function element_filtration2(
    feed::Water2, element::MembraneElement2; dt::Union{Float64, Nothing}=nothing
)::Tuple{Water2, Water2, Float64, Float64}
    C_feed  = feed.C
    Cf_feed = feed.Cf
    Q_feed  = feed.Q
    P_feed  = feed.P
    T_feed  = feed.T
    M_feed  = feed.m_avg

    # Membrane local height (Pristine height - cake thickness)
    local_height = element.height - element.cake_thickness
    # Membrane cross-section flux in m/s
    U_feed = Q_feed / element.width / local_height / 3600
    
    K      = element.spacer_resistance      # Spacer resistance coefficient
    k_fp   = element.k_fp                   # Fouling potential coefficient
    
    μ      = 2.414e-5 * 10^(247.8 / (T_feed + 273.15 - 140)) # Water viscosity coefficient [Pa·s]
    ρ      = 1e3                            # Water density (No temperature correction) [kg/m³]
    
    ϵₚ     = 0.36                           # Cake layer porosity [-]
    ρₚ     = 1.4e3                          # Cake layer particle density  [kg/m³]
    dₚ     = 20e-9                          # Cake layer particle diameter [m]
    
    A = 1 / (element.R_m + element.R_c) / μ  # Water permeability [m³/(m²·Pa·s)]
    B      = element.B                       # Salt permeability  [kg/(m²·Pa·s)]
    TCF_A  = exp(4140.0  * (1/(T_feed + 273.15) - 1/293.15)) # Temperature correction coefficient [-]
    TCF_B  = exp(-4968.0 * (1/(T_feed + 273.15) - 1/293.15)) # Temperature correction coefficient [-]
    A_corr = A / TCF_A                       # Corrected water permeability (Divided because TCF_A is applied to A⁻¹(Resistance))
    B_corr = B * TCF_B                       # Corrected salt permeability
    
    # All of the concentration units are in kg/m³
    # All of the flux units are in m/s
    # All of the pressure units are in Pa

    local err         ::Float64 = 1.0       # Error
    local idx_iter    ::Int64   = 0         # Iteration index
    local cm_cal      ::Float64 = C_feed    # Brine concentration
    local cm_guess    ::Float64 = 0.0
    local c_guess     ::Float64 = 0.0       # Bulk concentration
    local cp_guess    ::Float64 = 0.0       # Permeate concentration
    local v_w_guess   ::Float64 = 0.0       # Water flux through the membrane
    local v_s_guess   ::Float64 = 0.0       # Salt flux through the membrane
    local u_guess     ::Float64 = 0.0       # Brine flux over the membrane
    local osmo_p_guess::Float64 = 0.0       # Osmotic pressure
    
    # Main loop
    while err > 1e-4
        cm_guess     = cm_cal
        osmo_p_guess = osmo_p(cm_guess, T_feed, M_feed)
        v_w_guess    = max(P_feed - osmo_p_guess, 1e-10) * A_corr
        v_s_guess    = B_corr * cm_guess
        u_guess      = U_feed - v_w_guess * element.dx / local_height
        c_guess      = (U_feed * C_feed * local_height - v_s_guess * element.dx) / local_height / u_guess
        k_mt         = k_CEMT(u_guess, local_height, element.width, μ, ρ, ϵₚ, element.cake_thickness)
        cm_cal       = c_guess * exp(v_w_guess / k_mt)
        err          = abs((cm_cal - cm_guess) / cm_cal)
        idx_iter     += 1
        if idx_iter ≥ 100
            @assert false "\nFailed to converge. Iteration exceeded 100 times. Feed water: $(feed) | err: $(err)"
        end
    end
    
    @assert (u_guess > 0.0) "Brine flux is negative (u_guess: $(u_guess)). Backward flow isn't supported."
    cp_guess = v_s_guess / v_w_guess

    m_fc_in  = Cf_feed * U_feed * local_height                           # Foulant mass inflow (kg/m⋅s)
    m_fc_acc = min(
        k_fp * cake_mass(v_w_guess, u_guess, cm_cal) * element.dx,
        m_fc_in
        )                                                                # Foulant mass deposition rate (kg/m⋅s)
    m_fc_out = m_fc_in - m_fc_acc                                        # Foulant mass outflow (kg/m⋅s)
    cf_guess = m_fc_out / (u_guess * local_height)                       # Foulant concentration (kg/m³)
    
    

    Q_p  = v_w_guess * element.width * element.dx * 3600.0 # Permeate flowrate [m³/h]
    Q_b  = Q_feed - Q_p                                    # Brine flowrate    [m³/h]
    C_p  = cp_guess
    C_b  = c_guess
    # C_b  = cm_guess
    P_p  = 1e5                                             # Atmospheric pressure
    P_b  = P_feed - ((12 * K * μ * u_guess * element.dx) / local_height^2)
    Cf_p = 0.0                                             # Assume complete filtration of foulant (which is likely!)
    Cf_b = cf_guess

    # If dt isn't given, return ΔR_m which is d(R_m)/dt.
    # If dt is given, return ΔR_m which is actual increase of membrane resistance.
    if isnothing(dt)
        ΔCake_thickness = m_fc_acc / ρₚ / (1 - ϵₚ) / element.dx                   # Cake thickness increase rate (m/s)
        ΔR_c            = 180.0 * (1 - ϵₚ) / ρₚ / dₚ^2 / ϵₚ^3 * m_fc_acc           # Cake layer resistance increase rate (Carman-Kozeny equation in Darcy's law)
    else
        ΔCake_thickness = m_fc_acc / ρₚ / (1 - ϵₚ) / element.dx * dt              # Cake thickness increase (m)
        ΔR_c            = 180.0 * (1 - ϵₚ) / ρₚ / dₚ^2 / ϵₚ^3 * m_fc_acc * dt      # Cake layer resistance increase (Carman-Kozeny equation in Darcy's law)
    end

    permeate = Water2(Q_p, T_feed, C_p, Cf_p, P_p)
    brine    = Water2(Q_b, T_feed, C_b, Cf_b, P_b)

    return (brine, permeate, ΔCake_thickness, ΔR_c)
end

# ===================================================================
# Module and Vessel Filtration (Cascaded Elements)
# ===================================================================
# These functions call element filtration in sequence, with brine from
# one element becoming feed for the next (series flow arrangement)
# ===================================================================

"""
    module_filtration(feed::Water, membrane::MembraneModule; dt::Union{Float64, Nothing}=nothing)

Simulate filtration through entire membrane module (cascade of elements).
Automatically dispatches to element_filtration based on types.

# Arguments
- `feed::Water`: Inlet water stream
- `membrane::MembraneModule`: Membrane module containing element array
- `dt::Union{Float64, Nothing}=nothing`: Time step for fouling integration [s]

# Returns
- Tuple of `(brines, permeates, ΔR_ms)`
  - `brines::Array{Water}`: Brine stream from each element
  - `permeates::Array{Water}`: Permeate stream from each element
  - `ΔR_ms::Array{Float64}`: Resistance increase for each element

# Description
Simulates sequential filtration where:
1. Feed enters first element
2. Brine from element_i becomes feed for element_(i+1)
3. Permeates collected separately from each element
4. Final brine = brines[end]

# Type Compatibility
- Requires Water feed with MembraneElement array
- Type mismatch triggers assertion error

# Note
Arrays length = membrane.n_segments
Useful for analyzing spatial profiles along module length

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 55e5)
membrane = pristine_membrane(50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
brines, permeates, dRs = module_filtration(feed, membrane; dt=3600.0)
# brines has 50 elements, one per membrane segment
```
"""
function module_filtration(
    feed::Water, membrane::MembraneModule; dt::Union{Float64, Nothing}=nothing
)::Tuple{Array{Water}, Array{Water}, Array{Float64}}
    local next_feed::Water
    local brine::Water
    local permeate::Water
    local ΔR_m::Float64

    @assert (typeof(membrane.elements_array[1]) <: MembraneElement)
            "Water type and membrane type do not match (Element type $(typeof(membrane.elements_array[1])), Feed type $(typeof(feed)))." 
    
    brines    = []
    permeates = []
    ΔR_ms     = []
    
    next_feed = feed
    for (i, element) in enumerate(membrane.elements_array)
        # @debug "Processing segment #$(i)"
        brine, permeate, ΔR_m = element_filtration(next_feed, element; dt=dt)
        push!(brines, brine)
        push!(permeates, permeate)
        push!(ΔR_ms, ΔR_m)
        next_feed = brine
    end

    return brines, permeates, ΔR_ms
end

"""
    module_filtration(feed::Water2, membrane::MembraneModule; dt::Union{Float64, Nothing}=nothing)

Simulate module filtration with foulant tracking and cake layer formation.
Uses element_filtration2 for advanced fouling dynamics.

# Returns
- Tuple of `(brines, permeates, ΔCake_thicknesses, ΔR_cs)`
  - Returns cake thickness and resistance arrays instead of simple ΔR_m

# Note
Requires Water2 feed and MembraneElement2 array for type compatibility.
"""
function module_filtration(
    feed::Water2, membrane::MembraneModule; dt::Union{Float64, Nothing}=nothing
)::Tuple{Array{Water2}, Array{Water2}, Array{Float64}, Array{Float64}}
    local next_feed::Water2
    local brine::Water2
    local permeate::Water2
    local ΔCake_thickness::Float64
    local ΔR_c::Float64

    @assert (typeof(membrane.elements_array[1]) <: MembraneElement2) & (typeof(feed) <: Water2)
            "Water type and membrane type do not match (Element type $(typeof(membrane.elements_array[1])), Feed type $(typeof(feed)))." 
    
    brines    = []
    permeates = []
    ΔCake_thicknesses = []
    ΔR_cs     = []
    
    next_feed = feed
    for (i, element) in enumerate(membrane.elements_array)
        # @info "\tProcessing segment #$(i)"
        brine, permeate, ΔCake_thickness, ΔR_c = element_filtration2(next_feed, element; dt=dt)
        push!(brines, brine)
        push!(permeates, permeate)
        push!(ΔCake_thicknesses, ΔCake_thickness)
        push!(ΔR_cs, ΔR_c)
        next_feed = brine
    end

    # @info "Cake thickness increase (w.o. dt)" ΔCake_thicknesses

    return brines, permeates, ΔCake_thicknesses, ΔR_cs
end

"""
    vessel_filtration(feed::Water, vessel::PressureVessel; dt::Union{Float64, Nothing}=nothing)

Simulate filtration through entire pressure vessel (multiple modules in series).

# Arguments
- `feed::Water`: Inlet water stream
- `vessel::PressureVessel`: Pressure vessel containing module array
- `dt::Union{Float64, Nothing}=nothing`: Time step for fouling integration [s]

# Returns
- Tuple of `(brines_array, permeates_array, ΔR_ms_array)`
  - `brines_array::Array{Array{Water}}`: Nested array [module][element] of brines
  - `permeates_array::Vector{Vector{Water}}`: Nested array [module][element] of permeates
  - `ΔR_ms_array::Vector{Vector{Float64}}`: Nested array [module][element] of resistance changes

# Description
Cascades modules where brine from last element of module_i
becomes feed for first element of module_(i+1).

# Structure
- Outer array: modules (length = n_modules)
- Inner array: elements (length = n_segments per module)
- Total elements simulated = sum(n_segments)

# Common Usage
```julia
vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
feed = Water(10.0, 25.0, 35.0, 55e5)
brines, permeates, dRs = vessel_filtration(feed, vessel; dt=3600.0)

# Final outputs
final_brine = brines[end][end]  # Last element of last module
total_permeate = mix(reduce(vcat, permeates))  # Mix all permeates
```
"""
function vessel_filtration(
    feed::Water, vessel::PressureVessel; dt::Union{Float64, Nothing}=nothing
)::Tuple{Array{Array{Water}}, Vector{Vector{Water}}, Vector{Vector{Float64}}}
    local next_feed::Water
    local brines::Array{Water}
    local permeates::Array{Water}
    local ΔR_ms::Array{Float64}

    brines_array    = []
    permeates_array = []
    ΔR_ms_array     = []

    next_feed = feed
    for membrane in vessel.modules_array
        brines, permeates, ΔR_ms = module_filtration(next_feed, membrane; dt=dt)
        push!(brines_array, brines)
        push!(permeates_array, permeates)
        push!(ΔR_ms_array, ΔR_ms)
        next_feed = brines[end]
    end

    return brines_array, permeates_array, ΔR_ms_array
end

"""
    vessel_filtration(feed::Water2, vessel::PressureVessel; dt::Union{Float64, Nothing}=nothing)

Vessel filtration with foulant tracking and cake layer formation.
Returns nested arrays for brines, permeates, cake thickness changes, and resistance changes.

# Note
Requires Water2 feed and MembraneElement2 array in vessel for type compatibility.
"""
function vessel_filtration(
    feed::Water2, vessel::PressureVessel; dt::Union{Float64, Nothing}=nothing
)::Tuple{Array{Array{Water2}}, Array{Array{Water2}}, Vector{Vector{Float64}}, Vector{Vector{Float64}}}
    local next_feed::Water2
    local brines::Array{Water2}
    local permeates::Array{Water2}
    local ΔCake_thicknesses::Array{Float64}
    local ΔR_cs::Array{Float64}

    brines_array          = []
    permeates_array       = []
    ΔCake_thicknesses_array = []
    ΔR_cs_array           = []

    next_feed = feed
    for (i, membrane) in enumerate(vessel.modules_array)
        # @info "Processing module # $(i)"
        # @info "This feed: " next_feed
        brines, permeates, ΔCake_thicknesses, ΔR_cs = module_filtration(next_feed, membrane; dt=dt)
        push!(brines_array, brines)
        push!(permeates_array, permeates)
        push!(ΔCake_thicknesses_array, ΔCake_thicknesses)
        push!(ΔR_cs_array, ΔR_cs)
        next_feed = brines[end]
        # @info "Last brine:" brines[end]
    end

    return brines_array, permeates_array, ΔCake_thicknesses_array, ΔR_cs_array
end

end

