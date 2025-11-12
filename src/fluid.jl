# ===================================================================
# Water Stream Representations
# ===================================================================

"""
    Water(Q, T, C, P, m_avg)
    Water(Q, T, C, P)  # Assumes NaCl (m_avg = 58.44 g/mol)

Immutable structure representing a water stream in reverse osmosis processes.

# Fields
- `Q::Float64`: Volumetric flowrate [m³/h]
- `T::Float64`: Temperature [°C]
- `C::Float64`: Total dissolved solids (TDS) concentration [kg/m³]
- `P::Float64`: Pressure [Pa]
- `m_avg::Float64`: Average ion molecular weight [g/mol]

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 1e5)  # 10 m³/h, 25°C, 35 kg/m³ TDS, 1 bar
```
"""
struct Water
    Q::Float64      # Flowrate of water             [m³/h]
    T::Float64      # Temperature of water          [°C]
    C::Float64      # Concentration of water as TDS [kg/m³]
    P::Float64      # Pressure of water             [Pa]
    m_avg::Float64  # Average ion molecular weight  [g/mol]
end

"""
    Water2(Q, T, C, Cf, P, m_avg)
    Water2(Q, T, C, Cf, P)  # Assumes NaCl (m_avg = 58.44 g/mol)

Extended water stream structure with foulant tracking capability.
Used for simulations that model membrane fouling dynamics.

# Fields
- `Q::Float64`: Volumetric flowrate [m³/h]
- `T::Float64`: Temperature [°C]
- `C::Float64`: Salt concentration (TDS) [kg/m³]
- `Cf::Float64`: Foulant particle concentration [kg/m³]
- `P::Float64`: Pressure [Pa]
- `m_avg::Float64`: Average ion molecular weight [g/mol]

# Examples
```julia
feed = Water2(10.0, 25.0, 35.0, 0.1, 1e5)  # 0.1 kg/m³ foulant concentration
```
"""
struct Water2
    Q::Float64      # Flowrate of water             [m³/h]
    T::Float64      # Temperature of water          [°C]
    C::Float64      # Salt concentration of water as TDS [kg/m³]
    Cf::Float64     # Foulant concentration of water as TDS [kg/m³]
    P::Float64      # Pressure of water             [Pa]
    m_avg::Float64  # Average ion molecular weight  [g/mol]
end

# Convenience constructors assuming NaCl as primary salt
function Water(Q, T, C, P)
    Water(Q, T, C, P, 58.44) # Assume most of the ions are NaCl
end

function Water2(Q, T, C, Cf, P)
    Water2(Q, T, C, Cf, P, 58.44) # Assume most of the ions are NaCl
end

# ===================================================================
# Stream Mixing Operations
# ===================================================================

"""
    mix(water1::Water, water2::Water; pressure::Float64=1e5)

Mix two water streams using mass-weighted averaging.
Returns a new `Water` instance with combined properties.

# Arguments
- `water1::Water`: First water stream
- `water2::Water`: Second water stream
- `pressure::Float64=1e5`: Pressure of mixed stream [Pa] (default: 1 bar)

# Returns
- `Water`: Mixed water stream with flow-weighted average properties

# Note
Temperature, concentration, and molecular weight are calculated using
mass-weighted averaging to conserve mass and energy.

# Examples
```julia
stream1 = Water(5.0, 25.0, 35.0, 1e5)
stream2 = Water(3.0, 20.0, 40.0, 1e5)
mixed = mix(stream1, stream2; pressure=1e5)  # 8 m³/h combined flow
```
"""
function mix(water1::Water, water2::Water; pressure::Float64=1e5)
    if ((water1.Q == 0.0) & (water2.Q == 0.0))
        @warn "Flow rate of both Water instances are zero. Mixing result will not be reliable."
    end
    flow            = water1.Q + water2.Q
    temperature     = (water1.T * water1.Q + water2.T * water2.Q) / flow
    concentration   = (water1.C * water1.Q + water2.C * water2.Q) / flow
    m_avg           = (water1.C * water1.Q + water2.C * water2.Q) / (water1.C * water1.Q / water1.m_avg + water2.C * water2.Q / water2.m_avg)
    return Water(flow, temperature, concentration, pressure, m_avg)
end

"""
    mix(water1::Water2, water2::Water2; pressure::Float64=1e5)

Mix two Water2 streams including foulant concentration tracking.
"""
function mix(water1::Water2, water2::Water2; pressure::Float64=1e5)
    if (water1.Q == 0.0) & (water2.Q == 0.0)
        @warn "Flow rate of both Water2 instances are zero. Mixing result will not be reliable."
    end
    flow            = water1.Q + water2.Q
    temperature     = (water1.T * water1.Q + water2.T * water2.Q)   / flow
    concentration   = (water1.C * water1.Q + water2.C * water2.Q)   / flow
    Cf              = (water1.Cf * water1.Q + water2.Cf * water2.Q) / flow
    m_avg           = (water1.C * water1.Q + water2.C * water2.Q)   / (water1.C * water1.Q / water1.m_avg + water2.C * water2.Q / water2.m_avg)
    return Water2(flow, temperature, concentration, Cf, pressure, m_avg)
end

"""
    mix(water_array::Union{Array{Water}, Array{Water2}}; pressure::Float64=1e5)

Mix multiple water streams by reducing pairwise mixing operations.
Useful for combining permeate streams from multiple membrane elements.

# Examples
```julia
permeate_streams = [Water(1.0, 25.0, 0.5, 1e5) for _ in 1:10]
total_permeate = mix(permeate_streams)
```
"""
function mix(water_array::Union{Array{Water}, Array{Water2}}; pressure::Float64=1e5)
    return reduce(((water1, water2) -> mix(water1, water2; pressure)), water_array)
end


# ===================================================================
# Pressure Operations
# ===================================================================

"""
    pressurize(water::Water, pressure_to_apply::Float64)

Apply pressure to a water stream (e.g., via pump).
Returns a new `Water` instance with increased pressure.

# Arguments
- `water::Water`: Input water stream
- `pressure_to_apply::Float64`: Pressure increase [Pa]

# Returns
- `Water`: New water stream with pressure = water.P + pressure_to_apply

# Note
This function does not calculate energy consumption. For pump modeling
with efficiency considerations, use the `pump` function from Filtration module.

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 1e5)
pressurized = pressurize(feed, 50e5)  # Add 50 bar
```
"""
function pressurize(water::Water, pressure_to_apply::Float64)
    return Water(water.Q, water.T, water.C, water.P + pressure_to_apply)
end

"""
    pressurize(water::Water2, pressure_to_apply::Float64)

Pressurize Water2 stream (includes foulant concentration).
"""
function pressurize(water::Water2, pressure_to_apply::Float64)
    return Water2(water.Q, water.T, water.C, water.Cf, water.P + pressure_to_apply)
end

# ===================================================================
# Display Methods
# ===================================================================

"""
    Base.show(io::IO, water::Water)

Custom display format for Water streams showing key properties.
"""
function Base.show(io::IO, water::Water)
    @printf io "Q: %.2f C: %.2f P: %.2f T: %.2f" water.Q water.C water.P water.T
end

"""
    Base.show(io::IO, water::Water2)

Custom display format for Water2 streams including foulant concentration.
"""
function Base.show(io::IO, water::Water2)
    @printf io "Q: %.2f C: %.2f Cf: %.2f P: %.2f T: %.2f" water.Q water.C water.Cf water.P water.T
end


# ===================================================================
# Profiling and Data Export
# ===================================================================

"""
    unit_water(water::Water)

Convert Water properties to a NamedTuple with proper Unitful units.
Internal helper function for DataFrame creation.

# Returns
- NamedTuple with fields: Q, T, C, P, M_avg (all with appropriate units)
"""
function unit_water(water::Water)
    return (Q = water.Q*u"m^3/hr", T = water.T*u"°C", C = water.C*u"kg/m^3", P = water.P*u"Pa", M_avg =  water.m_avg*u"g/mol")
end

"""
    unit_water(water::Water2)

Convert Water2 properties to a NamedTuple with proper Unitful units.
"""
function unit_water(water::Water2)
    return (Q = water.Q*u"m^3/hr", T = water.T*u"°C", C = water.C*u"kg/m^3", Cf = water.Cf*u"kg/m^3", P = water.P*u"Pa", M_avg =  water.m_avg*u"g/mol")
end

"""
    profile_water(water::Union{Water, Water2})
    profile_water(waters::Union{Array{Water}, Array{Water2}})

Create a Unitful DataFrame from water stream(s) for analysis and visualization.

# Arguments
- `water`: Single water stream or array of water streams

# Returns
- `DataFrame`: DataFrame with columns for Q, T, C, P, M_avg (and Cf for Water2)
                All values include proper Unitful units

# Examples
```julia
feed = Water(10.0, 25.0, 35.0, 1e5)
df = profile_water(feed)

# For multiple streams
permeate_array = [Water(1.0, 25.0, 0.5, 1e5) for _ in 1:10]
profile = profile_water(permeate_array)
```
"""
function profile_water(water::Union{Water, Water2})
    water_tuple = unit_water(water)
    return DataFrame([water_tuple])
end

function profile_water(waters::Union{Array{Water}, Array{Water2}})
    # Gather field names from the struct type
    water_tuples = map(unit_water, waters)
    return DataFrame(water_tuples)
end