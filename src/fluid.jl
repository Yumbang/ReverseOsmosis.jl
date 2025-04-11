"""
`Water` struct describes water in reverse osmosis process.
"""
struct Water
    Q::Float64      # Flowrate of water             [m³/h]
    T::Float64      # Temperature of water          [°C]
    C::Float64      # Concentration of water as TDS [kg/m³]
    P::Float64      # Pressure of water             [Pa]
    m_avg::Float64  # Average ion molecular weight  [g/mol]
end

struct Water2
    Q::Float64      # Flowrate of water             [m³/h]
    T::Float64      # Temperature of water          [°C]
    C::Float64      # Salt concentration of water as TDS [kg/m³]
    Cf::Float64     # Foulant concentration of water as TDS [kg/m³]
    P::Float64      # Pressure of water             [Pa]
    m_avg::Float64  # Average ion molecular weight  [g/mol]
end

function Water(Q, T, C, P)
    Water(Q, T, C, P, 58.44) # Assume most of the ions are NaCl
end

function Water2(Q, T, C, Cf, P)
    Water2(Q, T, C, Cf, P, 58.44) # Assume most of the ions are NaCl
end

"""
Mix `Water` instances with `pressure` keyword argument
Returns Mixed `Water` instance.
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

function mix(water_array::Union{Array{Water}, Array{Water2}}; pressure::Float64=1e5)
    return reduce(((water1, water2) -> mix(water1, water2; pressure)), water_array)
end


"""
Pressurize an `Water(or 2)` with pressure_to_apply.
Returns new `Water(or 2)` by adding up the pressure.
"""
function pressurize(water::Water, pressure_to_apply::Float64)
    return Water(water.Q, water.T, water.C, water.P + pressure_to_apply)
end

function pressurize(water::Water2, pressure_to_apply::Float64)
    return Water2(water.Q, water.T, water.C, water.Cf, water.P + pressure_to_apply)
end

function Base.show(io::IO, water::Water)
    @printf io "Q: %.2f C: %.2f P: %.2f T: %.2f" water.Q water.C water.P water.T
end

function Base.show(io::IO, water::Water2)
    @printf io "Q: %.2f C: %.2f Cf: %.2f P: %.2f T: %.2f" water.Q water.C water.Cf water.P water.T
end


"""
Apply units and return a tuple of the water properties.
"""
function unit_water(water::Water)
    return (Q = water.Q*u"m^3/hr", T = water.T*u"°C", C = water.C*u"kg/m^3", P = water.P*u"Pa", M_avg =  water.m_avg*u"g/mol")
end

function unit_water(water::Water2)
    return (Q = water.Q*u"m^3/hr", T = water.T*u"°C", C = water.C*u"kg/m^3", Cf = water.Cf*u"kg/m^3", P = water.P*u"Pa", M_avg =  water.m_avg*u"g/mol")
end

"""
Create a *Unitful* DataFrame with Water or Water2 properties.
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