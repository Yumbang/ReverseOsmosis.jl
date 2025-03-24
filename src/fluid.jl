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

function Water(Q, T, C, P)
    Water(Q, T, C, P, 58.44) # Assume most of the ions are NaCl
end

"""
Mix `Water` instances with `pressure` keyword argument
Returns Mixed `Water` instance.
"""
function mix(water1::Water, water2::Water; pressure::Float64=1e-5)
    flow            = water1.Q + water2.Q
    temperature     = (water1.T * water1.Q + water2.T * water2.Q) / flow
    concentration   = (water1.C * water1.Q + water2.C * water2.Q) / flow
    m_avg           = (water1.C * water1.Q + water2.C * water2.Q) / (water1.C * water1.Q / water1.m_avg + water2.C * water2.Q / water2.m_avg)
    return Water(flow, temperature, concentration, pressure, m_avg)
end

function mix(water_array::Array{Water}; pressure::Float64)
    return reduce(((water1, water2) -> mix(water1, water2; pressure)), water_array)
end

"""
Pressurize an `Water` with pressure_to_apply.
Returns new `Water` by adding up the pressure.
"""
function pressurize(water::Water, pressure_to_apply::Float64)
    return Water(water.Q, water.T, water.C, water.P + pressure_to_apply)
end

function Base.show(io::IO, water::Water)
    @printf io "Q: %.2f C: %.2f P: %.2f T: %.2f" water.Q water.C water.P water.T
end

function unit_water(water::Water)
    return (Q = water.Q*u"m^3/hr", T = water.T*u"°C", C = water.C*u"kg/m^3", P = water.P*u"Pa")
end

function profile_water(water::Water)
    water_tuple = unit_water(water)
    return DataFrame([water_tuple])
end

function profile_water(waters::Array{Water})
    # Gather field names from the struct type
    water_tuples = map(unit_water, waters)
    return DataFrame(water_tuples)
end