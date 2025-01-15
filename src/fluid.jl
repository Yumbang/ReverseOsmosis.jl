module Fluid
using Printf, DataFrames

export Water, mix, pressurize, profile_water
"""
`Water` struct describes water in reverse osmosis process.
"""
struct Water
    Q::Float64      # Flowrate of water             [m³/h]
    T::Float64      # Temperature of water          [°C]
    C::Float64      # Concentration of water as TDS [kg/m³]
    P::Float64      # Pressure of water             [Pa]
end

"""
Mix two `Water` instances with `pressure` keyword argument
Returns Mixed `Water` instance.
"""
function mix(water1::Water, water2::Water; pressure::Float64)
    flow            = water1.Q + water2.Q
    temperature     = (water1.T * water1.Q + water2.T * water2.Q) / flow
    concentration   = (water1.C * water1.Q + water2.C * water2.Q) / flow
    pressure        = pressure
    return Water(flow, temperature, concentration, pressure)
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

function profile_water(water::Water)
    # Gather field names from the struct type
    names = fieldnames(Water)
    named_tuple = NamedTuple{names}(getfield(water, f) for f in names)
    return DataFrame(named_tuple)
end

# function profile_water(water_array::Array{Water})
#     return DataFrame
# end

end