module Filtration

export
    # Fluid
    Water, pressurize, profile_water, mix,
    # Membrane
    MembraneModule, MembraneElement, pristine_membrane,
    # Filtration
    osmo_p, element_filtration, pump

include("fluid.jl")
include("membrane.jl")

using .Fluid: Water, pressurize, profile_water, mix
using .Membrane: MembraneElement, MembraneModule, profile_membrane, pristine_membrane


"""
Calculates osmotic pressure of water based on van't Hoff equation.
"""
function osmo_p(concentration::Float64, temperature::Float64)::Float64
    return 2 / 58.44e3 * concentration * 1e3 * 8.3145e3 * (temperature + 273.15)
end

function osmo_p(water::Water)::Float64
    osmo_p(water.C, water.T)
end

"""
Pressurize an `Water` with pressure to apply [Pa], with keyword argument `efficiency`.
Returns pressurized `Water` and consumed power [W] as following:
Power = (Feed volume⋅Aplied pressure)/(efficiency)
"""
function pump(feed::Water, pressure_to_apply::Float64; efficiency::Float64=0.8)
    pressurized_water = pressurize(feed, pressure_to_apply)
    power_consumption = (feed.Q / 3600) * pressure_to_apply / efficiency
    return (pressurized_water, power_consumption)
end

"""
Model reverse osmosis of feed water on single membrane element.
Returns brine & permeate and resistance increase (i.e. fouling effect).
"""
function element_filtration(
    feed::Water, element::MembraneElement;
    K_setpoint=16.0, salt_rejection_setpoint=0.995, dt=nothing
)
    C_feed = feed.C
    Q_feed = feed.Q
    P_feed = feed.P
    T_feed = feed.T
    # Membrane cross-section flux in m/s
    U_feed = Q_feed / element.width / element.height / 3600

    K      = K_setpoint     # Spacer resistance coefficient
    k_fp   = element.k_fp   # Fouling potential coefficient
    μ      = 2.414e-5 * 10^(247.8 / (T_feed + 273.15 - 140)) # Water viscosity coefficient
    reject = salt_rejection_setpoint
    
    local err::Float64          = 1.0
    local idx_iter::Int64       = 0
    local c_cal::Float64        = C_feed
    local v_w_guess::Float64    = 0.0
    local u_guess::Float64      = 0.0
    local osmo_p_guess::Float64 = 0.0

    # Main loop
    while err > 1e-6
        c_guess = c_cal
        osmo_p_guess = osmo_p(c_guess, T_feed)
        @assert (osmo_p_guess - P_feed < 0.0) "Transmembrane pressure is negative (TMP: $(P_feed - osmo_p_guess)). Forward osmosis isn't currently supported."
        v_w_guess    = (P_feed - osmo_p_guess) / element.R_m
        u_guess      = U_feed - v_w_guess * element.dx / element.height
        c_cal        = ((C_feed * U_feed * element.height) - ((1-reject)*C_feed) * (v_w_guess * element.dx)) / (u_guess * element.height)
        err          = abs((c_cal - c_guess) / c_cal)
        idx_iter     += 1
        if idx_iter ≥ 100
            @error "Failed to converge. Iteration exceeded 100 times."
        end
    end

    @assert (v_w_guess > 0.0) "Cross membrane flux is negative (v_w_guess: $(v_w_guess)). Forward osmosis isn't currently supported."
    @assert (u_guess > 0.0) "Brine flux is negative (u_guess: $(u_guess)). Backward flow isn't currently supported."

    Q_p = v_w_guess * element.width * element.dx * 3600 # Permeate flowrate [m³/h]
    Q_b = Q_feed - Q_p                                      # Brine flowrate    [m³/h]
    C_p = (1-reject)*C_feed
    C_b = c_cal
    P_p = 1e-10
    P_b = P_feed - ((12 * K * μ * u_guess * element.dx) / element.height^2)

    # If dt isn't given, return ΔR_m which is d(R_m)/dt that will constitute RO ODE problem.
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

end

