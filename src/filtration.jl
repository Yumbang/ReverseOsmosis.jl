module Filtration

export
    # Filtration
    osmo_p, pump,
    element_filtration, module_filtration, vessel_filtration

using Unitful
using ..ReverseOsmosis: Water, pressurize, profile_water, mix
using ..ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel, profile_membrane, pristine_membrane


"""
Calculates osmotic pressure of water based on van't Hoff equation.

concentration   : Mass concentration (TDS) as kg/m³
m_avg           : Average ion molecular weight (if not given, assume most of the ions are NaCl)
"""
function osmo_p(concentration::Float64, temperature::Float64)::Float64
    return 2 / 58.44 * concentration * 8.3145e3 * (temperature + 273.15) # Assuming all of the ions are Na⁺ ions
end

function osmo_p(concentration::Float64, temperature::Float64, m_avg::Float64)::Float64
    return 2 / m_avg * concentration * 8.3145e3 * (temperature + 273.15) # Assuming all of the ions are Na⁺ ions
end

function osmo_p(water::Water)::Float64
    osmo_p(water.C, water.T, water.m_avg)
end

"""
Pressurize an `Water` with pressure to apply [Pa], with keyword argument `efficiency`.
Returns pressurized `Water` and consumed power [W] as following:
Power = (Feed volume⋅Aplied pressure)/(efficiency)
"""
function pump(feed::Water, pressure_to_apply::Float64; efficiency::Float64=0.8)
    @assert (0.0 ≤ efficiency ≤ 1.0) "Pump efficiency must be between 0 and 1."
    pressurized_water = pressurize(feed, pressure_to_apply)
    power_consumption = (feed.Q / 3600) * pressure_to_apply / efficiency * u"W"
    return (pressurized_water, power_consumption)
end

"""
Model reverse osmosis of feed water on single membrane element.
Returns brine & permeate and resistance increase (i.e. fouling effect).
"""
function element_filtration(
    feed::Water, element::MembraneElement; dt::Union{Float64, Nothing}=nothing
)::Tuple{Water, Water, Float64}
    C_feed = feed.C
    Q_feed = feed.Q
    P_feed = feed.P
    T_feed = feed.T
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
        osmo_p_guess = osmo_p(c_guess, T_feed, feed.m_avg)
        v_w_guess    = max(P_feed - osmo_p_guess, 0.0) / element.R_m
        u_guess      = U_feed - v_w_guess * element.dx / element.height
        c_cal        = ((C_feed * U_feed * element.height) - ((1-reject)*C_feed) * (v_w_guess * element.dx)) / (u_guess * element.height)
        err          = abs((c_cal - c_guess) / c_cal)
        idx_iter     += 1
        if idx_iter ≥ 100
            @error "Failed to converge. Iteration exceeded 100 times."
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

function module_filtration(
    feed::Water, membrane::MembraneModule; dt::Union{Float64, Nothing}=nothing
)::Tuple{Array{Water}, Array{Water}, Array{Float64}}
    local next_feed::Water
    local brine::Water
    local permeate::Water
    local ΔR_m::Float64

    brines    = []
    permeates = []
    ΔR_ms     = []
    
    next_feed = feed
    for element in membrane.elements_array
        brine, permeate, ΔR_m = element_filtration(next_feed, element; dt=dt)
        push!(brines, brine)
        push!(permeates, permeate)
        push!(ΔR_ms, ΔR_m)
        next_feed = brine
    end

    return brines, permeates, ΔR_ms
end

function vessel_filtration(
    feed::Water, vessel::PressureVessel; dt::Union{Float64, Nothing}=nothing
)
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

end

