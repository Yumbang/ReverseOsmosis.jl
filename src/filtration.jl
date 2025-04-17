module Filtration

export
    # Filtration
    osmo_p, pump,
    element_filtration, element_filtration2, module_filtration, vessel_filtration

using Unitful
using ..ReverseOsmosis: Water, Water2, pressurize, profile_water, mix
using ..ReverseOsmosis: MembraneElement, MembraneElement2, MembraneModule, PressureVessel, profile_membrane, pristine_membrane

"""
Pressurize an `Water` with pressure to apply [Pa], with keyword argument `efficiency`.
Returns pressurized `Water` and consumed power [W] as following:
Power = (Feed volume⋅Aplied pressure)/(efficiency)
"""
function pump(feed::Union{Water, Water2}, pressure_to_apply::Float64; efficiency::Float64=0.8)
    @assert (0.0 ≤ efficiency ≤ 1.0) "Pump efficiency must be between 0 and 1."
    pressurized_water = pressurize(feed, pressure_to_apply)
    power_consumption = (feed.Q / 3600) * pressure_to_apply / efficiency * u"W"
    return (pressurized_water, power_consumption)
end

# =================================================== Definition of base functions for element filtration algorithms ===================================================

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

function osmo_p(water::Union{Water, Water2})::Float64
    osmo_p(water.C, water.T, water.m_avg)
end


"""
Calculate cake mass deposition flux (kg/m² ̇s) on membrane surface.

v: Cross-membrane flux [m/s]
u: Brine velocity [m/s]
cf: Foulant concentration [kg/m³]

mf_c: Cake mass deposition flux [kg/m²⋅s]

15.0        : Empirical criteria flux [L/m²⋅h]
0.1         : Empirical criteria velocity [m/s]
(u/0.1)^0.4 : Empirical conversion factor [-]

https://doi.org/10.1016/j.memsci.2008.01.030
"""
function cake_mass(v, u, cf)
    v_crit = 15.0 / 3600.0 / 1e3 * (u / 0.1)^0.4
    mf_c = cf * max(0.0, v - v_crit)
    return mf_c
end


"""
Calculate Cake-Enhanced Mass Transfer parameter (k_CEMT) based Sherwood number

u: Brine velocity [m/s]
H: Membrane height [m]
W: Membrane width [m]
μ: Water viscosity [Pa.s]
ρ: Water density [kg/m³]
ϵ: Cake layer porosity [-]
δ: Cake layer thickness [m]

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

# =================================================== Definition of element filtration algorithms ===================================================
# All of the results are unit-safe, I promise 😤.
# Using Unitful is an attractive idea, but it resulted in about x100 execution time.
# Element filtration functions are called tremendously often, thus they must be extremely fast.

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
    while err > 1e-6
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
    
    @assert (u_guess   > 0.0) "Brine flux is negative (u_guess: $(u_guess)). Backward flow isn't supported."
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

# =================================================== Definition of module / vessel filtration functions ===================================================
# Both element_filtration and element_filtration2 algorithms are used inside dispatches of the functions.
# If feed is Water and element is MembraneElement, element_filtration is used. If feed is Water2 and element is MembraneElement2, element_filtration2 is used.

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

