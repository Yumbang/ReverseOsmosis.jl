module SemiBatch

export
    CirculationPipe, SemiBatchRO,
    process_semibatch_RO!

using DataFrames, Unitful, Printf

using DifferentialEquations: ODEProblem, solve, DP5
using DiffEqCallbacks:       SavedValues, SavingCallback

using ...ReverseOsmosis: Water, Water2, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneElement2, MembraneModule, PressureVessel, foul!

using ...ReverseOsmosis: pump, vessel_filtration


"""
Pipe full of fluid, with `length` and cross-section `area`.
A pipe is considered to be completely-stirred as concentration `C`.

This struct—designed as brine circulation line—is necessary to properly model brine circulation.
"""
struct CirculationPipe
    volume::Float64       # [m³]
    n_segments::Int64
end

"""
Simple semi-batch reverse osmosis process structure.
Feed is pressurized by a high-pressure pump and brine is circulated with circulation pump.

SemiBatchRO is 'stateful', it holds the pressure vessel and brine information.
If mode is `:CC`, the brine is mixed into feed, and purged out otherwise.

At the first evaluation, circulated_brine is `nothing`. Process is automatically evaluated without brine circulation (`:purge` mode).

Because of the property of semi-batch RO, using long `dt` is discouraged. Increasing the temporal resolution is recommended.
"""
mutable struct SemiBatchRO
    const high_pressure_pump_efficiency::Float64
    const circulation_pump_efficiency::Float64
    const pipe_to_junction::CirculationPipe
    const pipe_to_membrane::CirculationPipe
    pressure_vessel::PressureVessel
end

"""
Logger element for SBRO process (using element_filtration algorithm).
An element is passed as a parameter of ODE, and saved with callback function.
"""
struct SBROLoggerElement
    raw_feed::Water
    final_feed::Water
    raw_brine::Water
    junction_brine::Water
    permeate::Water
    ΔR_ms_array::Union{Nothing, Array{Array{Float64}}}
    power_feed::typeof(1.0u"W")
    power_circ::typeof(1.0u"W")
    u::Union{Nothing, Vector{Float64}}
    du::Union{Nothing, Vector{Float64}}
    time::Float64
end

mutable struct SBROParameters
    feed_concentration::Float64
    feed_temperature::Float64
    feed_pressure::Float64
    feed_molar_mass::Float64
    # brine_molar_mass::Float64
    flowrate_setpoint::Float64
    pressure_setpoint::Float64
    process::SemiBatchRO
    mode::Symbol
    logger::SBROLoggerElement
end

"""
Logger element for SBRO process (using element_filtration2 algorithm).
An element is passed as a parameter of ODE, and saved with callback function.
"""
struct SBROLoggerElement2
    raw_feed::Water2
    final_feed::Water2
    raw_brine::Water2
    junction_brine::Water2
    permeate::Water2
    ΔCake_thicknesses_array::Union{Nothing, Array{Array{Float64}}}
    ΔR_cs_array::Union{Nothing, Array{Array{Float64}}}
    power_feed::typeof(1.0u"W")
    power_circ::typeof(1.0u"W")
    u::Union{Nothing, Vector{Float64}}
    du::Union{Nothing, Vector{Float64}}
    time::Float64
end

mutable struct SBROParameters2
    feed_concentration::Float64
    feed_temperature::Float64
    feed_pressure::Float64
    feed_foulant::Float64
    feed_molar_mass::Float64
    # brine_molar_mass::Float64
    flowrate_setpoint::Float64
    pressure_setpoint::Float64
    process::SemiBatchRO
    mode::Symbol
    logger::SBROLoggerElement2
end


function sbro_save_func(u, t, integrator)
    return deepcopy(integrator.p.logger)
end


"""
Model junction that circulated brine and feed merge.
In most of the cases, both water are pressurized to same pressure setpoint.
If junction_loss is 0.0, there is no loss, assuming ideal junction.
"""
function junction(circulated_brine::Union{Water, Water2}, feed::Union{Water, Water2}; junction_loss::Float64)
    @assert (0.0 ≤ junction_loss ≤ 1.0)     "Junction loss must be between 0 and 1."
    merged_pressure = min(circulated_brine.P, feed.P) * (1.0-junction_loss)
    mixed_water     = mix(circulated_brine, feed; pressure=merged_pressure)
    return mixed_water
end

"""
ODE statement function of semi-batch RO using ElementFiltration algorithm.
Final feed means feed entering membrane vessel, after merged with circulated brine.
Unpressurized feed is considered to be unchanging.
Fresh brine means brine produced right after RO membrane.

u[1]: Final feed flowrate
u[2]: Final feed pressure
u[3]: Fresh brine flowrate
u[4]: Fresh brine pressure
u[5]: Fresh brine concentration
u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration
u[6+N] to u[5+N+M] : Junction -> Membrane pipe segments concentration
u[5+N+M+1] : Final feed average molar mass
"""
function SBRO!(du, u, p, t)
    # Raw feed info
    feed_concentration  = p.feed_concentration
    feed_temperature    = p.feed_temperature
    feed_pressure       = p.feed_pressure
    feed_molar_mass     = p.feed_molar_mass
    # brine_molar_mass    = p.brine_molar_mass
    # Setpoint parameters
    flowrate_setpoint   = p.flowrate_setpoint
    pressure_setpoint   = p.pressure_setpoint
    # Target process
    mode                = p.mode
    sbro                = p.process
    N = sbro.pipe_to_junction.n_segments
    M = sbro.pipe_to_membrane.n_segments
    pipe_to_junction_seg_volume = sbro.pipe_to_junction.volume / N
    pipe_to_membrane_seg_volume = sbro.pipe_to_membrane.volume / M

    if feed_pressure > pressure_setpoint
        @warn "Feed pressure ($(feed_pressure)) is already higher than setpoint ($(pressure_setpoint))"
    end

    final_feed = Water(u[1], feed_temperature, u[5+N+M], u[2], u[5+N+M+1])
    

    # RO modeling
    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    # fake_final_brine = brines_array[end][end]
    # final_brine      = Water(
    #     fake_final_brine.Q, fake_final_brine.T, fake_final_brine.C, fake_final_brine.P, brine_molar_mass
    # )
    final_brine      = brines_array[end][end]
    final_permeate   = mix(reduce(vcat, permeates_array); pressure=1e5)

    # ODE terms
    
    # Membrane to Junction modeling
    # u[3]: Fresh brine flowrate
    # u[4]: Fresh brine pressure
    # u[5]: Fresh brine concentration
    # u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration
    
    du[3] = 0.0
    du[4] = 0.0
    u[3]  = final_brine.Q
    u[4]  = final_brine.P
    if mode == :CC
        du[5] = 0.0
        u[5]  = final_brine.C
        # du[5] = (final_brine.C - u[5]) / (pipe_to_junction_seg_volume / u[3] * 3600.0)
        for n ∈ 1:N
            du[n+5] = (u[n+4] - u[n+5]) / (pipe_to_junction_seg_volume / u[3] * 3600.0)
        end
    else
        du[6:5+N] .= 0.0
    end

    
    # Junction to Membrane modeling
    # u[1]: Final feed flowrate
    # u[2]: Final feed pressure
    # u[3]: Fresh brine flowrate
    # u[4]: Fresh brine pressure
    # u[5]: Fresh brine concentration
    # u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration
    # u[6+N] to u[5+N+M] : Junction -> Membrane pipe segments concentration

    if mode == :CC
        fresh_feed_Q      = flowrate_setpoint - u[3]
        junction_brine    = Water(u[3], feed_temperature, u[5+N], u[4], feed_molar_mass)
        power_feed        = fresh_feed_Q / 3600 * max(0.0, pressure_setpoint-feed_pressure) / sbro.high_pressure_pump_efficiency * u"W"
        power_circ        = u[3]         / 3600 * (pressure_setpoint - u[4])                / sbro.circulation_pump_efficiency   * u"W"
    else
        fresh_feed_Q      = flowrate_setpoint
        junction_brine    = Water(0.0,  feed_temperature, u[5+N], u[4], feed_molar_mass)
        # junction_brine    = Water(u[3],  feed_temperature, u[5+N], u[4], feed_molar_mass)
        power_feed        = fresh_feed_Q / 3600 * max(0.0, pressure_setpoint-feed_pressure) / sbro.high_pressure_pump_efficiency * u"W"
        power_circ        = 0.0u"W"
    end

    fresh_feed          = Water(fresh_feed_Q, feed_temperature, feed_concentration, feed_pressure, feed_molar_mass)
    junction_mixed_feed = mix(fresh_feed, junction_brine; pressure=pressure_setpoint)
    
    du[1] = 0.0
    du[2] = 0.0
    du[5+N+M+1] = 0.0
    u[1]  = flowrate_setpoint
    u[2]  = pressure_setpoint
    u[5+N+M+1] = junction_mixed_feed.m_avg

    du[6+N] = (junction_mixed_feed.C - u[6+N]) / (pipe_to_membrane_seg_volume / u[1] * 3600.0)
    for m ∈ 1:M-1
        du[6+N+m] = (u[5+N+m] - u[6+N+m]) / (pipe_to_membrane_seg_volume / u[1] * 3600.0)
    end
    
    p.logger = SBROLoggerElement(
        fresh_feed, final_feed, final_brine, junction_brine, final_permeate,
        ΔR_ms_array, power_feed, power_circ, deepcopy(u), deepcopy(du), t
    )

    return nothing
end

"""
ODE statement function of semi-batch RO using ElementFiltration2 algorithm.
Final feed means feed entering membrane vessel, after merged with circulated brine.
Unpressurized feed is considered to be unchanging.
Fresh brine means brine produced right after RO membrane.

u[1]: Final feed flowrate
u[2]: Final feed pressure
u[3]: Fresh brine flowrate
u[4]: Fresh brine pressure
u[5]: Fresh brine concentration
u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration (salt)
u[6+N] to u[5+N+M] : Junction -> Membrane pipe segments concentration (salt)
u[5+N+M+1] : Final feed average molar mass
u[5+N+M+1+1] to u[5+N+M+1+n_foulant_segments]: Membrane -> Junction pipe segments concentration (foulant)
u[5+N+M+1+n_foulant_segments+1] to u[5+N+M+1+2*n_foulant_segments]: Junction -> Membrane pipe segments concentration (foulant)
"""
function SBRO2!(du, u, p, t)
    # Raw feed info
    feed_concentration  = p.feed_concentration
    feed_foulant        = p.feed_foulant
    feed_temperature    = p.feed_temperature
    feed_pressure       = p.feed_pressure
    feed_molar_mass     = p.feed_molar_mass
    # brine_molar_mass    = p.brine_molar_mass
    # Setpoint parameters
    flowrate_setpoint   = p.flowrate_setpoint
    pressure_setpoint   = p.pressure_setpoint
    # Target process
    mode                = p.mode
    sbro                = p.process
    N = sbro.pipe_to_junction.n_segments
    M = sbro.pipe_to_membrane.n_segments
    pipe_to_junction_seg_volume = sbro.pipe_to_junction.volume / N
    pipe_to_membrane_seg_volume = sbro.pipe_to_membrane.volume / M
    n_foulant_segments  = 3


    # final_feed = Water2(u[1], feed_temperature, u[5+N+M], u[5+N+M+1+2*n_foulant_segments], u[2], u[5+N+M+1])
    final_feed = Water2(u[1], feed_temperature, u[5+N+M], u[5+N+M+1+2*n_foulant_segments], u[2], feed_molar_mass)

    # RO modeling
    brines_array, permeates_array, ΔCake_thicknesses_array, ΔR_cs_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    final_brine      = brines_array[end][end]
    final_permeate   = mix(reduce(vcat, permeates_array); pressure=1e5)

    # ODE terms
    
    # Membrane to Junction modeling
    # u[3]: Fresh brine flowrate
    # u[4]: Fresh brine pressure
    # u[5]: Fresh brine concentration
    # u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration
    # u[5+N+M+1+1] to u[5+N+M+1+n_foulant_segments]: Membrane -> Junction pipe segments concentration (foulant)
    
    du[3] = 0.0
    du[4] = 0.0
    u[3]  = final_brine.Q
    u[4]  = final_brine.P
    if mode == :CC
        du[5] = 0.0
        u[5]  = final_brine.C
        # du[5] = (final_brine.C - u[5]) / (pipe_to_junction_seg_volume / u[3] * 3600.0)
        for n ∈ 1:N
            du[n+5] = (u[n+4] - u[n+5]) / (pipe_to_junction_seg_volume / u[3] * 3600.0)
        end
        
        du[5+N+M+1+1] = (final_brine.Cf - u[5+N+M+1 + 1]) / (sbro.pipe_to_junction.volume / n_foulant_segments / u[3] * 3600.0)
        for foulant_n ∈ 2:n_foulant_segments
            du[5+N+M+1+foulant_n] = (u[5+N+M+foulant_n] - u[5+N+M+1+foulant_n]) / (sbro.pipe_to_junction.volume / n_foulant_segments / u[3] * 3600.0)
        end
    else
        du[6:5+N] .= 0.0
        du[5+N+M+1 + 1:5+N+M+1 + n_foulant_segments] .= 0.0
    end

    
    # Junction to Membrane modeling
    # u[1]: Final feed flowrate
    # u[2]: Final feed pressure
    # u[3]: Fresh brine flowrate
    # u[4]: Fresh brine pressure
    # u[5]: Fresh brine concentration
    # u[6]   to u[5+N]   : Membrane -> Junction pipe segments concentration
    # u[6+N] to u[5+N+M] : Junction -> Membrane pipe segments concentration
    # u[5+N+M+1+1] to u[5+N+M+1+n_foulant_segments]: Membrane -> Junction pipe segments concentration (foulant)
    # u[5+N+M+1+n_foulant_segments+1] to u[5+N+M+1+2*n_foulant_segments]: Junction -> Membrane pipe segments concentration (foulant)

    if mode == :CC
        fresh_feed_Q      = flowrate_setpoint - u[3]
        junction_brine    = Water2(u[3], feed_temperature, u[5+N], u[5+N+M+1+n_foulant_segments], u[4], feed_molar_mass)
        power_feed        = fresh_feed_Q / 3600 * (pressure_setpoint-feed_pressure) / sbro.high_pressure_pump_efficiency * u"W"
        power_circ        = u[3]         / 3600 * (pressure_setpoint - u[4])        / sbro.circulation_pump_efficiency   * u"W"
    else
        fresh_feed_Q      = flowrate_setpoint
        junction_brine    = Water2(0.0,  feed_temperature, u[5+N], u[5+N+M+1+n_foulant_segments], u[4], feed_molar_mass)
        power_feed        = fresh_feed_Q / 3600 * (pressure_setpoint-feed_pressure) / sbro.high_pressure_pump_efficiency * u"W"
        power_circ        = 0.0u"W"
    end

    fresh_feed          = Water2(fresh_feed_Q, feed_temperature, feed_concentration, feed_foulant, feed_pressure, feed_molar_mass)
    junction_mixed_feed = mix(fresh_feed, junction_brine; pressure=pressure_setpoint)
    
    du[1] = 0.0
    du[2] = 0.0
    du[5+N+M+1] = 0.0
    u[1]  = flowrate_setpoint
    u[2]  = pressure_setpoint
    u[5+N+M+1] = junction_mixed_feed.m_avg

    du[6+N] = (junction_mixed_feed.C - u[6+N]) / (pipe_to_membrane_seg_volume / u[1] * 3600.0)
    for m ∈ 1:M-1
        du[6+N+m] = (u[5+N+m] - u[6+N+m]) / (pipe_to_membrane_seg_volume / u[1] * 3600.0)
    end

    du[5+N+M+1+n_foulant_segments+1] = (junction_mixed_feed.Cf - u[5+N+M+1+n_foulant_segments+1]) / (sbro.pipe_to_membrane.volume / n_foulant_segments / u[1] * 3600.0)
    for foulant_m ∈ 2:n_foulant_segments
        du[5+N+M+1+n_foulant_segments+foulant_m] = (u[5+N+M+n_foulant_segments+foulant_m] - u[5+N+M+1+n_foulant_segments+foulant_m]) / (sbro.pipe_to_junction.volume / n_foulant_segments / u[3] * 3600.0)
    end
    
    p.logger = SBROLoggerElement2(
        fresh_feed, final_feed, final_brine, junction_brine, final_permeate,
        ΔCake_thicknesses_array, ΔR_cs_array, power_feed, power_circ, deepcopy(u), deepcopy(du), t
    )

    return nothing
end

function process_semibatch_RO!(
    feed::Water, brine_M::typeof(1.0u"g/mol"), u₀::Union{Nothing, Vector{Float64}},
    flowrate_setpoint::typeof(1.0u"m^3/hr"), pressure_setpoint::typeof(1.0u"bar");
    process::SemiBatchRO, dt::Unitful.Time, mode::Symbol=:CC, fouling::Bool=true
)
    @assert (mode ∈ [:CC, :purge]) "Only CC mode and purging mode are supported."

    # Unit conversion: Feed water quality
    feed_temperature_°C     = feed.T
    feed_concentration_kgm3 = feed.C
    feed_pressure_Pa        = feed.P
    feed_molar_mass_gmol    = feed.m_avg
    # brine_molar_mass_gmol   = ustrip(uconvert(u"g/mol", brine_M))
    
    # Unit conversion: Feed water quality
    flowrate_setpoint_m3h   = ustrip(uconvert(u"m^3/hr", flowrate_setpoint))
    pressure_setpoint_Pa    = ustrip(uconvert(u"Pa", pressure_setpoint))
    
    # Unit conversion: Feed water quality
    dt_sec                  = ustrip(uconvert(u"s", dt))

    # Setup logging components and callback function
    dummy_water   = Water(0.0, feed_temperature_°C, feed_concentration_kgm3, 1e5)
    dummy_logger  = SBROLoggerElement(dummy_water, dummy_water, dummy_water, dummy_water, dummy_water, nothing, 0.0u"W", 0.0u"W", nothing, nothing, 0.0)
    logging_array = SavedValues(typeof(0.0), typeof(dummy_logger))
    log_callback  = SavingCallback(sbro_save_func, logging_array)

    # Setup parameters and initial condition
    if isnothing(u₀)
        u₀ = [flowrate_setpoint_m3h, 1e5, 0.0, 1e5, fill(feed_concentration_kgm3, process.pipe_to_junction.n_segments+process.pipe_to_membrane.n_segments+1)..., feed_molar_mass_gmol] # Configure u₀ with feed concentration, zero flow and atmospheric pressure.
    end

    params₀ = SBROParameters(
        feed_concentration_kgm3, feed_temperature_°C, feed_pressure_Pa, feed_molar_mass_gmol, # brine_molar_mass_gmol,  # Feed & brine info
        flowrate_setpoint_m3h, pressure_setpoint_Pa,        # Process setpoints info
        process, mode,                                      # Target process info
        dummy_logger                                        # Logger info
    )
    
    # Define and solve SBRO ODE.
    # DP5 solver turned out to be more efficient compared to Tsit5.
    sbro_ODE = ODEProblem(SBRO!, u₀, (0.0, dt_sec), params₀)
    sbro_sol = solve(sbro_ODE, DP5(); callback=log_callback)
    # sbro_sol = solve(sbro_ODE, Tsit5(); callback=log_callback, dtmax=1.0)

    # Store the last solution as initial condition of the next timestep
    next_u = deepcopy(sbro_sol.u[end])

    
    # Process the result into a dataframe
    time_profile             = DataFrame(:time=>logging_array.t[2:end]*u"s")
    raw_feed_profile         = profile_water(getfield.(logging_array.saveval, :raw_feed)[2:end])
    final_feed_profile       = profile_water(getfield.(logging_array.saveval, :final_feed)[2:end])
    raw_brine_profile        = profile_water(getfield.(logging_array.saveval, :raw_brine)[2:end])
    junction_brine_profile   = profile_water(getfield.(logging_array.saveval, :junction_brine)[2:end])
    permeate_profile         = profile_water(getfield.(logging_array.saveval, :permeate)[2:end])
    power_profile            = DataFrame(
                                    :HPP_power=> getfield.(logging_array.saveval, :power_feed)[2:end],
                                    :Circ_power=>getfield.(logging_array.saveval, :power_circ)[2:end]
                                )

    u_log  = getfield.(logging_array.saveval, :u)[2:end]
    du_log = getfield.(logging_array.saveval, :du)[2:end]

    # Rename dataframes before merging
    rename!((x -> "Raw_feed_"*x),         raw_feed_profile)
    rename!((x -> "Final_feed_"*x),       final_feed_profile)
    rename!((x -> "Raw_brine_"*x),        raw_brine_profile)
    rename!((x -> "Junction_brine_"*x),   junction_brine_profile)
    rename!((x -> "Permeate_"*x),         permeate_profile)
    
    # Merge dataframes into single result dataframe and add time difference (kinda little `dt`) column
    result_profile = hcat(
        time_profile, raw_feed_profile, final_feed_profile, 
        raw_brine_profile, junction_brine_profile, permeate_profile, power_profile
    )
    result_profile.time_diff = diff(logging_array.t)*u"s"
    result_profile.mode     .= mode

    # Process fouling if fouling is true
    if fouling
        ΔR_ms = getfield.(logging_array.saveval, :ΔR_ms_array)[2:end] .* ustrip.(result_profile.time_diff)
        for ΔR_m ∈ ΔR_ms
            foul!(process.pressure_vessel, ΔR_m)
        end
    end

    return result_profile, next_u, u_log, du_log
end

function process_semibatch_RO!(
    feed::Water2, brine_M::typeof(1.0u"g/mol"), u₀::Union{Nothing, Vector{Float64}},
    flowrate_setpoint::typeof(1.0u"m^3/hr"), pressure_setpoint::typeof(1.0u"bar");
    process::SemiBatchRO, dt::Unitful.Time, mode::Symbol=:CC, fouling::Bool=true
)
    @assert (mode ∈ [:CC, :purge]) "Only CC mode and purging mode are supported."

    # Unit conversion: Feed water quality
    feed_temperature_°C     = feed.T
    feed_concentration_kgm3 = feed.C
    feed_foulant_kgm3       = feed.Cf
    feed_pressure_Pa        = feed.P
    feed_molar_mass_gmol    = feed.m_avg
    # brine_molar_mass_gmol   = ustrip(uconvert(u"g/mol", brine_M))
    
    # Unit conversion: Feed water quality
    flowrate_setpoint_m3h   = ustrip(uconvert(u"m^3/hr", flowrate_setpoint))
    pressure_setpoint_Pa    = ustrip(uconvert(u"Pa", pressure_setpoint))
    
    # Unit conversion: Feed water quality
    dt_sec                  = ustrip(uconvert(u"s", dt))

    # Setup logging components and callback function
    dummy_water   = Water2(0.0, feed_temperature_°C, feed_concentration_kgm3, feed_foulant_kgm3, 1e5, feed_molar_mass_gmol)
    dummy_logger  = SBROLoggerElement2(dummy_water, dummy_water, dummy_water, dummy_water, dummy_water, nothing, nothing, 0.0u"W", 0.0u"W", nothing, nothing, 0.0)
    logging_array = SavedValues(typeof(0.0), typeof(dummy_logger))
    log_callback  = SavingCallback(sbro_save_func, logging_array)

    # Setup parameters and initial condition
    if isnothing(u₀)
        u₀ = [flowrate_setpoint_m3h, 1e5, 0.0, 1e5,     # Feed informations
            fill(feed_concentration_kgm3, process.pipe_to_junction.n_segments + process.pipe_to_membrane.n_segments + 1)..., # Fill circulation pipe with feed, if u₀ is nothing
            feed_molar_mass_gmol,           # Specify molar weight of the feed (will be depreciated in the next version)
            fill(feed_foulant_kgm3, 6)...]  # Fill circulation pipe with feed, if u₀ is nothing. Currently, foulant circulation is modeled only with 6 pipe segments.
    end
    params₀ = SBROParameters2(
        feed_concentration_kgm3, feed_temperature_°C, feed_pressure_Pa, feed_foulant_kgm3, feed_molar_mass_gmol, # brine_molar_mass_gmol,  # Feed & brine info
        flowrate_setpoint_m3h, pressure_setpoint_Pa,        # Process setpoints info
        process, mode,                                      # Target process info
        dummy_logger                                        # Logger info
    )
    
    # Define and solve SBRO ODE.
    # DP5 solver turned out to be more efficient compared to Tsit5.
    sbro_ODE2 = ODEProblem(SBRO2!, u₀, (0.0, dt_sec), params₀)
    sbro_sol = solve(sbro_ODE2, DP5(); callback=log_callback)
    # sbro_sol = solve(sbro_ODE, Tsit5(); callback=log_callback, dtmax=1.0)

    # Store the last solution as initial condition of the next timestep
    next_u = deepcopy(sbro_sol.u[end])

    
    # Process the result into a dataframe
    time_profile             = DataFrame(:time=>logging_array.t[2:end]*u"s")
    raw_feed_profile         = profile_water(getfield.(logging_array.saveval, :raw_feed)[2:end])
    final_feed_profile       = profile_water(getfield.(logging_array.saveval, :final_feed)[2:end])
    raw_brine_profile        = profile_water(getfield.(logging_array.saveval, :raw_brine)[2:end])
    junction_brine_profile   = profile_water(getfield.(logging_array.saveval, :junction_brine)[2:end])
    permeate_profile         = profile_water(getfield.(logging_array.saveval, :permeate)[2:end])
    power_profile            = DataFrame(
                                    :HPP_power=> getfield.(logging_array.saveval, :power_feed)[2:end],
                                    :Circ_power=>getfield.(logging_array.saveval, :power_circ)[2:end]
                                )

    u_log  = getfield.(logging_array.saveval, :u)[2:end]
    du_log = getfield.(logging_array.saveval, :du)[2:end]

    # Rename dataframes before merging
    rename!((x -> "Raw_feed_"*x),         raw_feed_profile)
    rename!((x -> "Final_feed_"*x),       final_feed_profile)
    rename!((x -> "Raw_brine_"*x),        raw_brine_profile)
    rename!((x -> "Junction_brine_"*x),   junction_brine_profile)
    rename!((x -> "Permeate_"*x),         permeate_profile)
    
    # Merge dataframes into single result dataframe and add time difference (kinda little `dt`) column
    result_profile = hcat(
        time_profile, raw_feed_profile, final_feed_profile, 
        raw_brine_profile, junction_brine_profile, permeate_profile, power_profile
    )
    result_profile.time_diff = diff(logging_array.t)*u"s"
    result_profile.mode     .= mode

    # Process fouling if fouling is true
    if fouling
        ΔCake_thicknesses = getfield.(logging_array.saveval, :ΔCake_thicknesses_array)[2:end] .* ustrip.(result_profile.time_diff)
        ΔR_cs = getfield.(logging_array.saveval, :ΔR_cs_array)[2:end] .* ustrip.(result_profile.time_diff)

        # for (ΔCake_thickness, ΔR_c) ∈ (ΔCake_thicknesses, ΔR_cs)
        #     @info "Average cake thickness increase" ΔCake_thickness
        #     foul!(process.pressure_vessel, ΔCake_thickness, ΔR_c)
        # end
        for module_idx ∈ 1:length(process.pressure_vessel.modules_array)
            foul!(process.pressure_vessel, ΔCake_thicknesses[module_idx], ΔR_cs[module_idx])
        end
    end

    return result_profile, next_u, u_log, du_log
end

end # SemiBatchODE