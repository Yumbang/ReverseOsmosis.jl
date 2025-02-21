module SemiBatch

export
    CirculationPipe, SemiBatchRO,
    process_semibatch_RO!

using DataFrames, Unitful, Printf

using DifferentialEquations: ODEProblem, solve, DP5
using DiffEqCallbacks:       SavedValues, SavingCallback

using ...ReverseOsmosis: Water, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, foul!

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
    const pipe::CirculationPipe
    pressure_vessel::PressureVessel
end

"""
Logger element for SBRO process.
An element is passed as a parameter of ODE, and saved with callback function.
"""
struct SBROLoggerElement
    raw_feed::Water
    final_feed::Water
    raw_brine::Water
    permeate::Water
    ΔR_ms_array::Union{Nothing, Array{Array{Float64}}}
    power_feed::typeof(1.0u"W")
    power_circ::typeof(1.0u"W")
    time::Float64
end

mutable struct SBROParameters
    feed_concentration::Float64
    feed_temperature::Float64
    feed_pressure::Float64
    flowrate_setpoint::Float64
    pressure_setpoint::Float64
    τ₁::Float64
    τ₂::Float64
    τ₃::Float64
    process::SemiBatchRO
    mode::Symbol
    logger::SBROLoggerElement
end

"""
Model junction that circulated brine and feed merge.
In most of the cases, both water are pressurized to same pressure setpoint.
If junction_loss is 0.0, there is no loss, assuming ideal junction.
"""
function junction(circulated_brine::Water, feed::Water; junction_loss::Float64)
    @assert (0.0 ≤ junction_loss ≤ 1.0)     "Junction loss must be between 0 and 1."
    merged_pressure = min(circulated_brine.P, feed.P) * (1.0-junction_loss)
    mixed_water     = mix(circulated_brine, feed; pressure=merged_pressure)
    return mixed_water
end


"""
ODE statement function of semi-batch RO.
Final feed means feed entering membrane vessel, after merged with circulated brine.
Unpressurized feed is considered to be unchanging.

u[1], u[2] : Flowrate, pressure of circulated brine (equal to those of final brine, ∵ water is assumed to be non-compressible)
u[2 + n]   : Concentration in n'th segment of pipe  (u[end] will be the circulated brine)
"""
function SBRO!(du, u, p, t)
    # Raw feed info
    feed_concentration  = p.feed_concentration
    feed_temperature    = p.feed_temperature
    feed_pressure       = p.feed_pressure
    # Setpoint parameters
    flowrate_setpoint   = p.flowrate_setpoint
    pressure_setpoint   = p.pressure_setpoint
    # Target process
    mode                = p.mode
    sbro                = p.process
    pipe_seg_volume     = sbro.pipe.volume / sbro.pipe.n_segments

    if mode == :CC
        # Unpressurized brine and feed
        unpressurized_brine = Water(u[1], feed_temperature, u[end], u[2])
        unpressurized_feed  = Water(max(flowrate_setpoint-u[1], 0), feed_temperature, feed_concentration, feed_pressure)
        # Pressurize feed and brine
        pressurized_feed,  power_feed = pump(unpressurized_feed, pressure_setpoint - unpressurized_feed.P; efficiency=0.8)
        pressurized_brine, power_circ = pump(unpressurized_brine, pressure_setpoint - unpressurized_brine.P; efficiency=0.8)
        # Final feed after mixing
        final_feed = junction(pressurized_brine, pressurized_feed; junction_loss=0.0)
    else
        # Unpressurized brine and feed
        unpressurized_brine = Water(u[2], feed_temperature, u[1], u[3])
        unpressurized_feed  = Water(flowrate_setpoint, feed_temperature, feed_concentration, feed_pressure)
        # Pressurize feed and (un)pressurized brine
        pressurized_feed,  power_feed = pump(unpressurized_feed, pressure_setpoint - unpressurized_feed.P; efficiency=0.8)
        pressurized_brine, power_circ = pump(unpressurized_brine, 0.0; efficiency=0.8) # No pressurization
        # Final feed without mixing
        final_feed = pressurized_feed
    end

    # RO modeling
    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    final_brine    = brines_array[end][end]
    final_permeate = mix(reduce(vcat, permeates_array); pressure=1e5)

    p.logger = SBROLoggerElement(
        unpressurized_feed, final_feed, final_brine, final_permeate,
        ΔR_ms_array, power_feed, power_circ, t
    )
    
    # ODE terms
    du[1] = 0.0
    du[2] = 0.0
    u[1]  = final_brine.Q
    u[2]  = final_brine.P

    du[3] = (final_brine.C - u[3]) / (pipe_seg_volume / u[1] * 3600.0)
    
    for n ∈ 2:sbro.pipe.n_segments
        du[n+2] = (u[n+1] - u[n+2]) / (pipe_seg_volume / u[1] * 3600.0)
    end

    return nothing
end

function sbro_save_func(u, t, integrator)
    return deepcopy(integrator.p.logger)
end

function process_semibatch_RO!(
    feed_T, feed_C, u₀,
    flowrate_setpoint, pressure_setpoint;
    process::SemiBatchRO, dt::Unitful.Time, mode::Symbol=:CC, fouling::Bool=true, τ::Tuple{Float64, Float64, Float64}=(1.0, 1.0, 1.0)
)
    @assert (mode ∈ [:CC, :purge]) "Only CC mode and purging mode are supported."

    # Unit conversion: Feed water quality
    feed_temperature_°C     = ustrip(uconvert(u"°C", feed_T))
    feed_concentration_kgm3 = ustrip(uconvert(u"kg/m^3", feed_C))
    
    # Unit conversion: Feed water quality
    flowrate_setpoint_m3h   = ustrip(uconvert(u"m^3/hr", flowrate_setpoint))
    pressure_setpoint_Pa    = ustrip(uconvert(u"Pa", pressure_setpoint))
    
    # Unit conversion: Feed water quality
    dt_sec                  = ustrip(uconvert(u"s", dt))

    # Setup logging components and callback function
    dummy_water   = Water(0.0, feed_temperature_°C, feed_concentration_kgm3, 1e5)
    dummy_logger  = SBROLoggerElement(dummy_water, dummy_water, dummy_water, dummy_water, nothing, 0.0u"W", 0.0u"W", 0.0)
    logging_array = SavedValues(typeof(0.0), typeof(dummy_logger))
    log_callback  = SavingCallback(sbro_save_func, logging_array)

    # Setup parameters and initial condition
    if isnothing(u₀)
        u₀ = [0.0, 1e5, fill(feed_concentration_kgm3, process.pipe.n_segments)...] # Configure u₀ with feed concentration, zero flow and atmospheric pressure.
    end

    params₀ = SBROParameters(
        feed_concentration_kgm3, feed_temperature_°C, 1e5,  # Feed info
        flowrate_setpoint_m3h, pressure_setpoint_Pa,        # Process setpoints info
        τ...,                                               # Lagging constants info
        process, mode,                                      # Target process info
        dummy_logger                                        # Logger info
    )
    
    # Define and solve SBRO ODE.
    # DP5 solver turned out to be more efficient compared to Tsit5.
    sbro_ODE = ODEProblem(SBRO!, u₀, (0.0, dt_sec), params₀)
    sbro_sol = solve(sbro_ODE, DP5(); callback=log_callback)

    # Store the last solution as initial condition of the next timestep
    next_u = deepcopy(sbro_sol.u[end])

    
    # Process the result into a dataframe
    time_profile             = DataFrame(:time=>logging_array.t[2:end]*u"s")
    raw_feed_profile         = profile_water(getfield.(logging_array.saveval, :raw_feed)[2:end])
    final_feed_profile       = profile_water(getfield.(logging_array.saveval, :final_feed)[2:end])
    raw_brine_profile        = profile_water(getfield.(logging_array.saveval, :raw_brine)[2:end])
    circulated_brine_profile = profile_water([Water(u[1], feed_temperature_°C, u[end], u[2]) for u in sbro_sol.u[2:end]])
    permeate_profile         = profile_water(getfield.(logging_array.saveval, :permeate)[2:end])
    power_profile            = DataFrame(
                                    :HPP_power=> getfield.(logging_array.saveval, :power_feed)[2:end],
                                    :Circ_power=>getfield.(logging_array.saveval, :power_circ)[2:end]
                                )

    # Rename dataframes before merging
    rename!((x -> "Raw_feed_"*x),         raw_feed_profile)
    rename!((x -> "Final_feed_"*x),       final_feed_profile)
    rename!((x -> "Raw_brine_"*x),        raw_brine_profile)
    rename!((x -> "Circulated_brine_"*x), circulated_brine_profile) # If mode is :purge, it is disposed brine's profile.
    rename!((x -> "Permeate_"*x),         permeate_profile)
    
    # Merge dataframes into single result dataframe and add time difference (kinda little `dt`) column
    result_profile = hcat(
        time_profile, raw_feed_profile, final_feed_profile, 
        raw_brine_profile, circulated_brine_profile, permeate_profile, power_profile
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

    return result_profile, next_u
end

end # SemiBatchODE