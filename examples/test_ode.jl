using Pkg
Pkg.activate("/Users/ybang_mac/Documents/병찬/Research/2025/Batch RO/Code/ReverseOsmosis.jl")
using ReverseOsmosis

using DataFrames, Unitful, Printf

using DifferentialEquations
using DiffEqCallbacks
using Plots
using BenchmarkTools

membrane_spec = (25, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
pristine_vessel = PressureVessel(7, membrane);

"""
Pipe full of fluid, with `volume` m³.
A pipe is considered to be completely-stirred.

This struct—designed as brine circulation line—is necessary to properly model brine circulation.
"""
struct Pipe
    volume::Float64       # [m³]
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
    const pipe::Pipe
    pressure_vessel::PressureVessel
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

struct SBROLoggerElement
    brine::Water
    permeate::Water
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
    process::SemiBatchRO
    logger::SBROLoggerElement
end

"""
ODE statement function of semi-batch RO.
Final feed means feed entering membrane vessel, after merged with circulated brine.
Unpressurized feed is considered to be unchanging.

u[1]: Pipe (circulated brine) concentration
u[2]: Pipe (circulated brine) flow rate
u[3]: Pipe (circulated brine) pressure
"""
function SBRO!(du, u, p, t)
    feed_concentration  = p.feed_concentration
    feed_temperature    = p.feed_temperature
    feed_pressure       = p.feed_pressure
    flowrate_setpoint   = p.flowrate_setpoint
    pressure_setpoint   = p.pressure_setpoint
    sbro                = p.process

    # Unpressurized brine and feed
    unpressurized_brine = Water(u[2], feed_temperature, u[1], u[3])
    unpressurized_feed  = Water(max(flowrate_setpoint-u[2], 0), feed_temperature, feed_concentration, feed_pressure)

    # Pressurize feed and brine
    pressurized_feed,  power_feed = pump(unpressurized_feed, pressure_setpoint - unpressurized_feed.P; efficiency=0.8)
    pressurized_brine, power_circ = pump(unpressurized_brine, pressure_setpoint - unpressurized_brine.P; efficiency=0.8)

    # Final feed after mixing
    final_feed = junction(pressurized_brine, pressurized_feed; junction_loss=0.0)

    # RO modeling
    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    final_brine    = brines_array[end][end]
    final_permeate = mix(reduce(vcat, permeates_array); pressure=1e-5)

    p.logger = SBROLoggerElement(final_brine, final_permeate, power_feed, power_circ, t)
    
    # ODE terms
    # du[1] = (final_brine.C * final_brine.Q - u[1] * final_brine.Q) / sbro.pipe.volume
    du[1] = (final_brine.C * final_brine.Q - u[1] * final_brine.Q) / sbro.pipe.volume
    du[2] = final_brine.Q - u[2]
    du[3] = final_brine.P - u[3]

    return nothing
end

pipe = Pipe(5.0)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(pristine_vessel))
feed = Water(100/84, 20.0, 300/1000, 1e5);
logger = SBROLoggerElement(feed, feed, 0.0u"W", 0.0u"W", 0.0)
logger_array = SBROLoggerElement[]

function save_func(u, t, integrator)
    return deepcopy(integrator.p.logger)
end

saved_values = SavedValues(typeof(0.0), typeof(logger))
cb = SavingCallback(save_func, saved_values)

u₀ = [feed.C, 500/84, 10e5]
# u₀ = [feed.C, 0.0, 0.0]
params = SBROParameters(
    300/1000, 20.0, 1e-5,
    1000/84, 12e5, sbro, logger
)
prob = ODEProblem(SBRO!, u₀, (0.0, 60.0), params)

sol = solve(prob, Tsit5(); callback=cb)

# @benchmark solve($prob, $Tsit5(); callback=$cb)
# @benchmark solve($prob, $DP5(); callback=$cb)


time_profile = DataFrame(:time=>saved_values.t[2:end]*u"s")
brine_profile = profile_water(getfield.(saved_values.saveval, :brine)[2:end])
permeate_profile = profile_water(getfield.(saved_values.saveval, :permeate)[2:end])
power_profile = DataFrame(
    :power_feed=>getfield.(saved_values.saveval, :power_feed)[2:end],
    :power_circ=>getfield.(saved_values.saveval, :power_circ)[2:end]
)

sol.u[2:end]

plot(time_profile.time, brine_profile.Q, label="Brine Q", title="Flowrate")
plot!(time_profile.time, [u[2] for u in sol.u[2:end]], label="Circulated brine Q")
plot!(time_profile.time, permeate_profile.Q, label="Permeate Q")
plot!(time_profile.time, [1000/84 for _ in time_profile.time], label="Merged feed Q")
plot!(time_profile.time, permeate_profile.Q+brine_profile.Q, label="Brine + Permeate Q")

plot(time_profile.time, power_profile.power_feed, label="HPP power use")
plot!(time_profile.time, power_profile.power_circ, label="Circulation pump power use")
plot!(time_profile.time, power_profile.power_feed+power_profile.power_circ, label="Total power use")

plot(time_profile.time, (power_profile.power_feed+power_profile.power_circ)./(permeate_profile.Q), label="SEC")

plot(time_profile.time, [u[1]*u"kg/m^3" for u in sol.u[2:end]], label="Circulated brine C", title="Concentration")
plot!(time_profile.time, brine_profile.C, label="Brine C")


plot(time_profile.time, permeate_profile.C.*permeate_profile.Q)

plot!(time_profile.time, brine_profile.C.*brine_profile.Q)


sol.u
plot(diff(sol.t), cumsum(diff(sol.t)))
# Concentration
plot(sol.t, [u[1] for u in sol.u], title="Concentration")
# Flow rate
plot(sol.t, [u[2] for u in sol.u], title="Flowrate")

# Pressure
plot(sol.t, [u[3] for u in sol.u], title="Pressure")

raw_feeds           = [feed for _ in sol.t]
# 파이프의 부피가 작으면 농축수의 농도가 거의 즉시 반영이 됨.
plot(sol.t, [u[1]*u[2] for u in sol.u], title="Salt flow")

circulated_brines   = [Water(u[2], feed.T, u[1], u[3]) for u in sol.u]
adjusted_feeds      = [Water(1000/84-brine.Q, feed.T, 300/1000, 1e-5) for brine in circulated_brines]
circulated_brine_profile = profile_water(circulated_brines);
adjusted_feed_profile = profile_water(adjusted_feeds);

final_feeds = [mix(row[1], row[2]) for row in eachslice(hcat(circulated_brines, adjusted_feeds), dims=1)]
final_feed_profile = profile_water(final_feeds);

plot(sol.t.*u"s", final_feed_profile.C)

sixty_sec_final_feed_profile = hcat(DataFrame(:time => sol.t.*u"s"), final_feed_profile)


plot(sixty_sec_final_feed_profile.time, sixty_sec_final_feed_profile.C .* sixty_sec_final_feed_profile.Q)
# plot!(sol.t.*u"s", cumsum(adjusted_feed_profile.C .* adjusted_feed_profile.Q))
plot!(sol.t.*u"s", (adjusted_feed_profile.C .* adjusted_feed_profile.Q .* sol.t))

plot(sixty_sec_final_feed_profile.time, sixty_sec_final_feed_profile.Q)

@show sol.t;