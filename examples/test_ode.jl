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
Pipe full of fluid, with `length` and cross-section `area`.
A pipe is considered to be completely-stirred as concentration `C`.

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
    time::Float64
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
    unpressurized_feed  = p.unpressurized_feed
    pressure_setpoint   = p.pressure_setpoint
    sbro                = p.process
    logger              = p.logger

    # Unpressurized brine and feed
    unpressurized_brine = Water(u[2], unpressurized_feed.T, u[1], u[3])

    # Pressurize feed and brine
    pressurized_feed,  power_feed = pump(unpressurized_feed, pressure_setpoint - unpressurized_feed.P; efficiency=0.8)
    pressurized_brine, power_circ = pump(unpressurized_brine, pressure_setpoint - unpressurized_brine.P; efficiency=0.8)

    # Final feed after mixing
    final_feed = junction(pressurized_brine, pressurized_feed; junction_loss=0.0)
    # println("Feed:  $(final_feed)")
    # RO modeling
    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    final_brine    = brines_array[end][end]
    final_permeate = mix(reduce(vcat, permeates_array); pressure=1e-5)

    # present_log = SBROLoggerElement(final_brine, final_permeate, t)
    # push!(logger, present_log)
    
    # ODE terms
    du[1] = (final_brine.C * final_brine.Q - u[1] * final_brine.Q) / sbro.pipe.volume
    du[2] = final_brine.Q - u[2]
    du[3] = final_brine.P - u[3]

    # u[2]  = final_brine.Q
    # u[3]  = final_brine.P
    # du[2] = 0.0
    # du[3] = 0.0

    return nothing
end

logger = []
pipe = Pipe(1.0)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(pristine_vessel))
feed = Water(100/84, 20.0, 300/1000, 1e5);
# u₀ = [feed.C, feed.Q, feed.P]
u₀ = [feed.C, 0.0, 0.0]
params = (unpressurized_feed=feed, pressure_setpoint=4e5, process=sbro, logger=logger)
prob = ODEProblem(SBRO!, u₀, (0.0, 60.0), params)

@benchmark sol = solve(prob, Tsit5())

sol.u
sol.t
# Concentration
plot(sol.t, [u[1] for u in sol.u], title="Concentration")
# Flow rate
plot(sol.t, [u[2] for u in sol.u], title="Flowrate")

# Pressure
plot(sol.t, [u[3] for u in sol.u], title="Pressure")

raw_feeds           = [feed for _ in sol.t]
# 파이프의 부피가 작으면 농축수의 농도가 거의 즉시 반영이 됨.
plot(sol.t, [u[1]*u[2] for u in sol.u], title="Salt flow")
plot!(sol.t, feed.C*feed.Q*sol.t)

circulated_brines   = [Water(u[2], feed.T, u[1], u[3]) for u in sol.u]
mixed_feeds         = [mix(feed, brine) for brine in circulated_brines]
@show profile_water(circulated_brines);
@show profile_water(mixed_feeds);


@show sol.t;
@show filtered_time = [log_element for log_element in logger if log_element.time in sol.t];