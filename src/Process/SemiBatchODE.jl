module SemiBatchODE

export
    SemiBatchRO, Pipe,
    process_semibatch_RO!

using DataFrames, Unitful, Printf

using DifferentialEquations

using ...ReverseOsmosis: Water, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, foul!

using ...ReverseOsmosis: pump, vessel_filtration


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

# """
# Model partial mixing behavior happening inside pipe. A pipe is considered to be completely-stirred.
# Without considering pipe's mixing-behavior, one will experience flushing phenomena that is unreasonably sensitive to `dt`.
# """
# function mix_pipe!(circulated_brine::Water, pipe::Pipe; dt_sec::Float64)::Water
#     local vol_displacement::Float64
#     # If the displacement volume is higher than pipe volume, unintended behavior is very likely to happen in brine mixing phase, also demonstrating sensitivity to `dt`.
#     if circulated_brine.Q * dt_sec * 3600 > pipe.volume
#         @warn @sprintf("Fluid displacement (.2f) in pipe is larger than pipe's volume (.2f). The result is not reliable.", feed.Q * dt_sec, pipe.volume)
#         vol_displacement = pipe.volume
#     else
#         vol_displacement = circulated_brine.Q * dt_sec * 3600
#     end

#     C_pipe_updated   = ((pipe.C * (pipe.volume - vol_displacement)) + (circulated_brine.C * vol_displacement)) / pipe.volume
#     effluent         = Water(circulated_brine.Q, circulated_brine.T, C_pipe_updated, circulated_brine.P)
#     pipe.C           = C_pipe_updated
    
#     return effluent
# end

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

    # Unpressurized brine and feed
    unpressurized_brine = Water(u[2], unpressurized_feed.T, u[1], u[3])

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
    # final_permeate = mix(permeates_array)
    
    # ODE terms
    du[1] = final_brine.Q - u[1]
    du[2] = (final_brine.C * final_brine.Q - u[1] * u[2]) / sbro.pipe.volume
    du[3] = final_brine.P - u[3]

    return nothing
end

end # SemiBatchODE