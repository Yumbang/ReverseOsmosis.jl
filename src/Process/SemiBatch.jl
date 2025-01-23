module SemiBatch

export
    SemiBatchRO, Pipe,
    process_semibatch_RO!

using DataFrames, Unitful, Printf

using ...ReverseOsmosis: Water, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, foul!

using ...ReverseOsmosis: pump, vessel_filtration


"""
Pipe full of fluid, with `length` and cross-section `area`.
A pipe is considered to be completely-stirred as concentration `C`.

This struct `Pipe`—designed to be used for brine circulation line—is necessary to properly model brine circulation.
"""
mutable struct Pipe
    const volume::Float64       # [m³]
    C::Union{Float64,Nothing}   # [kg/m³]
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
    pressure_vessel::PressureVessel
    pipe::Pipe
    circulated_brine::Union{Water,Nothing}
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
Model partial mixing behavior happening inside pipe. A pipe is considered to be completely-stirred.
Without considering pipe's mixing-behavior, one will experience flushing phenomena that is unreasonably sensitive to `dt`.
"""
function mix_pipe!(circulated_brine::Water, pipe::Pipe; dt_sec::Float64)::Water
    local vol_displacement::Float64
    # If the displacement volume is higher than pipe volume, unintended behavior is very likely to happen in brine mixing phase, also demonstrating sensitivity to `dt`.
    if circulated_brine.Q * dt_sec * 3600 > pipe.volume
        @warn @sprintf("Fluid displacement (.2f) in pipe is larger than pipe's volume (.2f). The result is not reliable.", feed.Q * dt_sec, pipe.volume)
        vol_displacement = pipe.volume
    else
        vol_displacement = circulated_brine.Q * dt_sec * 3600
    end

    C_pipe_updated   = ((pipe.C * (pipe.volume - vol_displacement)) + (circulated_brine.C * vol_displacement)) / pipe.volume
    effluent         = Water(circulated_brine.Q, circulated_brine.T, C_pipe_updated, circulated_brine.P)
    pipe.C           = C_pipe_updated
    
    return effluent
end

# TODO: Create constructor function to take `Unitful` process configuration variables (membrane spec, module number, etc.) as input.
# TODO: Re-consider position / keyword arguments. `mode` and `profile_process` don't necessarily seem to be keyword arguments.
"""
Process feed water with semi-batch RO process given.
dt is mandatory, and required to be reasonably short, `Unitful` time (e.g. dt = 1.0u"s"). 
""" 
function process_semibatch_RO!(
    unpressurized_feed::Water, pressure_setpoint::Unitful.Pressure;
    process::SemiBatchRO, dt::Unitful.Time=1u"s", mode::Symbol=:CC, fouling::Bool=true, profile_process::Bool=false
)
    @assert (mode ∈ [:CC, :purge]) "Only closed-circuit or purge mode is supported in semi-batch RO process."

    pressure_setpoint_Pa = Floar64(uconvert(u"Pa", pressure_setpoint)/u"Pa")
    if unpressurized_feed.P ≥ pressure_setpoint_Pa
        @warn "Feed pressure is higher than pressure setpoint. Setpoint is adjusted and the result is not reliable."
        pressure_setpoint_Pa = unpressurized_feed.P
    end

    local brines_array::Array{Array{Water}}
    local permeates_array::Array{Array{Water}}
    local ΔR_ms_array::Array{Array{Float64}}

    local pressurized_feed::Water
    local pressurized_brine::Water
    local final_feed::Water
    local final_power_consumption::Float64

    η_HPP   = process.high_pressure_pump_efficiency
    η_circ  = process.circulation_pump_efficiency
    l_junc  = 0.0
    dt_sec  = Float64(uconvert(u"s", dt)/u"s")

    pressure_to_apply_HPP = pressure_setpoint_Pa - unpressurized_feed.P
    pressurized_feed, power_consumption_HPP = pump(unpressurized_feed, pressure_to_apply_HPP; efficiency=η_HPP)
    
    if isnothing(process.circulated_brine) && (mode==:CC)
        # If it is first timestep (no circulated_brine defined), fill pipe and circulated brine with feed (assuming low-pressure filling phase).
        if isnothing(process.pipe.C)
            process.circulated_brine = pressurized_feed
            process.pipe.C           = pressurized_feed.C
        # If it is right after purging phase, 'flash-simulate' process to determine circulated brine while leaving pipe same.
        # Because of this mechanism, short dt is highly recommended.
        else
            brines_array, _, _ = vessel_filtration(final_feed, process.pressure_vessel; dt=dt_sec)
            process.circulated_brine = brines_array[end][end]
        end
    end
    
    # Mix or flow-through depending on operation mode.
    # If either mode is `:purge` flow-through the feed.
    # Else, pressurize feed and circulated brine into pressure setpoint and merge.
    if mode==:purge
        final_feed = pressurized_feed
        final_power_consumption = power_consumption_HPP
    else
        pressure_to_apply_circ = max(0.0, pressure_setpoint_Pa - process.circulated_brine)
        if process.circulated_brine ≥ pressure_setpoint_Pa
            @warn "Brine pressure is higher than pressure setpoint. Exceeding pressure is neglected and the result is not reliable."
        end
        pressurized_brine, power_consumption_circ = pump(process.circulated_brine, pressure_to_apply_circ; efficiency=η_circ)
        final_feed = junction(pressurized_feed, pressurized_brine; junction_loss=l_junc)
        final_power_consumption = power_consumption_HPP + power_consumption_circ
    end

    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
                                                    final_feed, process.pressure_vessel;
                                                    dt=dt_sec
                                                )
    if fouling
        foul!(process.pressure_vessel, ΔR_ms_array)
    end

    mixed_brine     = mix_pipe!(brines_array[end][end], process.pipe; dt_sec = dt_sec)

    final_brine     = mixed_brine
    final_permeate  = mix(reduce(vcat, permeates_array); pressure=1e-5)


    process.circulated_brine = mixed_brine

    if profile_process
        brines_profile      = vcat([profile_water(brines) for (i, brines) in enumerate(brines_array)]...)
        permeates_profile   = vcat([profile_water(permeates) for (i, permeates) in enumerate(permeates_array)]...)
    else
        brines_profile      = [nothing]
        permeates_profile   = [nothing]
    end

    return final_brine, final_permeate, brines_profile, permeates_profile, power_consumption
end

end # SemiBatch