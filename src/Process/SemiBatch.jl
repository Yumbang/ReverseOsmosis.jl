"""
Simple semi-batch reverse osmosis process structure.
Feed is pressurized by a high-pressure pump and brine is circulated with circulation pump.

SemiBatchRO is 'stateful', it holds the brine information.
If mode is `:CC`, the brine is mixed into feed, and purged out otherwise.

At the first evaluation, circulated_brine is `nothing`. Process is automatically evaluated without brine circulation (`:purge` mode).

Because of the property of semi-batch RO, using long `dt` is discouraged. Increasing the temporal resolution is recommended.
"""
mutable struct SemiBatchRO
    const high_pressure_pump_efficiency::Float64
    const circulation_pump_efficiency::Float64
    pressure_vessel::PressureVessel
    circulated_brine::Union{Water,Nothing}
end

function junction(water1::Water, water2::Water; junction_loss::Float64)
    @assert (0.0 ≤ junction_loss ≤ 1.0) "Junction loss must be between 0 and 1."
    merged_pressure = max(water1.P, water2.P) * (1.0-junction_loss)
    mixed_water     = mix(water1, water2; pressure=merged_pressure)
    return mixed_water
end

# TODO: Create constructor function to take `Unitful` process configuration variables (membrane spec, module number, etc.) as input.
# TODO: Reinforce docstrings to process_feed! function.
"""
Process feed water with semi-batch RO process given.
Currently, dt is mandatory, and required to be `Unitful` time (e.g. dt = 1.0u"minute"). 
"""
function ReverseOsmosisProcesses.process_feed!(
    unpressurized_feed::Water, pressure_setpoint::Unitful.Pressure;
    process::SemiBatchRO, dt::Unitful.Time=1u"minute", mode::Symbol=:CC, fouling::Bool=true, profile_process::Bool=false
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
    
    # Mix or flow-through depending on operation mode.
    # If either mode is `:purge` or it is first timestep (no circulated_brine defined), flow-through the feed.
    # Else, pressurize feed and circulated brine into pressure setpoint and merge.
    if (mode==:purge || isnothing(process.circulated_brine))
        final_feed = pressurized_feed
        final_power_consumption = power_consumption_HPP
    else
        pressure_to_apply_circ = max(0.0, pressure_setpoint_Pa - process.circulated_brine)
        if process.circulated_brine ≥ pressure_setpoint_Pa
            @warn "Brine pressure is higher than pressure setpoint. Exceeding pressure is neglected and the result is not reliable."
        end
        pressurized_brine, power_consumption_circ = pump(process.circulated_brine, pressure_to_apply_circ; efficiency=η_circ)
        # Currently, energy loss in junction is not considered.
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

    final_brine     = brines_array[end][end]
    final_permeate  = mix(reduce(vcat, permeates_array); pressure=1e-5)

    if mode==:purge
        process.circulated_brine = nothing
    else
        process.circulated_brine = final_brine
    end

    if profile_process
        brines_profile      = vcat([profile_water(brines) for (i, brines) in enumerate(brines_array)]...)
        permeates_profile   = vcat([profile_water(permeates) for (i, permeates) in enumerate(permeates_array)]...)
    else
        brines_profile      = [nothing]
        permeates_profile   = [nothing]
    end

    return final_brine, final_permeate, brines_profile, permeates_profile, power_consumption
end