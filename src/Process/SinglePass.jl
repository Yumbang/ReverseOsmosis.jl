"""
Simple single-pass reverse osmosis process structure.
Feed is pressurized by a high-pressure pump and filtered.
"""
mutable struct SinglePassRO
    const pump_efficiency::Float64
    pressure_vessel::PressureVessel
end

# TODO: Create constructor function to take `Unitful` process configuration variables (membrane spec, module number, etc.) as input.
# TODO: Reinforce docstrings to process_feed! function.
"""
Process feed water with single-pass RO process given.
Currently, dt is mandatory, and required to be `Unitful` time (e.g. dt = 1.0u"hr"). 
"""
function process_feed!(
    unpressurized_feed::Water, pressure_setpoint::Unitful.Pressure;
    process::SinglePassRO, dt::Unitful.Time, mode::Symbol=:forward, fouling::Bool=true, profile_process::Bool=false
)
    @assert mode==:forward "Only forward mode is supported in single-pass RO process."

    pressure_setpoint_Pa = Float64(uconvert(u"Pa", pressure_setpoint)/u"Pa")
    if unpressurized_feed.P ≥ pressure_setpoint_Pa
        @warn "Feed pressure is higher than pressure setpoint."
        pressure_setpoint_Pa = unpressurized_feed.P
    end

    local brines_array::Array{Array{Water}}
    local permeates_array::Array{Array{Water}}
    local ΔR_ms_array::Array{Array{Float64}}
    η_pump  = process.pump_efficiency
    dt_sec  = Float64(uconvert(u"s", dt)/u"s")

    # Feed pressurization and power consumption eseimation. `power_consumption` is `Unitful` power.
    pressure_to_apply = pressure_setpoint_Pa - unpressurized_feed.P
    pressurized_feed, power_consumption = pump(unpressurized_feed, pressure_to_apply; efficiency = η_pump)

    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
                                                    pressurized_feed, process.pressure_vessel;
                                                    dt=dt_sec
                                                )

    if fouling
        foul!(process.pressure_vessel, ΔR_ms_array)
    end

    final_brine     = brines_array[end][end]
    final_permeate  = mix(reduce(vcat, permeates_array); pressure=1e-5)

    if profile_process
        brines_profile      = vcat([profile_water(brines) for (i, brines) in enumerate(brines_array)]...)
        permeates_profile   = vcat([profile_water(permeates) for (i, permeates) in enumerate(permeates_array)]...)
    else
        brines_profile      = [nothing]
        permeates_profile   = [nothing]
    end

    return final_brine, final_permeate, brines_profile, permeates_profile, power_consumption
end