"""
    SinglePass

Module implementing single-pass reverse osmosis process models.

A single-pass RO system is the simplest configuration where:
- Feed → High-pressure pump → Membrane vessel → Permeate + Brine
- Brine is discarded (not recirculated)
- Common in municipal water treatment and seawater desalination

# Exported Types
- `SinglePassRO`: Process structure with pump and membrane vessel
- `process_singlepass_RO!`: Main simulation function (mutates membrane via fouling)
"""
module SinglePass

export
    SinglePassRO,
    process_singlepass_RO!

using DataFrames, Unitful

using ...ReverseOsmosis: Water, Water2, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneElement2, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, pristine_membrane2, foul!

using ...ReverseOsmosis: pump, vessel_filtration

"""
    SinglePassRO

Stateful structure representing a single-pass RO process.

# Fields
- `pump_efficiency::Float64`: High-pressure pump isentropic efficiency (0-1)
- `pressure_vessel::PressureVessel`: Membrane vessel (mutable, accumulates fouling)

# Description
This structure is stateful - the membrane vessel accumulates fouling over time.
Each call to `process_singlepass_RO!` modifies the membrane properties.

# Examples
```julia
vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
process = SinglePassRO(0.8, vessel)

# Simulate over time
for t in 1:100
    feed = Water(10.0, 25.0, 35.0, 1e5)
    brine, permeate, _, _, power = process_singlepass_RO!(
        feed, 55e5u"Pa"; process=process, dt=3600.0u"s"
    )
    # vessel now contains accumulated fouling
end
```
"""
mutable struct SinglePassRO
    const pump_efficiency::Float64
    pressure_vessel::PressureVessel
end

"""
    process_singlepass_RO!(unpressurized_feed, pressure_setpoint; process, dt, mode=:forward, fouling=true, profile_process=false)

Simulate one time step of single-pass RO operation.
Mutates the membrane vessel by applying fouling (if enabled).

# Arguments
- `unpressurized_feed::Union{Water, Water2}`: Feed water at atmospheric pressure
- `pressure_setpoint::Unitful.Pressure`: Target operating pressure (e.g., 55e5u"Pa")

# Keyword Arguments
- `process::SinglePassRO`: Process structure to simulate
- `dt::Unitful.Time`: Time step duration (e.g., 3600.0u"s" for 1 hour)
- `mode::Symbol=:forward`: Operating mode (currently only :forward supported)
- `fouling::Bool=true`: Whether to apply fouling to membranes
- `profile_process::Bool=false`: Whether to return detailed spatial profiles

# Returns
- Tuple of `(final_brine, final_permeate, brines_profile, permeates_profile, power_consumption)`
  - `final_brine`: Concentrate stream from vessel outlet
  - `final_permeate`: Mixed permeate from all elements
  - `brines_profile`: DataFrame of all brines (or [nothing] if profile_process=false)
  - `permeates_profile`: DataFrame of all permeates (or [nothing] if profile_process=false)
  - `power_consumption`: Electrical power consumption [W] (Unitful)

# Side Effects
- **Mutates `process.pressure_vessel`** by applying fouling (if fouling=true)
- Membrane resistance increases each call, affecting future performance

# Examples
```julia
# Basic usage
vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
process = SinglePassRO(0.8, vessel)
feed = Water(10.0, 25.0, 35.0, 1e5)

brine, permeate, _, _, power = process_singlepass_RO!(
    feed, 55e5u"Pa";
    process=process,
    dt=3600.0u"s",
    fouling=true
)

# With detailed profiling
brine, permeate, brine_df, perm_df, power = process_singlepass_RO!(
    feed, 55e5u"Pa";
    process=process,
    dt=3600.0u"s",
    profile_process=true
)
# brine_df contains spatial concentration/pressure profiles
```
"""
function process_singlepass_RO!(
    unpressurized_feed::Water, pressure_setpoint::Unitful.Pressure;
    process::SinglePassRO, dt::Unitful.Time, mode::Symbol=:forward, fouling::Bool=true, profile_process::Bool=false
)
    @assert (mode == :forward) "Only forward mode is supported in single-pass RO process."

    pressure_setpoint_Pa = Float64(uconvert(u"Pa", pressure_setpoint)/u"Pa")
    if unpressurized_feed.P ≥ pressure_setpoint_Pa
        @warn "Feed pressure is higher than pressure setpoint."
        pressure_setpoint_Pa = unpressurized_feed.P
    end

    local brines_array::Array{Array{Water}}
    local permeates_array::Array{Array{Water}}
    local ΔR_ms_array::Vector{Vector{Float64}}
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
        brines_profile      = vcat([profile_water(brines) for brines in brines_array]...)
        permeates_profile   = vcat([profile_water(permeates) for permeates in permeates_array]...)
    else
        brines_profile      = [nothing]
        permeates_profile   = [nothing]
    end

    return final_brine, final_permeate, brines_profile, permeates_profile, power_consumption
end

"""
    process_singlepass_RO!(unpressurized_feed::Water2, pressure_setpoint; ...)

Process single-pass RO with foulant tracking (Water2 version).
Applies both cake thickness and resistance changes to membranes.

# Note
Requires Water2 feed and MembraneElement2 array in vessel.
Returns same outputs as Water version but uses element_filtration2 internally.
"""
function process_singlepass_RO!(
    unpressurized_feed::Water2, pressure_setpoint::Unitful.Pressure;
    process::SinglePassRO, dt::Unitful.Time, mode::Symbol=:forward, fouling::Bool=true, profile_process::Bool=false
)
    @assert (mode == :forward) "Only forward mode is supported in single-pass RO process."

    pressure_setpoint_Pa = Float64(uconvert(u"Pa", pressure_setpoint)/u"Pa")
    if unpressurized_feed.P ≥ pressure_setpoint_Pa
        @warn "Feed pressure is higher than pressure setpoint."
        pressure_setpoint_Pa = unpressurized_feed.P
    end

    local brines_array::Array{Array{Water}}
    local permeates_array::Array{Array{Water}}
    local ΔR_ms_array::Array{Array{Float64}}
    η_pump  = process.pump_efficiency
    dt_sec  = dt |> u"s" |> ustrip

    # Feed pressurization and power consumption eseimation. `power_consumption` is `Unitful` power.
    pressure_to_apply = pressure_setpoint_Pa - unpressurized_feed.P
    pressurized_feed, power_consumption = pump(unpressurized_feed, pressure_to_apply; efficiency = η_pump)

    brines_array, permeates_array, ΔCake_thicknesses_array, ΔR_cs_array = vessel_filtration(
                                                    pressurized_feed, process.pressure_vessel;
                                                    dt=dt_sec
                                                )

    if fouling
        foul!(process.pressure_vessel, ΔCake_thicknesses_array, ΔR_cs_array)
    end

    final_brine     = brines_array[end][end]
    final_permeate  = mix(reduce(vcat, permeates_array); pressure=1e-5)

    if profile_process
        brines_profile      = vcat([profile_water(brines) for brines in brines_array]...)
        permeates_profile   = vcat([profile_water(permeates) for permeates in permeates_array]...)
    else
        brines_profile      = [nothing]
        permeates_profile   = [nothing]
    end

    return final_brine, final_permeate, brines_profile, permeates_profile, power_consumption
end

end # SinglePass