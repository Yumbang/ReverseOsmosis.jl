module SemiBatch

export
    CirculationPipe, CirculationPipe2, empty_pipe, empty_pipe2, fill_pipe!, update_pipe!,
    SemiBatchRO, process_semibatch_RO!

using DataFrames, Unitful, Printf

using Interpolations: linear_interpolation

using ...ReverseOsmosis: Water, Water2, pressurize, profile_water, mix

using ...ReverseOsmosis: MembraneElement, MembraneElement2, MembraneModule, PressureVessel, foul!

using ...ReverseOsmosis: pump, vessel_filtration


"""
Pipe full of fluid, with `length` and cross-section `area`.
A pipe is considered to be completely-stirred as concentration `C`.

This struct—designed as brine circulation line—is necessary to properly model brine circulation.
"""
mutable struct CirculationPipe
    const n_segments::Int64     # [-]
    const area::Float64         # [m²] 
    const length::Float64       # [m]
    C_array::Vector{Float64}    # [kg/m³]
end

mutable struct CirculationPipe2
    const n_segments_C::Int64   # [-]
    const n_segments_Cf::Int64  # [-]
    const area::Float64         # [m²]
    const length::Float64       # [m]
    C_array::Vector{Float64}    # [kg/m³]
    Cf_array::Vector{Float64}   # [kg/m³]
end

function empty_pipe(n_segments::Int64, area::Float64, length::Float64)
    return CirculationPipe(n_segments, area, length, zeros(Float64, n_segments))
end

function empty_pipe2(n_segments_C::Int64, n_segments_Cf::Int64, area::Float64, length::Float64)
    return CirculationPipe2(n_segments_C, n_segments_Cf, area, length, zeros(Float64, n_segments_C), zeros(Float64, n_segments_Cf))
end

function fill_pipe!(water::Water; pipe::CirculationPipe)
    pipe.C_array .= water.C
    return nothing
end

function fill_pipe!(water::Water2; pipe::CirculationPipe2)
    pipe.C_array  .= water.C
    pipe.Cf_array .= water.Cf
    return nothing
end

function update_pipe!(influent::Water, dt::Float64; pipe::CirculationPipe)::Water
    # Flux calculation
    u                   = influent.Q / 3600 / pipe.area # [m/s]
    dx                  = u * dt                        # [m]

    if dx ≈ 0.0
        println("dx is too small. Skip updating pipe... $(influent)")
        return Water(influent.Q, influent.T, pipe.C_array[end], influent.P, influent.m_avg)
    end

    # Pipe segments' concentration mapping
    len_seg = pipe.length / pipe.n_segments          # [m]
    initial_centroid  = len_seg/2.0                  # [m]
    last_centroid     = pipe.length - len_seg/2.0    # [m]
    centroids         = range(initial_centroid, last_centroid, length=pipe.n_segments)
    shifted_centroids = centroids .+ dx              # [m]
    C                 = influent.C                   # [kg/m³]

    sample_centroids  = vcat([0.0, dx], shifted_centroids)
    sample_concs      = vcat([C, C]   , pipe.C_array)
    
    # Updated segments' concentration calculation
    concentration_map = linear_interpolation(sample_centroids, sample_concs)
    updated_C_array   = concentration_map.(centroids)

    # Effluent water concentration calculation
    eff_segs   = dx/len_seg
    n_eff_segs = floor(Int64, eff_segs)
    remain_seg = eff_segs - n_eff_segs
    
    C_eff = (remain_seg * pipe.C_array[end-n_eff_segs] + sum(pipe.C_array[end-n_eff_segs+1:end])) / eff_segs

    pipe.C_array = updated_C_array

    return Water(influent.Q, influent.T, C_eff, influent.P, influent.m_avg)
end

function update_pipe!(influent::Water2, dt::Float64; pipe::CirculationPipe2)::Water2
    # Flux calculation
    u                   = influent.Q / 3600 / pipe.area # [m/s]
    dx                  = u * dt                        # [m]

    if dx == 0.0
        println("dx is too small. Skip updating pipe... (Influent: $(influent))")
        return Water2(influent.Q, influent.T, pipe.C_array[end], pipe.Cf_array[end], influent.P, influent.m_avg)
    end
    # ======================= C =======================
    # Pipe segments' concentration mapping
    len_seg_C  = pipe.length / pipe.n_segments_C       # [m]
    initial_centroid_C  = len_seg_C/2.0                # [m]
    last_centroid_C     = pipe.length - len_seg_C/2.0  # [m]
    centroids_C         = range(initial_centroid_C, last_centroid_C, length=pipe.n_segments_C)
    shifted_centroids_C = centroids_C .+ dx            # [m]
    C                   = influent.C                   # [kg/m³]

    # Updated segments' concentration calculation
    sample_centroids_C  = vcat([0.0, dx], shifted_centroids_C)
    sample_concs        = vcat([C, C]   , pipe.C_array)
    
    concentration_map = linear_interpolation(sample_centroids_C, sample_concs)
    updated_C_array   = concentration_map.(centroids_C)

    # Effluent water concentration calculation
    eff_segs_C   = dx/len_seg_C
    n_eff_segs_C = floor(Int64, eff_segs_C)
    remain_seg_C = eff_segs_C - n_eff_segs_C
    
    C_eff        = (remain_seg_C * pipe.C_array[end-n_eff_segs_C] +
                     sum(pipe.C_array[end-n_eff_segs_C+1:end])) / eff_segs_C

    pipe.C_array = updated_C_array

    # ======================= Cf =======================
    # Pipe segments' concentration mapping
    len_seg_Cf  = pipe.length / pipe.n_segments_Cf      # [m]
    initial_centroid_Cf  = len_seg_Cf/2.0                  # [m]
    last_centroid_Cf     = pipe.length - len_seg_Cf/2.0    # [m]
    centroids_Cf         = range(initial_centroid_Cf, last_centroid_Cf, length=pipe.n_segments_Cf)
    shifted_centroids_Cf = centroids_Cf .+ dx              # [m]
    Cf                   = influent.Cf                   # [kg/m³]

    # Updated segments' concentration calculation
    sample_centroids_Cf  = vcat([0.0, dx], shifted_centroids_Cf)
    sample_concs_f       = vcat([Cf, Cf]   , pipe.Cf_array)
    
    concentration_f_map = linear_interpolation(sample_centroids_Cf, sample_concs_f)
    updated_Cf_array    = concentration_f_map.(centroids_Cf)

    # Effluent water concentration calculation
    eff_segs_Cf   = dx/len_seg_Cf
    n_eff_segs_Cf = floor(Int64, eff_segs_Cf)
    remain_seg_Cf = eff_segs_Cf - n_eff_segs_Cf
    
    Cf_eff = (remain_seg_Cf * pipe.Cf_array[end-n_eff_segs_Cf] + sum(pipe.Cf_array[end-n_eff_segs_Cf+1:end])) / eff_segs_Cf

    pipe.Cf_array = updated_Cf_array

    return Water2(influent.Q, influent.T, C_eff, Cf_eff, influent.P, influent.m_avg)
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
    const n_vessels::Int64
    pipe_to_junction::Union{CirculationPipe, CirculationPipe2}
    pipe_to_membrane::Union{CirculationPipe, CirculationPipe2}
    pressure_vessel::PressureVessel
end


function process_semibatch_RO!(
    fresh_feed::Water, last_raw_brine::Water, flowrate_setpoint::typeof(0.0u"m^3/hr"), pressure_setpoint::typeof(0.0u"bar");
    Δt::Unitful.Time, mode::Symbol, process::SemiBatchRO, fouling::Bool=true
    )
    local last_brine::Water

    dt       = 0.1u"s" # HARD CODED dt
    dt_sec   = ustrip(dt)
    t        = 0.0
    iter_num = (Δt / dt) |> NoUnits |> Int64

    @assert (dt < Δt) "Δt ($(Δt)) must be longer than dt ($(dt))"
    @assert (((Δt*10.0) % (dt*10.0)) ≈ 0.0u"s") "Δt ($(Δt)) must be divisible with dt ($(dt)) <mod $(Δt % dt)>"
    @assert (mode ∈ [:CC, :purge]) "Only :CC and :purge mode are supported ($(mode))"

    flowrate_setpoint_m3h = flowrate_setpoint |> u"m^3/hr" |> ustrip
    pressure_setpoint_Pa  = pressure_setpoint |> u"Pa"     |> ustrip
    
    time_log         = Float64[]
    fresh_feed_log   = Water[]
    junction_eff_log = Water[]
    final_feed_log   = Water[]
    raw_brine_log    = Water[]
    disp_brine_log   = Water[]
    permeate_log     = Water[]
    power_circ_log   = typeof(0.0u"W")[]
    power_HPP_log    = typeof(0.0u"W")[]

    last_brine = last_raw_brine


    for iter ∈ 1:iter_num
        t += dt_sec

        if mode == :CC
            # Simulate raw brine circulation
            junction_effluent = update_pipe!(last_brine, dt_sec; pipe = process.pipe_to_junction)

            # Define fresh feed
            fresh_feed_Q = max(0.0, flowrate_setpoint_m3h - junction_effluent.Q)
            fresh_feed   = Water(fresh_feed_Q, fresh_feed.T, fresh_feed.C, fresh_feed.P, fresh_feed.m_avg)

            # Calculate pressure to apply
            effluent_pressure_to_apply = max(0.0, pressure_setpoint_Pa - junction_effluent.P)
            feed_pressure_to_apply     = max(0.0, pressure_setpoint_Pa - fresh_feed.P)
    
            # Pressurize & calculate power consumption
            pressurized_junction_effluent, power_circ = pump(junction_effluent, effluent_pressure_to_apply;
                                                             process.circulation_pump_efficiency)
            pressurized_fresh_feed, power_HPP         = pump(fresh_feed, feed_pressure_to_apply;
                                                             process.high_pressure_pump_efficiency)
    
            # Mix Water variables.
            # Both water are assumed to be in lower pressure compared to pressure_setpoint. The result may not be reliable otherwise.
            junction_feed = mix(pressurized_junction_effluent, pressurized_fresh_feed; pressure=pressure_setpoint_Pa)
    
            # Simulate junction feed circulation
            final_feed = update_pipe!(junction_feed, dt_sec; pipe = process.pipe_to_membrane)
        else
            # Raw brine is not circulated
            junction_effluent = Water(0.0, last_brine.T, process.pipe_to_junction.C_array[end], last_brine.P, last_brine.m_avg)

            # Define fresh feed
            fresh_feed_Q = flowrate_setpoint_m3h
            fresh_feed   = Water(fresh_feed_Q, fresh_feed.T, fresh_feed.C, fresh_feed.P, fresh_feed.m_avg)

            # Calculate pressure to apply
            effluent_pressure_to_apply = 0.0
            feed_pressure_to_apply     = max(0.0, pressure_setpoint_Pa - fresh_feed.P)

            # Pressurize & calculate power consumption
            pressurized_junction_effluent, power_circ = pump(junction_effluent, effluent_pressure_to_apply;
                                                            process.circulation_pump_efficiency)
            pressurized_fresh_feed, power_HPP         = pump(fresh_feed, feed_pressure_to_apply;
                                                            process.high_pressure_pump_efficiency)

            # Mix Water variables.
            # Both water are assumed to be in lower pressure compared to pressure_setpoint. The result may not be reliable otherwise.
            junction_feed = pressurized_fresh_feed

            # Simulate junction feed circulation
            final_feed = update_pipe!(junction_feed, dt_sec; pipe = process.pipe_to_membrane)
        end
        
        # Simulate RO process
        # SBRO process is assumed to have n_vessels of vessels in parallel, and we do not discriminate those (∵ computation time).
        split_final_feed = Water(final_feed.Q / process.n_vessels, final_feed.T, final_feed.C, final_feed.P, final_feed.m_avg)
        RO_result        = vessel_filtration(split_final_feed, process.pressure_vessel; dt=dt_sec)

        brines_array    = RO_result[1]
        permeates_array = RO_result[2]
        fouling_info    = RO_result[3] # ΔR_ms_array

        split_brine    = brines_array[end][end]
        raw_brine      = Water(split_brine.Q * process.n_vessels, split_brine.T,
                                split_brine.C, split_brine.P, split_brine.m_avg)

        split_permeate = mix(reduce(vcat, permeates_array))
        final_permeate = Water(split_permeate.Q * process.n_vessels, split_permeate.T,
                                split_permeate.C, split_permeate.P, split_permeate.m_avg)

        if fouling
            foul!(process.pressure_vessel, fouling_info...)
        end

        if mode == :CC
          # raw_brine  = raw_brine
            disp_brine = Water(0.0, raw_brine.T, raw_brine.C, raw_brine.P, raw_brine.m_avg)
        else
            disp_brine = raw_brine
            raw_brine  = Water(0.0, raw_brine.T, raw_brine.C, raw_brine.P, raw_brine.m_avg)
        end

        push!(time_log,         t)
        push!(fresh_feed_log,   fresh_feed)
        push!(junction_eff_log, junction_effluent)
        push!(final_feed_log,   final_feed)
        push!(raw_brine_log,    raw_brine)
        push!(disp_brine_log,   disp_brine)
        push!(permeate_log,     final_permeate)
        push!(power_circ_log,   power_circ)
        push!(power_HPP_log,    power_HPP)
        
        last_brine = raw_brine
    end

    time_profile         = DataFrame(:time => time_log .* u"s")
    fresh_feed_profile   = profile_water(fresh_feed_log)
    junction_eff_profile = profile_water(junction_eff_log)    
    final_feed_profile   = profile_water(final_feed_log)
    raw_brine_profile    = profile_water(raw_brine_log)
    disp_brine_profile   = profile_water(disp_brine_log)
    permeate_profile     = profile_water(permeate_log)
    power_profile        = DataFrame(:HPP_power=>power_HPP_log, :Circ_power=>power_circ_log)

    rename!((x -> "Fresh_feed_"*x),        fresh_feed_profile)
    rename!((x -> "Junction_effluent_"*x), junction_eff_profile)
    rename!((x -> "Final_feed_"*x),        final_feed_profile)
    rename!((x -> "Raw_brine_"*x),         raw_brine_profile)
    rename!((x -> "Disposed_brine_"*x),    disp_brine_profile)
    rename!((x -> "Permeate_"*x),          permeate_profile)

    result_profile = hcat(time_profile,
                          fresh_feed_profile, 
                          junction_eff_profile, 
                          final_feed_profile, 
                          raw_brine_profile, 
                          disp_brine_profile, 
                          permeate_profile, 
                          power_profile)

    result_profile.time_diff .= dt
    result_profile.mode      .= mode

    return result_profile, last_brine
end

function process_semibatch_RO!(
    fresh_feed::Water2, last_raw_brine::Water2, flowrate_setpoint::typeof(0.0u"m^3/hr"), pressure_setpoint::typeof(0.0u"bar");
    # fresh_feed::Water2, last_raw_brine::Water2, last_final_permeate::Water2, pressure_setpoint::typeof(0.0u"bar");
    Δt::Unitful.Time, mode::Symbol, process::SemiBatchRO, fouling::Bool=true
    )
    local last_brine::Water2
    local last_permeate::Water2

    dt       = 0.1u"s" # HARD CODED dt
    dt_sec   = ustrip(dt)
    t        = 0.0
    iter_num = (Δt / dt) |> NoUnits |> Int64


    @assert (dt < Δt) "Δt ($(Δt)) must be longer than dt ($(dt))"
    @assert (((Δt*10.0) % (dt*10.0)) ≈ 0.0u"s") "Δt ($(Δt)) must be divisible with dt ($(dt)) <mod $(Δt % dt)>"
    @assert (mode ∈ [:CC, :purge]) "Only :CC and :purge mode are supported ($(mode))"

    flowrate_setpoint_m3h = flowrate_setpoint |> u"m^3/hr" |> ustrip
    pressure_setpoint_Pa  = pressure_setpoint |> u"Pa"     |> ustrip
    
    time_log         = Float64[]
    fresh_feed_log   = Water2[]
    junction_eff_log = Water2[]
    final_feed_log   = Water2[]
    raw_brine_log    = Water2[]
    disp_brine_log   = Water2[]
    permeate_log     = Water2[]
    power_circ_log   = typeof(0.0u"W")[]
    power_HPP_log    = typeof(0.0u"W")[]
    
    pipe_j_profile_log = []
    pipe_m_profile_log = []

    last_brine    = last_raw_brine
    # last_permeate = last_final_permeate

    for iter ∈ 1:iter_num
        t += dt_sec

        if mode == :CC
            # Simulate raw brine circulation
            junction_effluent = update_pipe!(last_brine, dt_sec; pipe = process.pipe_to_junction)

            # Define fresh feed
            fresh_feed_Q = max(0.0, flowrate_setpoint_m3h - junction_effluent.Q)
            # fresh_feed_Q = last_permeate.Q
            fresh_feed   = Water2(fresh_feed_Q, fresh_feed.T, fresh_feed.C, fresh_feed.Cf, fresh_feed.P, fresh_feed.m_avg)

            # Calculate pressure to apply
            effluent_pressure_to_apply = max(0.0, pressure_setpoint_Pa - junction_effluent.P)
            feed_pressure_to_apply     = max(0.0, pressure_setpoint_Pa - fresh_feed.P)
    
            # Pressurize & calculate power consumption
            pressurized_junction_effluent, power_circ = pump(junction_effluent, effluent_pressure_to_apply;
                                                             efficiency = process.circulation_pump_efficiency)
            pressurized_fresh_feed, power_HPP         = pump(fresh_feed, feed_pressure_to_apply;
                                                             efficiency = process.high_pressure_pump_efficiency)
    
            # Mix Water2 variables.
            # Both water are assumed to be in lower pressure compared to pressure_setpoint. The result may not be reliable otherwise.
            junction_feed = mix(pressurized_junction_effluent, pressurized_fresh_feed; pressure = pressure_setpoint_Pa)
    
            # Simulate junction feed circulation
            final_feed = update_pipe!(junction_feed, dt_sec; pipe = process.pipe_to_membrane)
        else
            # Raw brine is not circulated
            junction_effluent = Water2(0.0, last_brine.T, process.pipe_to_junction.C_array[end], process.pipe_to_junction.Cf_array[end], last_brine.P, last_brine.m_avg)

            # Define fresh feed
            fresh_feed_Q = flowrate_setpoint_m3h
            # fresh_feed_Q = last_permeate.Q
            fresh_feed   = Water2(fresh_feed_Q, fresh_feed.T, fresh_feed.C, fresh_feed.Cf, fresh_feed.P, fresh_feed.m_avg)

            # Calculate pressure to apply
            effluent_pressure_to_apply = 0.0
            feed_pressure_to_apply     = max(0.0, pressure_setpoint_Pa - fresh_feed.P)

            # Pressurize & calculate power consumption
            pressurized_junction_effluent, power_circ = pump(junction_effluent, effluent_pressure_to_apply;
                                                            efficiency = process.circulation_pump_efficiency)
            pressurized_fresh_feed, power_HPP         = pump(fresh_feed, feed_pressure_to_apply;
                                                            efficiency = process.high_pressure_pump_efficiency)

            # Mix Water2 variables.
            # Both water are assumed to be in lower pressure compared to pressure_setpoint. The result may not be reliable otherwise.
            junction_feed = pressurized_fresh_feed

            # Simulate junction feed circulation
            final_feed = update_pipe!(junction_feed, dt_sec; pipe = process.pipe_to_membrane)
        end
        
        # Simulate RO process
        # SBRO process is assumed to have n_vessels of vessels in parallel, and we do not discriminate those (∵ computation time).
        split_final_feed = Water2(final_feed.Q / process.n_vessels, final_feed.T, final_feed.C, final_feed.Cf, final_feed.P, final_feed.m_avg)
        RO_result        = vessel_filtration(split_final_feed, process.pressure_vessel; dt=dt_sec)

        brines_array    = RO_result[1]
        permeates_array = RO_result[2]
        fouling_info    = RO_result[3:4] # ΔCake_thicknesses_array, ΔR_cs_array

        split_brine    = brines_array[end][end]
        raw_brine      = Water2(split_brine.Q * process.n_vessels, split_brine.T,
                                split_brine.C, split_brine.Cf, split_brine.P, split_brine.m_avg)

        split_permeate = mix(reduce(vcat, permeates_array))
        final_permeate = Water2(split_permeate.Q * process.n_vessels, split_permeate.T,
                                split_permeate.C, split_permeate.Cf, split_permeate.P, split_permeate.m_avg)

        if fouling
            foul!(process.pressure_vessel, fouling_info...)
        end

        if mode == :CC
          # raw_brine  = raw_brine
            disp_brine = Water2(0.0, raw_brine.T, raw_brine.C, raw_brine.Cf, raw_brine.P, raw_brine.m_avg)
        else
            disp_brine = raw_brine
            raw_brine  = Water2(0.0, raw_brine.T, raw_brine.C, raw_brine.Cf, raw_brine.P, raw_brine.m_avg)
        end

        push!(time_log,         t)
        push!(fresh_feed_log,   fresh_feed)
        push!(junction_eff_log, junction_effluent)
        push!(final_feed_log,   final_feed)
        push!(raw_brine_log,    raw_brine)
        push!(disp_brine_log,   disp_brine)
        push!(permeate_log,     final_permeate)
        push!(power_circ_log,   power_circ)
        push!(power_HPP_log,    power_HPP)

        push!(pipe_j_profile_log, process.pipe_to_junction.C_array)
        push!(pipe_m_profile_log, process.pipe_to_membrane.C_array)
        
        last_brine    = raw_brine
        # last_permeate = final_permeate
    end

    time_profile         = DataFrame(:time => time_log .* u"s")
    fresh_feed_profile   = profile_water(fresh_feed_log)
    junction_eff_profile = profile_water(junction_eff_log)    
    final_feed_profile   = profile_water(final_feed_log)
    raw_brine_profile    = profile_water(raw_brine_log)
    disp_brine_profile   = profile_water(disp_brine_log)
    permeate_profile     = profile_water(permeate_log)
    power_profile        = DataFrame(:HPP_power=>power_HPP_log, :Circ_power=>power_circ_log)

    rename!((x -> "Fresh_feed_"*x),        fresh_feed_profile)
    rename!((x -> "Junction_effluent_"*x), junction_eff_profile)
    rename!((x -> "Final_feed_"*x),        final_feed_profile)
    rename!((x -> "Raw_brine_"*x),         raw_brine_profile)
    rename!((x -> "Disposed_brine_"*x),    disp_brine_profile)
    rename!((x -> "Permeate_"*x),          permeate_profile)

    result_profile = hcat(time_profile,
                          fresh_feed_profile, 
                          junction_eff_profile, 
                          final_feed_profile, 
                          raw_brine_profile, 
                          disp_brine_profile, 
                          permeate_profile, 
                          power_profile)

    result_profile.time_diff .= dt
    result_profile.mode      .= mode

    # return result_profile, last_brine, last_permeate
    # return result_profile, last_brine
    return result_profile, last_brine, pipe_j_profile_log, pipe_m_profile_log
end
end # SemiBatchODE

# using Plots
# using BenchmarkTools
# using Base.FastMath
# using ReverseOsmosis

# # Pipes construction
# n_seg_C   = 100
# n_seg_Cf  = 10

# A_j = π * (40.0u"mm"/2)^2 |> u"m^2" |> ustrip
# L_j = 2.0u"m" |> ustrip
# pipe_j = empty_pipe2(n_seg_C, n_seg_Cf, A_j, L_j)

# A_m = π * (25.0u"mm"/2)^2 |> u"m^2" |> ustrip
# L_m = 5.5u"m" |> ustrip
# pipe_m = empty_pipe2(n_seg_C, n_seg_Cf, A_m, L_m)

# fresh_water = Water2(1.0, 20.0, 25e-3, 1e-5, 1e5)
# fill_pipe!(fresh_water; pipe = pipe_j)
# fill_pipe!(fresh_water; pipe = pipe_m)

# # Pressure vessel construction
# begin
#     membrane_segments   = 100
#     spec_membrane_channel_H = 34.0u"mil" |> u"m" |> ustrip
#     spec_membrane_Area  = 78.0u"ft^2"  |> u"m^2" |> ustrip
#     spec_membrane_Len   = 40.0u"inch"  |> u"m"   |> ustrip
#     spec_membrane_Width = spec_membrane_Area / spec_membrane_Len

#     spec_permeate_Q     = 9.1u"m^3/d"
#     spec_membane_Area   = 78.0u"ft^2"
#     spec_water_C        = 2000.0u"mg/L"
#     temp_water_Cf       = 1e-2u"mg/L"
#     spec_applied_P      = 15.5u"bar"
#     spec_dummy_water    = Water(1e-10, 25.0, spec_water_C |> u"kg/m^3" |> ustrip, 1e5)
#     spec_dummy_water2   = Water2(1e-10, 25.0, spec_water_C |> u"kg/m^3" |> ustrip, temp_water_Cf |> u"kg/m^3" |> ustrip, 1e5)
#     spec_net_P          = spec_applied_P - osmo_p(spec_dummy_water) * u"Pa" |> u"bar"
#     membrane_A          = spec_permeate_Q / (spec_membane_Area * spec_net_P) |> u"L/m^2/hr/bar"
#     membrane_B          = 0.75625u"L/m^2/hr" |> u"m/s" |> ustrip
#     membrane_R_m        = 1 / membrane_A |> u"Pa*s/m" |> ustrip
#     spec_water_μ        = 2.414e-5 * 10^(247.8 / (20.0 + 273.15 - 140)) * u"Pa*s"
#     membrane_R_m2       = 1 / membrane_A / spec_water_μ |> u"m^-1" |> ustrip

#     k_fp = 0.6741807585817046e9
#     K    = 16.0

#     n_modules = 2
#     average_feed_molar_mass = 44.48828244517122u"g/mol"
#     EC_TO_TDS = 0.6 * 1e-3 * u"kg/m^3/μS*cm"
# end

# membrane         = pristine_membrane2(25, spec_membrane_Len, spec_membrane_channel_H, spec_membrane_Width, membrane_R_m2/2.0, membrane_B/1.6, 0.34, 16.0*3.2)
# vessel           = PressureVessel(2, membrane)

# η_HPP     = 0.8
# η_circ    = 0.8
# n_vessels = 4

# sbro = SemiBatchRO(η_HPP, η_circ, n_vessels, deepcopy(pipe_j), deepcopy(pipe_m), vessel)

# Δt = 30.0u"minute"
# t  = 0.0u"s"

# sbro_result = process_semibatch_RO!(fresh_water, fresh_water, 4.5u"m^3/hr", 5.0u"bar";
#                                     Δt = Δt, mode = :CC, process = sbro, fouling = true)

# # @benchmark sbro_result = process_semibatch_RO!($fresh_water, $fresh_water, 4.5u"m^3/hr", 5.0u"bar";
# #                                                 Δt = 10.0u"s", mode = :CC, process = $sbro, fouling = true)   # Time  (mean ± σ):   4.687 ms ± 935.736 μs
# # @benchmark sbro_result = process_semibatch_RO!($fresh_water, $fresh_water, 4.5u"m^3/hr", 5.0u"bar";
# #                                                 Δt = 10.0u"s", mode = :CC, process = $sbro, fouling = true)   # Time  (mean ± σ):   4.631 ms ± 849.650 μs ) (FastMath version)


# result_df = sbro_result[1]
# last_brine = sbro_result[2]

# plot(result_df.time, result_df.Final_feed_Q)
# plot!(result_df.time, result_df.Permeate_Q)

# plot(sbro.pipe_to_junction.C_array)
# plot!(sbro.pipe_to_membrane.C_array)
# plot!(pipe_m.C_array)

# # @benchmark update_pipe!(new_feed, dt; pipe=pipe)           # 2.643 μs per an evaluation on average
# # @benchmark @fastmath update_pipe!(new_feed, dt; pipe=pipe) # 2.574 μs per an evaluation on average

# feed_log = []
# eff_log  = []
# salt_pipe_log = []

# anim = @animate for i in 1:1000
#     global t
#     global new_feed
#     if i % 25 == 0
#         new_feed = Water2(1.0*rand(), 20.0, 300e-3*rand(), 1e-5, 1e5)
#     end
#     push!(feed_log, new_feed)
#     plot(pipe.C_array, label="Solute concentration", 
#          xlabel="Segment number", ylabel="Concentration (kg/m³)",
#          title="t = $(round(t, digits=1)) s", ylims=[0.0, 300e-3])
#     effluent = update_pipe!(new_feed, dt; pipe=pipe)
#     salt_pipe = ([c for c ∈ pipe.C_array] |> mean) * A * L
#     push!(salt_pipe_log, salt_pipe)
#     push!(eff_log, effluent)
#     t += dt
# end

# mp4(anim, "pipe_concentration.mp4", fps=30)

# salt_feed_log = [feed.C * feed.Q * dt / 3600.0 for feed ∈ feed_log]
# salt_eff_log  = [eff.C * eff.Q * dt / 3600.0 for eff ∈ eff_log]

# diff_salt_pipe = vcat(0.0, diff(salt_pipe_log))

# salt_feed_plot  = plot(salt_feed_log)
# salt_eff_plot   = plot(salt_eff_log)
# salt_pipe_plot  = plot(diff_salt_pipe)
# salt_total_plot = plot((salt_feed_log .- salt_eff_log .- diff_salt_pipe)[2:end])

# plot(salt_feed_plot, salt_eff_plot, salt_pipe_plot, salt_total_plot, layout=(2,2), link=:all)
# feed_Q = [feed.Q for feed ∈ feed_log]
# eff_Q  = [eff.Q for eff ∈ eff_log]
# scatter(feed_Q, eff_Q)
# using Interpolations

# xs = 0:0.5:100
# ys = xs.^2

# lin_interp = linear_interpolation(xs, ys)   # 150ns per evaluation on average
# new_xs = xs .+ 0.25
# new_ys = lin_interp.($new_xs[1:end-1])      # 700ns per evaluation on average

# scatter(xs, ys, label="Original")
# scatter!(new_xs[1:end-1], new_ys, label="Interpolated")

# len_seg = 0.001
# dx = 0.0276

