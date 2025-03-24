begin
    using Pkg
    Pkg.activate("/Users/ybang_mac/Documents/병찬/Research/2025/Batch RO/Code/ReverseOsmosis.jl")
    using Revise
    using ReverseOsmosis
    using DataFrames
    using BenchmarkTools
    using Unitful
    using Plots
    using CSV
    using ProgressMeter
end


pipe_to_junction_volume = 0.0026998061866787284 / 4
pipe_to_membrane_volume = 0.007037167544041136  / 4
n_segments  = 20
pipe_to_junction = CirculationPipe(pipe_to_junction_volume, n_segments)
pipe_to_membrane = CirculationPipe(pipe_to_membrane_volume, n_segments)
membrane_spec = (100, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
vessel = PressureVessel(2, membrane)
sbro = SemiBatchRO(0.8, 0.8, pipe_to_junction, pipe_to_membrane, deepcopy(vessel))

feed_T = 14.0u"°C"
feed_C = 45.0u"mg/L"

flowrate_setpoint = 6.0u"m^3/hr" / 4.0
pressure_setpoint_low   = 3.0u"bar"
pressure_setpoint_high  = 3.0u"bar"
pressure_setpoint_purge = 2.0u"bar"

dt      = 10.0u"s"
mode    = :CC
fouling = true



CC_duration    = 30u"minute"
CC_index       = uconvert(NoUnits, CC_duration / dt) |> Int64

purge_duration = 1u"minute"
purge_index    = uconvert(NoUnits, purge_duration / dt) |> Int64

next_u = nothing
result_dfs = []
u_log = []
@showprogress for cycle ∈ 1:3
    global next_u
    for i ∈ 1:CC_index
        pressure = pressure_setpoint_high*i/CC_index + pressure_setpoint_low*(CC_index-i)/CC_index
        # pressure = pressure_setpoint_high
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure;
            process=sbro, dt=dt, mode=:CC, fouling=fouling
        )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
        push!(u_log, copy(next_u))
        # push!(result_dfs, last(result_df,1))
    end
    next_u[1] = flowrate_setpoint |> ustrip
    for i ∈ 1:purge_index
        pressure = pressure_setpoint_purge
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure;
            process=sbro, dt=dt, mode=:purge, fouling=fouling
        )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
        push!(u_log, copy(next_u))
        # push!(result_dfs, last(result_df,1))
    end
    next_u[1] = flowrate_setpoint |> ustrip
end

total_df        = reduce(vcat, result_dfs)
total_df.time   = cumsum(total_df.time_diff)

scatter(u_log[10][end-2*n_segments+1:end-n_segments])
scatter(u_log[10][end-n_segments+1:end])

scatter(u_log[10][end-2*n_segments+1:end])

total_df

plotting_time_unit = u"minute"
save_dir = "/Users/ybang_mac/Documents/병찬/Research/2025/Batch RO/Code/results"

begin
    # Permeate, brine, circulated brine Q
    flowrate_plot = plot(size=(750,500), legend=:right, label="Permeate Q", title="Flowrate", framestyle=:box)
    plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_Q, label="Permeate Q", xlabel="Time", ylabel="Flowrate", leftpad=25Plots.mm)
    plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_Q, label="Raw brine Q")
    scatter!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Junction_brine_Q, label="Circulated brine Q", markersize=1, markerstrokewidth=0.01)

    # Permeate, brine, circulated brine C
    conc_plot = plot(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_C, label="Permeate C", framestyle=:box, title="Concentration", xlabel="Time", ylabel="Concentration", legend=:bottomright)
    plot!(conc_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_C, label="Raw brine C")
    scatter!(conc_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Junction_brine_C, label="Circulated (disposed) brine C", markersize=1, markerstrokewidth=0.01)
    plot!(conc_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Final_feed_C, label="Final feed C")

    # Permeate, brine, circulated brine P
    pressure_plot = plot(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_P.|>u"bar", label="Permeate P", framestyle=:box, title="Pressure", xlabel="Time", ylabel="Pressure")
    scatter!(pressure_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_P.|>u"bar", label="Raw brine P", markersize=1, markerstrokewidth=0.01)
    scatter!(pressure_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Junction_brine_P.|>u"bar", label="Circulated (disposed) brine P", markersize=1, markerstrokewidth=0.01)

    cumulative_permeate_plot = plot(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_Q .* total_df.time_diff .|>u"m^3" |> cumsum, label="Cumulative permeate quantity", framestyle=:box, title="Cumulate permeate quantity", xlabel="Time", ylabel="Quantity")

    summary_plot = plot(flowrate_plot, conc_plot, pressure_plot, cumulative_permeate_plot, layout=(2,2), size=(950, 800), link=:x)
    display(summary_plot)
    savefig(summary_plot, joinpath(save_dir, "Water-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
end


# Cycle-wise recovery
cycle_dfs = []
begin
    grouped_total_df = groupby(total_df, :cycle)
    
    for cycle_df ∈ grouped_total_df
        cycle_df_CC             = copy(filter((x -> (x.mode==:CC)), cycle_df))
        cycle_df_purge          = copy(filter((x -> (x.mode==:purge)), cycle_df))
        cycle_raw_feed_quantity = cumsum(cycle_df_CC.Raw_feed_Q.*cycle_df_CC.time_diff)[end]
        cycle_permeate_quantity = cumsum(cycle_df_CC.Permeate_Q.*cycle_df_CC.time_diff)
        cycle_purged_feed_quantity     = cumsum(cycle_df_purge.Raw_feed_Q.*cycle_df_purge.time_diff)[end]
        cycle_purged_permeate_quantity = cumsum(cycle_df_purge.Permeate_Q.*cycle_df_purge.time_diff)[end]
        cycle_recovery          = uconvert.(NoUnits, (cycle_permeate_quantity.+cycle_purged_permeate_quantity)./(cycle_raw_feed_quantity+cycle_purged_feed_quantity)) .* 100.0
        cycle_df_CC.recovery    = cycle_recovery
        push!(cycle_dfs, cycle_df_CC)
    end
    cycle_df = reduce(vcat, cycle_dfs)

    recovery_plot = scatter(
        uconvert.(plotting_time_unit,cycle_df.time),
        cycle_df.recovery,label="Recovery",
        framestyle=:box,
        markersize=2, markerstrokewidth=0.01,
        ylabel="Recovery rate (%)", xlabel="Time"
        )
    ylims!(recovery_plot, -5.0, 105.0)
    
    pressure_plot = scatter(
        uconvert.(plotting_time_unit,cycle_df.time),
        cycle_df.Final_feed_P.|>u"bar",label="Feed pressure",
        framestyle=:box,
        markersize=2, markerstrokewidth=0.01,
        ylabel="Pressure", xlabel="Time"
        )
    scatter!(
            pressure_plot, cycle_df.time .|> plotting_time_unit, cycle_df.Raw_brine_P.|>u"bar",label="Brine pressure",
            markersize=2, markerstrokewidth=0.01
    )
    scatter!(
            pressure_plot, cycle_df.time .|> plotting_time_unit, osmo_p.(cycle_df.Final_feed_C .|> ustrip, cycle_df.Final_feed_T .|> ustrip) * u"Pa" .|>u"bar",label="Feed osmotic pressure",
            markersize=2, markerstrokewidth=0.01
    )

    differential_pressure_plot = scatter(
        cycle_df.time .|> plotting_time_unit, (cycle_df.Final_feed_P .- cycle_df.Raw_brine_P).|>u"bar",label="Differential pressure",
        markersize=2, markerstrokewidth=0.01
    )
    osmotic_pressure = osmo_p.(cycle_df.Final_feed_C .|> ustrip, cycle_df.Final_feed_T .|> ustrip) .* u"Pa"
    scatter!(
            differential_pressure_plot, cycle_df.time .|> plotting_time_unit, cycle_df.Final_feed_P .- osmotic_pressure .|>u"bar",label="Transmembrane pressure",
            markersize=2, markerstrokewidth=0.01,
            ylabel="Pressure", xlabel="Time"
    )

    # Merging plots altogether    
    merged_plot = plot(recovery_plot, pressure_plot, differential_pressure_plot, layout=(2,2), link=:x, size=(950, 800))
    display(merged_plot)
    savefig(merged_plot, joinpath(save_dir, "Recovery-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
end

# Salt quantity balance
begin
    purge_df = copy(filter((x->(x.mode == :purge)), total_df))
    purge_salt_flux = (purge_df.Raw_feed_Q .* purge_df.Raw_feed_C - purge_df.Circulated_brine_Q .* purge_df.Circulated_brine_C - purge_df.Permeate_Q .* purge_df.Permeate_C) .* purge_df.time_diff
    purge_flux_df = DataFrame(Dict(:time => purge_df.time, :salt_flux => purge_salt_flux, :time_diff => purge_df.time_diff))

    CC_df = copy(filter((x->(x.mode == :CC)), total_df))
    CC_salt_flux    = (CC_df.Raw_feed_Q .* CC_df.Raw_feed_C - CC_df.Permeate_Q .* CC_df.Permeate_C) .* CC_df.time_diff
    CC_flux_df = DataFrame(Dict(:time => CC_df.time, :salt_flux => CC_salt_flux, :time_diff => CC_df.time_diff))

    tmp_df = vcat(purge_flux_df, CC_flux_df)
    sort!(tmp_df, :time)

    salt_quantity_plot = plot(tmp_df.time .|> plotting_time_unit, tmp_df.salt_flux |> cumsum .|> u"kg", xlabel="Time", ylabel="Salt mass")
    display(salt_quantity_plot)
    savefig(salt_quantity_plot, joinpath(save_dir, "Salt_balance-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
end

# Water quantity balance
begin
    purge_df = copy(filter((x->(x.mode == :purge)), total_df))
    purge_water_flux = (purge_df.Raw_feed_Q - purge_df.Circulated_brine_Q - purge_df.Permeate_Q) .* purge_df.time_diff
    purge_flux_df = DataFrame(Dict(:time => purge_df.time, :water_flux => purge_water_flux, :time_diff => purge_df.time_diff))

    CC_df = copy(filter((x->(x.mode == :CC)), total_df))
    CC_water_flux    = (CC_df.Raw_feed_Q - CC_df.Permeate_Q) .* CC_df.time_diff
    CC_flux_df = DataFrame(Dict(:time => CC_df.time, :water_flux => CC_water_flux, :time_diff => CC_df.time_diff))

    tmp_df = vcat(purge_flux_df, CC_flux_df)
    sort!(tmp_df, :time)

    water_quantity_plot = plot(tmp_df.time .|> plotting_time_unit, tmp_df.water_flux |> cumsum .|> u"m^3", xlabel="Time", ylabel="Water volume")
    display(water_quantity_plot)
    savefig(water_quantity_plot, joinpath(save_dir, "Water_balance-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
end

# Membrane resistance increase
begin
    pristine_vessel_profile = profile_vessel(vessel)
    fouled_vessel_profile   = profile_vessel(sbro.pressure_vessel)
    resistance_plot = plot(cumsum(pristine_vessel_profile.dx), pristine_vessel_profile.R_m, label="Pristine membrane")
    plot!(resistance_plot, cumsum(fouled_vessel_profile.dx), fouled_vessel_profile.R_m, fillrange=pristine_vessel_profile.R_m, fillcolor=:black, fillstyle=:\, label="Fouled membrane")
    display(resistance_plot)
    savefig(resistance_plot, joinpath(save_dir, "Membrane_resistance-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
end


# Energy
begin
    energy_plot = plot(total_df.time  .|> plotting_time_unit, total_df.HPP_power .* total_df.time_diff |> cumsum .|> u"kW*hr", label="HPP energy consumption")
    plot!(energy_plot, total_df.time  .|> plotting_time_unit, total_df.Circ_power.* total_df.time_diff |> cumsum .|> u"kW*hr", label="Circulation energy consumption")
    plot!(energy_plot, total_df.time  .|> plotting_time_unit, (total_df.HPP_power + total_df.Circ_power) .* total_df.time_diff |> cumsum .|> u"kW*hr", label="Total energy consumption")
    display(energy_plot)
    savefig(energy_plot, joinpath(save_dir, "Energy_consumption-$(CC_duration):$(purge_duration)-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png"))
    # SEC = ((total_df.HPP_power + total_df.Circ_power) .* total_df.time_diff |> cumsum .|> u"kW*hr") ./ (total_df.Permeate_Q .* total_df.time_diff .|>u"m^3" |> cumsum) .|> u"kW*hr/m^3"
    # plot(total_df.time  .|> plotting_time_unit, SEC)
end