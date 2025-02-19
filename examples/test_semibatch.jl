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

pipe_volume = 0.05
pipe = CirculationPipe(pipe_volume)
membrane_spec = (25, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
vessel = PressureVessel(7, membrane)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(vessel))

feed_T = 20.0u"°C"
feed_C = 0.5u"kg/m^3"

flowrate_setpoint = 20.0u"m^3/hr"
pressure_setpoint_low  = 7.0u"bar"
pressure_setpoint_high = 13.0u"bar"

dt      = 60.0u"s"
mode    = :CC
fouling = true

next_u = nothing
result_dfs = []
@showprogress for cycle ∈ 1:5
    global next_u
    for i ∈ 1:80
        pressure = pressure_setpoint_high*i/80 + pressure_setpoint_low*(80-i)/80
        # pressure = pressure_setpoint_high
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure;
            process=sbro, dt=dt, mode=:CC, fouling=fouling
        )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
        # push!(result_dfs, last(result_df,1))
    end
    
    for i ∈ 1:1
        pressure = pressure_setpoint_high
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure;
            process=sbro, dt=dt, mode=:purge, fouling=fouling
            )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
        # push!(result_dfs, last(result_df,1))
    end
end

total_df        = reduce(vcat, result_dfs)
total_df.time   = cumsum(total_df.time_diff)

plotting_time_unit = u"hr"

begin
    # Permeate, brine, circulated brine Q
    flowrate_plot = plot(size=(750,500), legend=:right, label="Permeate Q", title="Flowrate", framestyle=:box)
    plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_Q, label="Permeate Q", xlabel="Time", ylabel="Flowrate", leftpad=25Plots.mm)
    plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_Q, label="Raw brine Q")
    scatter!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Circulated_brine_Q, label="Circulated brine Q", markersize=1, markerstrokewidth=0.1)

    # Permeate, brine, circulated brine C
    conc_plot = plot(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_C, label="Permeate C", framestyle=:box, title="Concentration", xlabel="Time", ylabel="Concentration")
    plot!(conc_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_C, label="Raw brine C")
    scatter!(conc_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Circulated_brine_C, label="Circulated (disposed) brine C", markersize=1, markerstrokewidth=0.1)

    # Permeate, brine, circulated brine P
    pressure_plot = plot(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_P.|>u"bar", label="Permeate P", framestyle=:box, title="Pressure", xlabel="Time", ylabel="Pressure")
    plot!(pressure_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_P.|>u"bar", label="Raw brine P")
    scatter!(pressure_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Circulated_brine_P.|>u"bar", label="Circulated (disposed) brine P", markersize=1, markerstrokewidth=0.1)

    cumulative_permeate_plot = scatter(uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_Q .* total_df.time_diff .|>u"m^3" |> cumsum, label="Cumulative permeate quantity", framestyle=:box, title="Cumulate permeate quantity", xlabel="Time", ylabel="Quantity", markersize=1, markerstrokewidth=0.1)

    summary_plot = plot(flowrate_plot, conc_plot, pressure_plot, cumulative_permeate_plot, layout=(2,2), size=(950, 800), link=:x)
    display(summary_plot)
    savefig(summary_plot, "Water-80:2-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png")
end

# Cycle-wise recovery
cycle_dfs = []
begin
    grouped_total_df = groupby(total_df, :cycle)
    
    for cycle_df ∈ grouped_total_df
        cycle_df_CC             = copy(filter((x -> (x.mode==:CC)), cycle_df))
        cycle_raw_feed_quantity = cumsum(cycle_df_CC.Raw_feed_Q.*cycle_df_CC.time_diff)[end]
        cycle_permeate_quantity = cumsum(cycle_df_CC.Permeate_Q.*cycle_df_CC.time_diff)
        cycle_recovery          = uconvert.(NoUnits, cycle_permeate_quantity./(cycle_raw_feed_quantity+pipe_volume*u"m^3")) .* 100.0
        cycle_df_CC.recovery    = cycle_recovery
        push!(cycle_dfs, cycle_df_CC)
    end
    cycle_df = reduce(vcat, cycle_dfs)

    recovery_plot = scatter(
        uconvert.(plotting_time_unit,cycle_df.time),
        cycle_df.recovery,label="Recovery",
        framestyle=:box,
        markersize=2, markerstrokewidth=0.1,
        ylabel="Recovery rate (%)", xlabel="Time"
        )
    ylims!(recovery_plot, -5.0, 105.0)
    
    pressure_plot = scatter(
        uconvert.(plotting_time_unit,cycle_df.time),
        cycle_df.Final_feed_P.|>u"bar",label="Feed pressure",
        framestyle=:box,
        markersize=2, markerstrokewidth=0.1,
        ylabel="Pressure", xlabel="Time"
        )
    scatter!(
            pressure_plot, cycle_df.time .|> plotting_time_unit, cycle_df.Raw_brine_P.|>u"bar",label="Brine pressure",
            markersize=2, markerstrokewidth=0.1
        )
    # ylims!(pressure_plot, (0.9*pressure_setpoint_low)|> u"bar" |> ustrip, (1.1*pressure_setpoint_high)|> u"bar" |> ustrip)
    
    merged_plot = plot(recovery_plot, pressure_plot, layout=(2,1), link=:x)
    display(merged_plot)
    savefig(merged_plot, "Recovery-80:2-$(pressure_setpoint_low)-to-$(pressure_setpoint_high).png")
end
cycle_df


# Membrane resistance increase
pristine_vessel_profile = profile_vessel(vessel)
fouled_vessel_profile   = profile_vessel(sbro.pressure_vessel)
resistance_plot = plot(cumsum(pristine_vessel_profile.dx), pristine_vessel_profile.R_m, label="Pristine membrane")
plot!(resistance_plot, cumsum(fouled_vessel_profile.dx), fouled_vessel_profile.R_m, fillrange=pristine_vessel_profile.R_m, fillcolor=:black, fillstyle=:/, label="Fouled membrane")