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

pipe = CirculationPipe(0.23)
membrane_spec = (25, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
vessel = PressureVessel(7, membrane)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(vessel))

feed_T = 20.0u"°C"
feed_C = 1.0u"kg/m^3"

flowrate_setpoint = 20.0u"m^3/hr"
pressure_setpoint = 10.0u"bar"

dt      = 60.0u"s"
mode    = :CC
fouling = true

next_u = nothing
result_dfs = []
@showprogress for cycle ∈ 1:5
    global next_u
    for i ∈ 1:80
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure_setpoint;
            process=sbro, dt=dt, mode=mode, fouling=fouling
        )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
    end
    
    for i ∈ 1:2
        result_df, next_u = process_semibatch_RO!(
            feed_T, feed_C, next_u,
            flowrate_setpoint, pressure_setpoint;
            process=sbro, dt=dt, mode=:purge, fouling=fouling
            )
        result_df.cycle .= cycle
        push!(result_dfs, result_df)
    end
end

total_df        = reduce(vcat, result_dfs)
total_df.time   = cumsum(total_df.time_diff)

plotting_time_unit = u"d"

flowrate_plot = plot(size=(750,500), legend=:right, label="Permeate Q", title="Flowrate", framestyle=:box)
plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Permeate_Q, label="Permeate Q", xlabel="Time", ylabel="Flowrate", leftpad=25Plots.mm)
plot!(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_Q, label="Raw brine Q")
scatter(flowrate_plot, uconvert.(plotting_time_unit,total_df.time), total_df.Circulated_brine_Q, label="Circulated brine Q", markersize=1, markerstrokewidth=0.1)

plot(uconvert.(plotting_time_unit,total_df.time), total_df.Raw_brine_C, label="Raw brine C", framestyle=:box)
plot!(uconvert.(plotting_time_unit,total_df.time), total_df.Circulated_brine_C, label="Circulated (disposed) brine C")

pristine_vessel_profile = profile_vessel(vessel)
fouled_vessel_profile   = profile_vessel(sbro.pressure_vessel)

plot(cumsum(pristine_vessel_profile.dx), pristine_vessel_profile.R_m, label="Pristine membrane")
plot!(cumsum(fouled_vessel_profile.dx), fouled_vessel_profile.R_m, fillrange=pristine_vessel_profile.R_m, fillcolor=:black, fillstyle=:/, label="Fouled membrane")