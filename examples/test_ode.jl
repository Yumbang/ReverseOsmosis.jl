using Pkg
Pkg.activate("/Users/ybang_mac/Documents/병찬/Research/2025/Batch RO/Code/ReverseOsmosis.jl")
using ReverseOsmosis

using DataFrames, Unitful, Printf

using DifferentialEquations: ODEProblem, Tsit5, DP5, solve
using DiffEqCallbacks
using Plots
using BenchmarkTools

# membrane_spec = (25, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane_spec = (25, 1.016, 7e-4, 37/1.016, 6.0e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
pristine_vessel = PressureVessel(7, membrane);

"""
Pipe full of fluid, with `volume` m³.
A pipe is considered to be completely-stirred.

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
    final_feed::Water
    power_feed::typeof(1.0u"W")
    power_circ::typeof(1.0u"W")
    time::Float64
end

mutable struct SBROParameters
    feed_concentration::Float64
    feed_temperature::Float64
    feed_pressure::Float64
    flowrate_setpoint::Float64
    pressure_setpoint::Float64
    τ₁::Float64
    τ₂::Float64
    τ₃::Float64
    process::SemiBatchRO
    logger::SBROLoggerElement
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
    # Raw feed info
    feed_concentration  = p.feed_concentration
    feed_temperature    = p.feed_temperature
    feed_pressure       = p.feed_pressure
    # Setpoint parameters
    flowrate_setpoint   = p.flowrate_setpoint
    pressure_setpoint   = p.pressure_setpoint
    # Time lagging constants
    τ₁                  = p.τ₁    
    τ₂                  = p.τ₂
    τ₃                  = p.τ₃    
    # Target process
    sbro                = p.process

    
    # Unpressurized brine and feed
    unpressurized_brine = Water(u[2], feed_temperature, u[1], u[3])
    unpressurized_feed  = Water(max(flowrate_setpoint-u[2], 0), feed_temperature, feed_concentration, feed_pressure)

    # Pressurize feed and brine
    pressurized_feed,  power_feed = pump(unpressurized_feed,  pressure_setpoint - unpressurized_feed.P;  efficiency=0.8)
    pressurized_brine, power_circ = pump(unpressurized_brine, pressure_setpoint - unpressurized_brine.P; efficiency=0.8)

    # Final feed after mixing
    final_feed = junction(pressurized_brine, pressurized_feed; junction_loss=0.0)

    # RO modeling
    brines_array, permeates_array, ΔR_ms_array = vessel_filtration(
        final_feed, sbro.pressure_vessel; dt=nothing
    )
    
    final_brine    = brines_array[end][end]
    final_permeate = mix(reduce(vcat, permeates_array); pressure=1e-5)

    p.logger = SBROLoggerElement(final_brine, final_permeate, final_feed, power_feed, power_circ, t)
    
    # ODE terms
    # Since (Pipe volume / brine Q) is in hour unit, multiply by 3600 to convert into sec unit.
    du[1] = (final_brine.C - u[1]) / (sbro.pipe.volume/final_brine.Q * 3600) / τ₁ 
    du[2] = (final_brine.Q - u[2]) / (sbro.pipe.volume/final_brine.Q * 3600) / τ₂ 
    du[3] = (final_brine.P - u[3]) / (sbro.pipe.volume/final_brine.Q * 3600) / τ₃ 

    return nothing
end

pipe = Pipe(0.23)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(pristine_vessel))
feed = Water(100/84, 20.0, 300/1000, 1e5);
feed_flowrate_setpoint = 20.0
initial_pressure_setpoint = 10e5
logger = SBROLoggerElement(feed, feed, feed, 0.0u"W", 0.0u"W", 0.0)
logger_array = SBROLoggerElement[]

function save_func(u, t, integrator)
    return deepcopy(integrator.p.logger)
end



next_u = nothing
result_log = []
u_log = []
for i ∈ 1:180
    saved_values = SavedValues(typeof(0.0), typeof(logger))
    cb = SavingCallback(save_func, saved_values)

    # 이전 시행의 결과로 초기조건을 맞추기 위한 코드
    if isnothing(next_u)
        # 1e5 Pa: Atmospheric pressure
        u₀ = [feed.C, 0.0, 1e5]
    else
        u₀ = next_u
    end

    flowrate_setpoint = feed_flowrate_setpoint
    
    # pressure_setpoint = initial_pressure_setpoint + i*0.01e5
    pressure_setpoint = initial_pressure_setpoint

    # 순서대로 농도, 수온, (가압 전)압력, 유량 setpoint, 압력 setpoint, 대상 공정, logging용 변수
    params = SBROParameters(
        1000/1000, flowrate_setpoint, 1e-5,
        feed_flowrate_setpoint, pressure_setpoint, 1.0, 1.0, 1.0, sbro, logger
    )
    
    prob = ODEProblem(SBRO!, u₀, (0.0, 10.0), params)

    sol = solve(prob, Tsit5(); callback=cb)

    # 마지막 시행값을 다음 시행 초기값으로 사용하기 위해 저장
    next_u = deepcopy(sol.u[end])
    push!(u_log, deepcopy(sol.u[end]))

    # Logging 결과값 정리
    time_profile = DataFrame(:time=>saved_values.t[2:end]*u"s")
    brine_profile = profile_water(getfield.(saved_values.saveval, :brine)[2:end])
    rename!((x -> "Brine_"*x), brine_profile)
    permeate_profile = profile_water(getfield.(saved_values.saveval, :permeate)[2:end])
    rename!((x -> "Permeate_"*x), permeate_profile)
    final_feed_profile = profile_water(getfield.(saved_values.saveval, :final_feed)[2:end])
    rename!((x -> "Final_feed_"*x), final_feed_profile)
    circulated_brine_profile = profile_water([Water(u[2], feed.T, u[1], u[3]) for u in sol.u[2:end]])
    rename!((x -> "Circulated_brine_"*x), circulated_brine_profile)
    power_profile = DataFrame(
        :power_feed=>getfield.(saved_values.saveval, :power_feed)[2:end],
        :power_circ=>getfield.(saved_values.saveval, :power_circ)[2:end]
    )
    result_profile = hcat(time_profile, permeate_profile, brine_profile, circulated_brine_profile, final_feed_profile)
    result_profile.time_diff = vcat(diff(result_profile.time), [0.0u"s"])
    push!(result_log, result_profile)

end
# for문에서 축적된 결과값을 병합 후 처리
total_result = vcat(result_log...)
total_result.time = cumsum(total_result.time_diff)

# 누적 permeate, 누적 raw feed 유량 계산
cumulative_permeate = cumsum(uconvert.(u"m^3",total_result.Permeate_Q.*total_result.time_diff))
cumulative_raw_feed = cumsum(uconvert.(u"m^3",(feed_flowrate_setpoint*u"m^3/hr".-total_result.Brine_Q).*total_result.time_diff))
begin
    time_data = ustrip.(total_result.time)./60.0
    cum_perm = ustrip.(cumulative_permeate)
    cum_raw = ustrip.(cumulative_raw_feed)

    # Also, prepare other data by stripping units:
    brine_Q = ustrip.(total_result.Brine_Q)
    final_feed_Q = ustrip.(total_result.Final_feed_Q)
    permeate_Q = ustrip.(total_result.Permeate_Q)
    
    brine_C = ustrip.(total_result.Brine_C)
    final_feed_C = ustrip.(total_result.Final_feed_C)
    permeate_Q = ustrip.(total_result.Permeate_Q)
    permeate_C = ustrip.(total_result.Permeate_C)
    final_feed_P = ustrip.(total_result.Final_feed_P)
    
    brine_P = ustrip.(total_result.Brine_P)

    # Compute mass fluxes (C * Q) for brine and final feed.
    mass_flux_brine = ustrip.(total_result.Brine_C .* total_result.Brine_Q)
    mass_flux_final = ustrip.(total_result.Final_feed_C .* total_result.Final_feed_Q)

    # Compute the recovery rate, defined here as:
    # (Permeate_Q * Permeate_C) / (Final_feed_Q * Final_feed_C)
    # (Be sure these quantities are consistent with your intended definition.)
    recovery_rate = ustrip.(cumulative_permeate./(cumulative_raw_feed[end].+0.23*u"m^3").*100)

    # === Create a Grid of Subplots ===
    plt = plot(layout = (3, 3), size = (1200, 900), right_margin = 1Plots.mm, left_margin = 10Plots.mm)

    # Subplot 1: Cumulative Permeate vs. Time
    plot!(plt[1], time_data, cum_perm,
        xlabel = "Time (min)", ylabel = "Cumulative Permeate (m³)",
        title = "Cumulative Permeate", legend = false)

    # Subplot 2: Cumulative Raw Feed vs. Time
    plot!(plt[2], time_data, cum_raw,
        xlabel = "Time (min)", ylabel = "Cumulative Raw Feed (m³)",
        title = "Cumulative Raw Feed", legend = false)

    # Subplot 3: Flow Rates (Brine Q and Final Feed Q) vs. Time
    plot!(plt[3], time_data, brine_Q,
        xlabel = "Time (min)", ylabel = "Flow Rate (m³/hr)",
        title = "Flow Rates", label = "Brine Q", legend=:right)
    plot!(plt[3], time_data, final_feed_Q, label = "Feed Q setpoint")
    plot!(plt[3], time_data, permeate_Q, label = "Permeate Q")
    ylims!(plt[3], -0.5, 1.05*feed_flowrate_setpoint)

    # Subplot 4: Concentrations (Brine C and Final Feed C) vs. Time
    plot!(plt[4], time_data, brine_C.*1000,
        xlabel = "Time (min)", ylabel = "Concentration (ppm)",
        title = "Concentrations", label = "Brine C")
    plot!(plt[4], time_data, final_feed_C.*1000, label = "Final Feed C")
    ylims!(plt[4], -100.0, maximum(brine_C).*1100)
    # hline!(plt[4], [1000.0], label = "Raw Feed C")

    # Subplot 5: Mass Fluxes (Concentration × Flow Rate) vs. Time
    plot!(plt[5], time_data, mass_flux_brine,
        xlabel = "Time (min)", ylabel = "Salt Mass Flux (kg/s)",
        title = "Mass Fluxes", label = "Brine Mass Flux")
    plot!(plt[5], time_data, mass_flux_final, label = "Final Feed Mass Flux")

    # Subplot 6: Recovery Rate vs. Time
    plot!(plt[6], time_data, recovery_rate,
        xlabel = "Time (min)", ylabel = "Recovery Rate",
        title = "Recovery Rate", legend = false)
    hline!(plt[6], [100.0])
    ylims!(plt[6], -5.0, 105.0)

    # Subplot 7: Permeate Concentration vs. Time
    plot!(plt[7], time_data, ustrip.(total_result.Permeate_C),
        xlabel = "Time (min)", ylabel = "Permeate Concentration (kg/m³)",
        title = "Permeate Concentration", legend = false)

    # Subplot 8 and 9 are left empty (or can be used for additional plots)
    plot!(plt[8], time_data, ustrip.(final_feed_P)./1e5,
        xlabel = "Time (min)", ylabel = "Feed Pressure (bar)",
        title = "Feed Prsesure", legend = false)
    plot!(plt[9], time_data, brine_P, title = "Brine pressure (Pa)")

    display(plt)
    savefig(plt, "summary_plot-adjustedODE.pdf")
end


plot(total_result.time, cumulative_permeate)

plot(total_result.time, cumulative_raw_feed)

plot(total_result.time, total_result.Brine_Q)
plot!(total_result.time, total_result.Final_feed_Q)

plot(total_result.time, total_result.Brine_C)
plot!(total_result.time, total_result.Final_feed_C)

plot(total_result.time, total_result.Brine_C.*total_result.Brine_Q)
plot!(total_result.time, total_result.Final_feed_C.*total_result.Final_feed_Q)

plot(total_result.time, total_result.Permeate_Q.*total_result.Permeate_C./(total_result.Final_feed_C.*total_result.Final_feed_Q))

plot(total_result.time, total_result.Permeate_C)


plot(
    total_result.time, Float64.(cumulative_permeate./(cumulative_raw_feed[end].+1.0*u"m^3")).*100, title="Recovery rate (%)", xlabel="Time"
)
hline!([100.0])
savefig("Recovery.png")

plot(
    total_result.time, total_result.Brine_P, title="Brine pressure", xlabel="Time"
)
savefig("Brine pressure.png")

plot(
    total_result.time, total_result.Final_feed_P, title="Feed pressure", xlabel="Time"
)
savefig("Feed pressure.png")

total_result