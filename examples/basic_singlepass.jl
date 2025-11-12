# Basic Single-Pass RO Example
# Simulates a seawater desalination system with fouling over time

using ReverseOsmosis
using Unitful
using DataFrames
using Printf

println("="^60)
println("Single-Pass Reverse Osmosis Simulation")
println("="^60)

# =============================================================================
# 1. Create Membrane Pressure Vessel
# =============================================================================
println("\n[1/4] Creating pressure vessel...")

# Typical industrial configuration: 7 modules in series
# Module specs based on Dow Filmtec SW30HR-380 (seawater RO)
vessel = pristine_vessel(
    7,          # 7 modules in series (typical)
    50,         # 50 axial segments per module (for spatial resolution)
    1.016,      # 1.016 m = 40 inches (standard module length)
    7e-4,       # 0.7 mm channel height
    37.0,       # 37 m² membrane area per module
    6.8e10,     # Membrane resistance [Pa·s/m]
    0.67e9,     # Fouling potential coefficient [Pa·s/m²]
    16.0,       # Spacer resistance coefficient
    0.995       # 99.5% salt rejection
)

println("  ✓ Created 7-module pressure vessel (350 total elements)")

# =============================================================================
# 2. Create Single-Pass RO Process
# =============================================================================
println("\n[2/4] Initializing process...")

# Create process with 80% pump efficiency
process = SinglePassRO(0.8, vessel)

println("  ✓ Process initialized with 80% pump efficiency")

# =============================================================================
# 3. Define Operating Conditions
# =============================================================================
println("\n[3/4] Setting operating conditions...")

# Seawater feed at atmospheric pressure
feed_flowrate = 10.0        # m³/h
feed_temperature = 25.0     # °C
feed_TDS = 35.0             # kg/m³ (typical seawater)
feed_pressure = 1e5         # Pa (atmospheric)

# Operating parameters
operating_pressure = 55e5   # Pa (55 bar, typical for seawater RO)
time_step = 3600.0          # seconds (1 hour)
simulation_hours = 100      # Simulate 100 hours

println("  Feed: $(feed_flowrate) m³/h, $(feed_TDS) kg/m³ TDS, $(feed_temperature)°C")
println("  Operating pressure: $(operating_pressure/1e5) bar")
println("  Simulation duration: $(simulation_hours) hours")

# =============================================================================
# 4. Run Simulation
# =============================================================================
println("\n[4/4] Running simulation...")

# Storage for results
results = DataFrame(
    time_hr = Float64[],
    permeate_flow = Float64[],
    permeate_TDS = Float64[],
    brine_flow = Float64[],
    brine_TDS = Float64[],
    recovery_pct = Float64[],
    power_kW = Float64[]
)

# Time-stepping loop
for hour in 1:simulation_hours
    # Create feed water at current conditions
    feed = Water(feed_flowrate, feed_temperature, feed_TDS, feed_pressure)

    # Simulate one hour
    brine, permeate, _, _, power = process_singlepass_RO!(
        feed,
        operating_pressure * u"Pa";
        process=process,
        dt=time_step * u"s",
        fouling=true,           # Enable fouling
        profile_process=false   # Don't store detailed profiles
    )

    # Calculate metrics
    recovery = (permeate.Q / feed.Q) * 100

    # Store results
    push!(results, (
        hour,
        permeate.Q,
        permeate.C,
        brine.Q,
        brine.C,
        recovery,
        ustrip(u"kW", power)
    ))

    # Progress indicator every 10 hours
    if hour % 10 == 0
        @printf("  Hour %3d: Permeate %.2f m³/h (%.1f%% recovery), Power %.1f kW\n",
                hour, permeate.Q, recovery, ustrip(u"kW", power))
    end
end

println("  ✓ Simulation complete")

# =============================================================================
# 5. Results Summary
# =============================================================================
println("\n" * "="^60)
println("RESULTS SUMMARY")
println("="^60)

# Initial performance
println("\n📊 Initial Performance (Hour 1):")
@printf("  Permeate production: %.2f m³/h\n", results.permeate_flow[1])
@printf("  Permeate TDS: %.2f kg/m³ (%.0f ppm)\n", results.permeate_TDS[1], results.permeate_TDS[1] * 1000)
@printf("  Recovery rate: %.1f%%\n", results.recovery_pct[1])
@printf("  Salt rejection: %.2f%%\n", (1 - results.permeate_TDS[1] / feed_TDS) * 100)
@printf("  Power consumption: %.1f kW\n", results.power_kW[1])
@printf("  Specific energy: %.2f kWh/m³\n", results.power_kW[1] / results.permeate_flow[1])

# Final performance
println("\n📊 Final Performance (Hour $(simulation_hours)):")
@printf("  Permeate production: %.2f m³/h\n", results.permeate_flow[end])
@printf("  Permeate TDS: %.2f kg/m³ (%.0f ppm)\n", results.permeate_TDS[end], results.permeate_TDS[end] * 1000)
@printf("  Recovery rate: %.1f%%\n", results.recovery_pct[end])
@printf("  Salt rejection: %.2f%%\n", (1 - results.permeate_TDS[end] / feed_TDS) * 100)
@printf("  Power consumption: %.1f kW\n", results.power_kW[end])
@printf("  Specific energy: %.2f kWh/m³\n", results.power_kW[end] / results.permeate_flow[end])

# Performance decline due to fouling
flux_decline = (1 - results.permeate_flow[end] / results.permeate_flow[1]) * 100
println("\n📉 Fouling Impact:")
@printf("  Flux decline: %.1f%%\n", flux_decline)
@printf("  Production loss: %.2f m³/h\n", results.permeate_flow[1] - results.permeate_flow[end])

# Total production
total_permeate = sum(results.permeate_flow) # m³
total_energy = sum(results.power_kW)        # kWh
avg_specific_energy = total_energy / total_permeate

println("\n📈 Cumulative Totals ($(simulation_hours) hours):")
@printf("  Total permeate produced: %.1f m³\n", total_permeate)
@printf("  Total energy consumed: %.1f kWh\n", total_energy)
@printf("  Average specific energy: %.2f kWh/m³\n", avg_specific_energy)

# Membrane condition
println("\n🧪 Membrane Condition:")
vessel_profile = profile_vessel(process.pressure_vessel)
initial_R_m = 6.8e10
final_R_m_avg = mean(vessel_profile.R_m)
final_R_m_max = maximum(vessel_profile.R_m)
resistance_increase_avg = ((final_R_m_avg / initial_R_m) - 1) * 100
resistance_increase_max = ((final_R_m_max / initial_R_m) - 1) * 100

@printf("  Average resistance increase: %.1f%%\n", resistance_increase_avg)
@printf("  Maximum resistance increase: %.1f%%\n", resistance_increase_max)
@printf("  (Lead elements foul more than tail elements)\n")

println("\n" * "="^60)
println("✓ Example completed successfully!")
println("="^60)

# Optional: Save results to CSV
# using CSV
# CSV.write("singlepass_results.csv", results)
# println("\n💾 Results saved to: singlepass_results.csv")

# Optional: Visualization
# using Plots
# plot(results.time_hr, results.permeate_flow,
#      xlabel="Time [hours]",
#      ylabel="Permeate Flow [m³/h]",
#      title="Membrane Performance Decline Due to Fouling",
#      label="Permeate Production",
#      linewidth=2)
# savefig("fouling_decline.png")
# println("📊 Plot saved to: fouling_decline.png")
