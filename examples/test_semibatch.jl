using Pkg
Pkg.activate("/Users/ybang_mac/Documents/병찬/Research/2025/Batch RO/Code/ReverseOsmosis.jl")
using Revise
using ReverseOsmosis
using DataFrames
using BenchmarkTools
using Unitful
using Plots

pipe = CirculationPipe(5.0)
membrane_spec = (100, 1.016, 7e-4, 37/1.016, 6.798858730956808e10, 0.6741807585817046e9, 16.0, 0.995)
membrane = pristine_membrane(membrane_spec...);
vessel = PressureVessel(7, membrane)
sbro = SemiBatchRO(0.8, 0.8, pipe, deepcopy(vessel))

feed_T = 20.0
feed_C = 1.0
feed_P = 1e-5
u₀     = nothing

flowrate_setpoint = 20.0
pressure_setpoint = 10.0

dt      = 10.0u"s"
mode    = :CC
fouling = true