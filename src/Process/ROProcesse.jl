module ReverseOsmosisProcesses

export
    SinglePassRO,
    process_feed!

using DataFrames, Unitful

using ..ReverseOsmosis: Water, pressurize, profile_water, mix

using ..ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, foul!

using ..ReverseOsmosis: pump, vessel_filtration

# TODO: Create CLI interface to construct new process.

include("./SinglePass.jl")
include("./SemiBatch.jl")

end # module ReverseOsmosisProcesses