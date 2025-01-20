module ReverseOsmosisProcesses

export
    SinglePassRO, SemiBatchRO,
    process_feed!

using DataFrames, Unitful

using ..ReverseOsmosis: Water, pressurize, profile_water, mix

using ..ReverseOsmosis: MembraneElement, MembraneModule, PressureVessel,
        profile_membrane, pristine_membrane, foul!

using ..ReverseOsmosis: pump, vessel_filtration

# TODO: Create CLI interface to construct new process.


"""
RO process codes redefine process_feed! function to implement their feed process mechanism.
By process_feed! itself, receiving AbstractROProcess, will rise ImplementationError.
"""
function process_feed!()
    error("ImplementationError: process_feed! not yet implemented.")
end

include("./SinglePass.jl")
include("./SemiBatch.jl")

end # module ReverseOsmosisProcesses