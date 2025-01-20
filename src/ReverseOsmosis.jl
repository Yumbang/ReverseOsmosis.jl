module ReverseOsmosis

using Printf, DataFrames, Unitful

export
    # Fluid
    Water, pressurize, profile_water, mix,
    # Membrane
    MembraneElement, MembraneModule, PressureVessel,
    pristine_membrane, pristine_vessel,
    profile_membrane, profile_vessel,
    # Filtration
    osmo_p, pump,
    element_filtration, module_filtration, vessel_filtration,
    # Process models
    SinglePassRO, process_singlepass_RO!,
    SemiBatchRO, process_semibatch_RO!

# Basic elements
include("fluid.jl")
include("membrane.jl")

# Filtration algorithm
include("filtration.jl")
using .Filtration

# Reverse osmosis process implementations

# Single-pass reverse osmosis process
include("Process/ROProcess.jl")
using .ReverseOsmosisProcesses

end # module ReverseOsmosis
