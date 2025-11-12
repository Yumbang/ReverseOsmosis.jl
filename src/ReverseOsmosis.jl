module ReverseOsmosis

using Printf, DataFrames, Unitful

export
    # Fluid
    Water, Water2, pressurize, profile_water, mix,
    # Membrane
    MembraneElement, MembraneElement2, MembraneModule, PressureVessel,
    pristine_membrane, pristine_membrane2, pristine_vessel, pristine_vessel2,
    profile_membrane, profile_vessel, foul!,
    # Filtration
    osmo_p, pump,
    element_filtration, element_filtration2, module_filtration, vessel_filtration,
    # Process models
    # Single-pass
    SinglePassRO, process_singlepass_RO!

# Basic elements
include("fluid.jl")
include("membrane.jl")

# # Register units (lmh, lmhbar)
# include("ro_units.jl")
# Unitful.register(ROUnits)

# Filtration algorithm
include("filtration.jl")
using .Filtration

# Reverse osmosis process implementations

# Single-pass reverse osmosis process
include("Process/ROProcess.jl")
using .ReverseOsmosisProcesses

end # module ReverseOsmosis
