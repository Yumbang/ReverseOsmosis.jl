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
    SinglePassRO, process_feed!

include("fluid.jl")
include("membrane.jl")

include("filtration.jl")
using .Filtration


include("Process/SinglePass.jl")
using .SinglePassReverseOsmosis

end # module ReverseOsmosis
