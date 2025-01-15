module ReverseOsmosis

export
    # Fluid
    Water, pressurize, profile_water, mix,
    # Membrane
    MembraneModule, MembraneElement, pristine_membrane,
    # Filtration
    osmo_p, element_filtration, pump

include("filtration.jl")

using .Filtration

end # module ReverseOsmosis
