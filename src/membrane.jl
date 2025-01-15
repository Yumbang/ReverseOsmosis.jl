module Membrane
using DataFrames

export MembraneElement, MembraneModule, pristine_membrane, profile_membrane

mutable struct MembraneElement
    # Vessel shape information
    length::Float64
    height::Float64
    width::Float64
    dx::Float64
    # Membrane property information
    R_m::Float64
    k_fp::Float64
end


mutable struct MembraneModule
    n_segments::Int64                       # Currently, single-segment module is not supported.
    elements_array::Array{MembraneElement}
end

function Base.show(io::IO, membrane::MembraneModule)
    print(io, "Membrane module with $(membrane.n_segments) elements")
end

function pristine_membrane(n_segments::Int64, length::Float64, height::Float64, width::Float64, R_m::Float64, k_fp::Float64)
    # Construct membrane specification and elements array
    pristine_element        = MembraneElement(
                                        length, height, width, length/n_segments,
                                        R_m, k_fp
                                        )
    pristine_elements_array = [deepcopy(pristine_element) for _ in range(1,n_segments)]
    # Compose the variables into MembraneModule
    return MembraneModule(n_segments, pristine_elements_array)
end

function profile_membrane(membrane::MembraneModule)
    # Gather field names from the struct type
    names = fieldnames(MembraneElement)
    # Convert each element of arr to a NamedTuple
    named_tuples = [
        NamedTuple{names}(getfield(s, f) for f in names) for s in membrane.elements_array
    ]
    # Turn the array of NamedTuples into a DataFrame
    return DataFrame(named_tuples)
end

end