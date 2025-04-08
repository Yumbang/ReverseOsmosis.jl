# =================================================== Definition of membrane elements / modules / vessels ===================================================
mutable struct MembraneElement
    # Vessel shape information
    const height::Float64
    const width::Float64
    const dx::Float64
    const spacer_resistance::Float64
    # Membrane property information
    R_m::Float64
    k_fp::Float64
    salt_rejection::Float64
end

mutable struct MembraneElement2
    # Vessel shape information
    const height::Float64
    const width::Float64
    const dx::Float64
    const spacer_resistance::Float64
    # Membrane property information
    R_m::Float64 # Pristine membrane resistance
    R_c::Float64 # Cake resistance
    B::Float64   # Salt permeability
    cake_thickness::Float64
    k_fp::Float64
end

mutable struct MembraneModule
    const n_segments::Int64                       # Currently, single-segment module is not supported.
    const length::Float64
    const width::Float64
    elements_array::Union{Array{MembraneElement}, Array{MembraneElement2}}
end

mutable struct PressureVessel
    const n_modules::Int64
    const n_segments::Array{Int64}
    const length::Array{Float64}
    const width::Array{Float64}
    modules_array::Array{MembraneModule}
end

# =================================================== Definition of Base.show for membrane modules / vessels ===================================================

function Base.show(io::IO, membrane::MembraneModule)
    @printf io "Membrane module with %.2f m length, %.2f m rolled width and %d %s" membrane.length membrane.width membrane.n_segments typeof(membrane.elements_array[1])
end

function Base.show(io::IO, vessel::PressureVessel)
    print(io, "\t⇓\n")
    for (i, membrane) in enumerate(vessel.modules_array)
        print(io, "#$(lpad(i, 2))\t")
        show(io, membrane)
        print(io, "\n")
    end
    print(io, "\t⇓\n")
end




# =================================================== Constructors or functions to construct membrane modules / vessels ===================================================
"""
Construct pristine membrane with uniform dx length.
"""
function pristine_membrane(
    n_segments::Int64, length::Float64, height::Float64, width::Float64,
    R_m::Float64, k_fp::Float64, spacer_resistance::Float64, salt_rejection::Float64
)
    # Construct membrane specification and elements array
    pristine_element        = MembraneElement(
                                        height, width, length/n_segments, spacer_resistance,
                                        R_m, k_fp, salt_rejection
                                        )
    pristine_elements_array = [deepcopy(pristine_element) for _ in range(1,n_segments)]
    # Compose the variables into MembraneModule
    return MembraneModule(n_segments, length, width, pristine_elements_array)
end

function pristine_membrane2(
    n_segments::Int64, length::Float64, height::Float64, width::Float64,
    R_m::Float64, B::Float64,
    k_fp::Float64, spacer_resistance::Float64
)
    # Pristine membrane with no cake resistance
    R_c_pristine            = 0.0
    cake_thickness_pristine = 0.0

    # Construct membrane specification and elements array
    pristine_element        = MembraneElement2(
                               height, width, length/n_segments, spacer_resistance,
                               R_m, R_c_pristine, B, cake_thickness_pristine, k_fp
                               )
    pristine_elements_array = [deepcopy(pristine_element) for _ in range(1,n_segments)]

    # Compose the variables into MembraneModule
    return MembraneModule(n_segments, length, width, pristine_elements_array)
end


"""
Construct pristine pressure vessel with n_modules uniform membrane modules.
"""
function pristine_vessel(
    n_modules::Int64, n_segments::Int64,
    length::Float64, height::Float64, width::Float64,
    R_m::Float64, k_fp::Float64, spacer_resistance::Float64, salt_rejection::Float64
)
    pristine_membrane       = pristine_membrane(n_segments, length, height, width, R_m, k_fp, spacer_resistance, salt_rejection)

    n_modules               = n_modules
    n_segments_array        = [n_segments for _ in range(1,n_modules)]
    length_array            = [length for _ in range(1,n_modules)]
    width_array             = [width  for _ in range(1,n_modules)]
    pristine_modules_array  = [deepcopy(pristine_membrane) for _ in range(1,n_modules)]
    
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, pristine_modules_array)
end

"""
Construct pressure vessel with 1-D array of membrane modules.
"""
function PressureVessel(membranes::Vector{MembraneModule})
    n_modules               = length(membranes)
    n_segments_array        = [membrane.n_segments for membrane in membranes]
    length_array            = [membrane.length for membrane in membranes]
    width_array             = [membrane.width  for membrane in membranes]
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, membranes)
end

"""
Construct pressure vessel with n_modules and a membrane module.
"""
function PressureVessel(n_modules::Int64, membrane::MembraneModule)
    n_modules               = n_modules
    n_segments_array        = [membrane.n_segments for _ in range(1,n_modules)]
    length_array            = [membrane.length for _ in range(1,n_modules)]
    width_array             = [membrane.width  for _ in range(1,n_modules)]
    modules_array           = [deepcopy(membrane) for _ in range(1,n_modules)]
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, modules_array)
end


# =================================================== Profile functions to profile membrane modules / vessels into Unitful DataFrame ===================================================
function unit_element(element::MembraneElement)
    return (
        height = element.height * u"m", width = element.width * u"m",
        dx = element.dx * u"m", spacer_resistance = element.spacer_resistance,
        R_m = element.R_m * u"Pa*s/m", k_fp = element.k_fp * u"Pa*s/m^2",
        salt_rejection = element.salt_rejection
    )
end

function unit_element(element::MembraneElement2)
    return (
        height = element.height * u"m", width = element.width * u"m",
        dx = element.dx * u"m", spacer_resistance = element.spacer_resistance,
        R_m = element.R_m * u"m^-1", R_c = element.R_c * u"m^-1", k_fp = element.k_fp,
        B = element.B * u"kg/m"
    )
end

function profile_element(element::Union{MembraneElement, MembraneElement2})
    element_tuple = unit_element(element)
    return DataFrame([element_tuple])
end

function profile_membrane(membrane::MembraneModule)
    element_tuples = map(unit_element, membrane.elements_array)
    # Turn the array of NamedTuples into a DataFrame
    return DataFrame(element_tuples)
end

function profile_vessel(vessel::PressureVessel)
    return reduce(vcat, [profile_membrane(membrane) for membrane in vessel.modules_array])
end

# =================================================== Functions to manipulate membrane element / modules / vessels (e.g. fouling) ===================================================
"""
Elemental fouling function for membrane element and Water
"""
function foul!(element::MembraneElement, ΔR_m::Float64)
    element.R_m += ΔR_m
    return nothing
end

function foul!(membrane::MembraneModule, ΔR_ms::Array{Float64})
    @assert length(membrane.elements_array) == length(ΔR_ms) "Size of membrane elements and ΔR_ms do not match."
    for element_idx in range(1, length(membrane.elements_array))
        foul!(membrane.elements_array[element_idx], ΔR_ms[element_idx])
    end
    return nothing
end

function foul!(vessel::PressureVessel, ΔR_ms_array::Vector{Vector{Float64}})
    @assert length(vessel.modules_array) == length(ΔR_ms_array) "Size of membrane modules and ΔR_ms array do not match."
    for module_idx in range(1, length(vessel.modules_array))
        foul!(vessel.modules_array[module_idx], ΔR_ms_array[module_idx])
    end
    return nothing
end

"""
Elemental fouling function for membrane element2 and Water2
"""
function foul!(element::MembraneElement2, ΔCake_thickness::Float64, ΔR_c::Float64)
    element.cake_thickness += ΔCake_thickness
    element.R_c += ΔR_c
    return nothing
end

function foul!(membrane::MembraneModule, ΔCake_thicknesses::Array{Float64}, ΔR_cs::Array{Float64})
    @assert length(membrane.elements_array) == length(ΔCake_thicknesses) "Size of membrane elements and ΔCake_thicknesses do not match."
    @assert length(membrane.elements_array) == length(ΔR_cs) "Size of membrane elements and ΔR_cs do not match."
    for element_idx in range(1, length(membrane.elements_array))
        foul!(membrane.elements_array[element_idx], ΔCake_thicknesses[element_idx], ΔR_cs[element_idx])
    end
    return nothing
end

function foul!(vessel::PressureVessel, ΔCake_thicknesses_array::Vector{Vector{Float64}}, ΔR_cs_array::Vector{Vector{Float64}})
    @assert length(vessel.modules_array) == length(ΔCake_thicknesses_array) "Size of membrane modules and ΔCake_thicknesses array do not match."
    @assert length(vessel.modules_array) == length(ΔR_cs_array) "Size of membrane modules and ΔR_cs array do not match."
    for module_idx in range(1, length(vessel.modules_array))
        foul!(vessel.modules_array[module_idx], ΔCake_thicknesses_array[module_idx], ΔR_cs_array[module_idx])
    end
    return nothing
end