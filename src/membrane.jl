# ===================================================================
# Membrane Data Structures
# ===================================================================
# Hierarchical structure: Element → Module → Vessel
# Elements are discretized membrane segments for spatial resolution
# ===================================================================

"""
    MembraneElement

Mutable structure representing a single discretized membrane segment.
Uses salt rejection coefficient model for simple fouling simulations.

# Fields
## Geometry (immutable)
- `height::Float64`: Channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `dx::Float64`: Axial length of this segment [m]
- `spacer_resistance::Float64`: Feed channel spacer resistance coefficient [-]

## Membrane Properties (mutable for fouling)
- `R_m::Float64`: Membrane hydraulic resistance [Pa·s/m]
- `k_fp::Float64`: Fouling potential coefficient [Pa·s/m²]
- `salt_rejection::Float64`: Salt rejection rate [-], typically 0.99+

# Note
This structure uses a simplified salt rejection model. For more detailed
salt transport modeling including concentration polarization effects,
use `MembraneElement2` instead.
"""
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

"""
    MembraneElement2

Advanced membrane element with solution-diffusion model and cake layer fouling.
Tracks both membrane resistance and deposited cake layer properties.

# Fields
## Geometry (immutable)
- `height::Float64`: Channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `dx::Float64`: Axial length of this segment [m]
- `spacer_resistance::Float64`: Feed channel spacer resistance coefficient [-]

## Membrane Properties (mutable for fouling)
- `R_m::Float64`: Intrinsic membrane resistance [m⁻¹]
- `R_c::Float64`: Cake layer resistance [m⁻¹]
- `B::Float64`: Salt permeability coefficient [kg/m²·s]
- `cake_thickness::Float64`: Accumulated cake layer thickness [m]
- `k_fp::Float64`: Fouling propensity coefficient [-]

# Note
This structure implements:
- Solution-diffusion model for salt transport
- Cake-enhanced concentration polarization (CEMT)
- Carman-Kozeny equation for cake resistance
- Temperature-corrected permeabilities
"""
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

"""
    MembraneModule

Collection of membrane elements arranged in series along flow direction.
Represents a single spiral-wound membrane module.

# Fields
- `n_segments::Int64`: Number of discretized elements (≥ 2 recommended)
- `length::Float64`: Total module length [m]
- `width::Float64`: Membrane rolled width [m]
- `elements_array::Union{Array{MembraneElement}, Array{MembraneElement2}}`:
   Array of membrane elements in flow sequence

# Note
Higher n_segments provides better spatial resolution of concentration
and pressure profiles, but increases computational cost.
"""
mutable struct MembraneModule
    const n_segments::Int64
    const length::Float64
    const width::Float64
    elements_array::Union{Array{MembraneElement}, Array{MembraneElement2}}
end

"""
    PressureVessel

Collection of membrane modules arranged in series.
Represents a complete pressure vessel with multiple modules.

# Fields
- `n_modules::Int64`: Number of modules in series
- `n_segments::Array{Int64}`: Number of segments per module
- `length::Array{Float64}`: Length of each module [m]
- `width::Array{Float64}`: Width of each module [m]
- `modules_array::Array{MembraneModule}`: Array of membrane modules

# Note
In industrial RO systems, multiple modules are placed in series within
a single pressure vessel to maximize membrane area and recovery.
"""
mutable struct PressureVessel
    const n_modules::Int64
    const n_segments::Array{Int64}
    const length::Array{Float64}
    const width::Array{Float64}
    modules_array::Array{MembraneModule}
end

# ===================================================================
# Display Methods
# ===================================================================

"""
    Base.show(io::IO, membrane::MembraneModule)

Custom display format for membrane modules showing key specifications.
"""
function Base.show(io::IO, membrane::MembraneModule)
    @printf io "Membrane module with %.2f m length, %.2f m rolled width and %d %s" membrane.length membrane.width membrane.n_segments typeof(membrane.elements_array[1])
end

"""
    Base.show(io::IO, vessel::PressureVessel)

Custom display format for pressure vessels showing module arrangement.
"""
function Base.show(io::IO, vessel::PressureVessel)
    print(io, "\t⇓\n")
    for (i, membrane) in enumerate(vessel.modules_array)
        print(io, "#$(lpad(i, 2))\t")
        show(io, membrane)
        print(io, "\n")
    end
    print(io, "\t⇓\n")
end

# ===================================================================
# Constructor Functions
# ===================================================================

"""
    pristine_membrane(n_segments, length, height, width, R_m, k_fp, spacer_resistance, salt_rejection)

Construct a pristine (unfouled) membrane module with uniform element discretization.
Uses MembraneElement structure with salt rejection coefficient model.

# Arguments
- `n_segments::Int64`: Number of axial segments (higher = better resolution)
- `length::Float64`: Total module length [m]
- `height::Float64`: Feed channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `R_m::Float64`: Initial membrane resistance [Pa·s/m]
- `k_fp::Float64`: Fouling potential coefficient [Pa·s/m²]
- `spacer_resistance::Float64`: Feed spacer resistance coefficient [-]
- `salt_rejection::Float64`: Salt rejection rate [-], typically 0.990-0.999

# Returns
- `MembraneModule`: Pristine membrane module with uniform properties

# Examples
```julia
# Typical seawater RO membrane (Dow Filmtec SW30HR-380)
membrane = pristine_membrane(
    50,         # 50 segments
    1.016,      # 1.016 m length (40 inches)
    7e-4,       # 0.7 mm channel height
    37.0,       # 37 m² membrane area
    6.8e10,     # Membrane resistance
    0.67e9,     # Fouling potential
    16.0,       # Spacer resistance
    0.995       # 99.5% salt rejection
)
```
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

"""
    pristine_membrane2(n_segments, length, height, width, R_m, B, k_fp, spacer_resistance)

Construct a pristine membrane module using MembraneElement2 structure.
Implements solution-diffusion model with cake layer fouling capability.

# Arguments
- `n_segments::Int64`: Number of axial segments
- `length::Float64`: Total module length [m]
- `height::Float64`: Feed channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `R_m::Float64`: Intrinsic membrane resistance [m⁻¹]
- `B::Float64`: Salt permeability coefficient [kg/m²·s]
- `k_fp::Float64`: Fouling propensity coefficient [-]
- `spacer_resistance::Float64`: Feed spacer resistance coefficient [-]

# Returns
- `MembraneModule`: Pristine membrane module (R_c = 0, cake_thickness = 0)

# Note
This constructor creates a pristine membrane with no initial cake layer.
Cake resistance and thickness accumulate during simulation via fouling.

# Examples
```julia
membrane = pristine_membrane2(
    50,      # 50 segments
    1.016,   # 1.016 m length
    7e-4,    # 0.7 mm channel height
    37.0,    # 37 m² membrane area
    1.5e15,  # Membrane resistance [m⁻¹]
    3.5e-8,  # Salt permeability [kg/m²·s]
    0.8,     # Fouling propensity
    16.0     # Spacer resistance
)
```
"""
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
    pristine_vessel(n_modules, n_segments, length, height, width, R_m, k_fp, spacer_resistance, salt_rejection)

Construct a pristine pressure vessel with uniform membrane modules.
Uses MembraneElement structure with salt rejection model.

# Arguments
- `n_modules::Int64`: Number of modules in series
- `n_segments::Int64`: Number of segments per module
- `length::Float64`: Length per module [m]
- `height::Float64`: Feed channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `R_m::Float64`: Initial membrane resistance [Pa·s/m]
- `k_fp::Float64`: Fouling potential coefficient [Pa·s/m²]
- `spacer_resistance::Float64`: Feed spacer resistance coefficient [-]
- `salt_rejection::Float64`: Salt rejection rate [-]

# Returns
- `PressureVessel`: Pristine vessel with n_modules identical modules

# Examples
```julia
# 7-element pressure vessel (typical industrial configuration)
vessel = pristine_vessel(
    7,       # 7 modules
    50,      # 50 segments per module
    1.016,   # Module length
    7e-4,    # Channel height
    37.0,    # Membrane area
    6.8e10, 0.67e9, 16.0, 0.995
)
```
"""
function pristine_vessel(
    n_modules::Int64, n_segments::Int64,
    length::Float64, height::Float64, width::Float64,
    R_m::Float64, k_fp::Float64, spacer_resistance::Float64, salt_rejection::Float64
)
    pristine_module         = pristine_membrane(n_segments, length, height, width, R_m, k_fp, spacer_resistance, salt_rejection)

    n_modules               = n_modules
    n_segments_array        = [n_segments for _ in range(1,n_modules)]
    length_array            = [length for _ in range(1,n_modules)]
    width_array             = [width  for _ in range(1,n_modules)]
    pristine_modules_array  = [deepcopy(pristine_module) for _ in range(1,n_modules)]
    
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, pristine_modules_array)
end

"""
    pristine_vessel2(n_modules, n_segments, length, height, width, R_m, B, k_fp, spacer_resistance)

Construct a pristine pressure vessel using MembraneElement2 structure.
Implements solution-diffusion model with cake layer fouling capability.

# Arguments
- `n_modules::Int64`: Number of modules in series
- `n_segments::Int64`: Number of segments per module
- `length::Float64`: Length per module [m]
- `height::Float64`: Feed channel height [m]
- `width::Float64`: Membrane rolled width [m]
- `R_m::Float64`: Intrinsic membrane resistance [m⁻¹]
- `B::Float64`: Salt permeability coefficient [kg/m²·s]
- `k_fp::Float64`: Fouling propensity coefficient [-]
- `spacer_resistance::Float64`: Feed spacer resistance coefficient [-]

# Returns
- `PressureVessel`: Pristine vessel with n_modules identical modules
"""
function pristine_vessel2(
    n_modules::Int64, n_segments::Int64,
    length::Float64, height::Float64, width::Float64,
    R_m::Float64, B::Float64,
    k_fp::Float64, spacer_resistance::Float64
)
    pristine_membrane       = pristine_membrane2(n_segments, length, height, width, R_m, B, k_fp, spacer_resistance)

    n_modules               = n_modules
    n_segments_array        = [n_segments for _ in range(1,n_modules)]
    length_array            = [length for _ in range(1,n_modules)]
    width_array             = [width  for _ in range(1,n_modules)]
    pristine_modules_array  = [deepcopy(pristine_membrane) for _ in range(1,n_modules)]

    return PressureVessel(n_modules, n_segments_array, length_array, width_array, pristine_modules_array)
end

"""
    PressureVessel(membranes::Vector{MembraneModule})

Construct a pressure vessel from an array of membrane modules.
Allows heterogeneous module configurations within a single vessel.

# Arguments
- `membranes::Vector{MembraneModule}`: Array of membrane modules in flow sequence

# Returns
- `PressureVessel`: Vessel containing the provided modules

# Note
This constructor is useful for modeling vessels with mixed membrane types
or degraded/fouled modules at different positions.
"""
function PressureVessel(membranes::Vector{MembraneModule})
    n_modules               = length(membranes)
    n_segments_array        = [membrane.n_segments for membrane in membranes]
    length_array            = [membrane.length for membrane in membranes]
    width_array             = [membrane.width  for membrane in membranes]
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, membranes)
end

"""
    PressureVessel(n_modules::Int64, membrane::MembraneModule)

Construct a pressure vessel with n_modules copies of a single membrane module.
Convenience constructor for uniform vessel configurations.

# Arguments
- `n_modules::Int64`: Number of modules in series
- `membrane::MembraneModule`: Template membrane module to replicate

# Returns
- `PressureVessel`: Vessel with n_modules deep copies of the template

# Examples
```julia
module = pristine_membrane(50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
vessel = PressureVessel(7, module)  # 7-element vessel
```
"""
function PressureVessel(n_modules::Int64, membrane::MembraneModule)
    n_modules               = n_modules
    n_segments_array        = [membrane.n_segments for _ in range(1,n_modules)]
    length_array            = [membrane.length for _ in range(1,n_modules)]
    width_array             = [membrane.width  for _ in range(1,n_modules)]
    modules_array           = [deepcopy(membrane) for _ in range(1,n_modules)]
    return PressureVessel(n_modules, n_segments_array, length_array, width_array, modules_array)
end

# ===================================================================
# Profiling Functions for Data Export
# ===================================================================

"""
    unit_element(element::MembraneElement)

Convert membrane element properties to NamedTuple with Unitful units.
Internal helper function for DataFrame creation.
"""
function unit_element(element::MembraneElement)
    return (
        height = element.height * u"m", width = element.width * u"m",
        dx = element.dx * u"m", spacer_resistance = element.spacer_resistance,
        R_m = element.R_m * u"Pa*s/m", k_fp = element.k_fp * u"Pa*s/m^2",
        salt_rejection = element.salt_rejection
    )
end

"""
    unit_element(element::MembraneElement2)

Convert MembraneElement2 properties to NamedTuple with Unitful units.
"""
function unit_element(element::MembraneElement2)
    return (
        height = element.height * u"m", width = element.width * u"m",
        dx = element.dx * u"m", spacer_resistance = element.spacer_resistance, cake_thickness = element.cake_thickness * u"m",
        R_m = element.R_m * u"m^-1", R_c = element.R_c * u"m^-1", k_fp = element.k_fp,
        B = element.B * u"kg/m"
    )
end

"""
    profile_element(element::Union{MembraneElement, MembraneElement2})

Convert single membrane element to Unitful DataFrame.
"""
function profile_element(element::Union{MembraneElement, MembraneElement2})
    element_tuple = unit_element(element)
    return DataFrame([element_tuple])
end

"""
    profile_membrane(membrane::MembraneModule)

Create DataFrame with spatial profile of membrane module properties.
Each row represents one membrane element in the flow direction.

# Arguments
- `membrane::MembraneModule`: Membrane module to profile

# Returns
- `DataFrame`: DataFrame with element-by-element membrane properties

# Note
Useful for visualizing spatial distribution of fouling, concentration, etc.
"""
function profile_membrane(membrane::MembraneModule)
    element_tuples = map(unit_element, membrane.elements_array)
    # Turn the array of NamedTuples into a DataFrame
    return DataFrame(element_tuples)
end

"""
    profile_vessel(vessel::PressureVessel)

Create DataFrame with complete spatial profile of pressure vessel.
Concatenates profiles from all modules in the vessel.

# Arguments
- `vessel::PressureVessel`: Pressure vessel to profile

# Returns
- `DataFrame`: DataFrame with all elements from all modules

# Examples
```julia
vessel = pristine_vessel(7, 50, 1.016, 7e-4, 37.0, 6.8e10, 0.67e9, 16.0, 0.995)
profile = profile_vessel(vessel)  # 350 rows (7 modules × 50 segments)
```
"""
function profile_vessel(vessel::PressureVessel)
    return reduce(vcat, [profile_membrane(membrane) for membrane in vessel.modules_array])
end

# ===================================================================
# Fouling Functions (Mutate Membrane Properties)
# ===================================================================

"""
    foul!(element::MembraneElement, ΔR_m::Float64)

Apply fouling to a single membrane element by increasing resistance.
Mutates the element in-place.

# Arguments
- `element::MembraneElement`: Element to foul
- `ΔR_m::Float64`: Resistance increase [Pa·s/m]

# Note
This function is typically called automatically during simulation.
ΔR_m can be a rate (when dt=nothing) or absolute change (when dt is specified).
"""
function foul!(element::MembraneElement, ΔR_m::Float64)
    element.R_m += ΔR_m
    return nothing
end

"""
    foul!(membrane::MembraneModule, ΔR_ms::Array{Float64})

Apply spatially-distributed fouling to all elements in a membrane module.

# Arguments
- `membrane::MembraneModule`: Module to foul
- `ΔR_ms::Array{Float64}`: Array of resistance increases for each element [Pa·s/m]

# Note
Array length must match number of membrane elements.
"""
function foul!(membrane::MembraneModule, ΔR_ms::Array{Float64})
    @assert length(membrane.elements_array) == length(ΔR_ms) "Size of membrane elements and ΔR_ms do not match."
    for element_idx in range(1, length(membrane.elements_array))
        foul!(membrane.elements_array[element_idx], ΔR_ms[element_idx])
    end
    return nothing
end

"""
    foul!(vessel::PressureVessel, ΔR_ms_array::Vector{Vector{Float64}})

Apply spatially-distributed fouling to all modules in a pressure vessel.

# Arguments
- `vessel::PressureVessel`: Vessel to foul
- `ΔR_ms_array::Vector{Vector{Float64}}`: Nested array of resistance increases
                                           Outer array: modules, Inner array: elements

# Note
This function is automatically called during simulation with fouling=true.
"""
function foul!(vessel::PressureVessel, ΔR_ms_array::Vector{Vector{Float64}})
    @assert length(vessel.modules_array) == length(ΔR_ms_array) "Size of membrane modules and ΔR_ms array do not match."
    for module_idx in range(1, length(vessel.modules_array))
        foul!(vessel.modules_array[module_idx], ΔR_ms_array[module_idx])
    end
    return nothing
end

"""
    foul!(element::MembraneElement2, ΔCake_thickness::Float64, ΔR_c::Float64)

Apply cake layer fouling to MembraneElement2.
Updates both cake thickness and cake resistance.

# Arguments
- `element::MembraneElement2`: Element to foul
- `ΔCake_thickness::Float64`: Cake thickness increase [m]
- `ΔR_c::Float64`: Cake resistance increase [m⁻¹]

# Note
For MembraneElement2, fouling includes both cake deposition and
resistance change calculated via Carman-Kozeny equation.
"""
function foul!(element::MembraneElement2, ΔCake_thickness::Float64, ΔR_c::Float64)
    element.cake_thickness += ΔCake_thickness
    element.R_c += ΔR_c
    return nothing
end

"""
    foul!(membrane::MembraneModule, ΔCake_thicknesses::Array{Float64}, ΔR_cs::Array{Float64})

Apply spatially-distributed cake layer fouling to membrane module.
"""
function foul!(membrane::MembraneModule, ΔCake_thicknesses::Array{Float64}, ΔR_cs::Array{Float64})
    @assert length(membrane.elements_array) == length(ΔCake_thicknesses) "Size of membrane elements and ΔCake_thicknesses do not match."
    @assert length(membrane.elements_array) == length(ΔR_cs) "Size of membrane elements and ΔR_cs do not match."
    for element_idx in range(1, length(membrane.elements_array))
        foul!(membrane.elements_array[element_idx], ΔCake_thicknesses[element_idx], ΔR_cs[element_idx])
    end
    return nothing
end

"""
    foul!(vessel::PressureVessel, ΔCake_thicknesses_array, ΔR_cs_array)

Apply spatially-distributed cake layer fouling to pressure vessel.
For use with MembraneElement2 and Water2 simulations.
"""
function foul!(vessel::PressureVessel, ΔCake_thicknesses_array::Vector{Vector{Float64}}, ΔR_cs_array::Vector{Vector{Float64}})
    @assert length(vessel.modules_array) == length(ΔCake_thicknesses_array) "Size of membrane modules and ΔCake_thicknesses array do not match."
    @assert length(vessel.modules_array) == length(ΔR_cs_array) "Size of membrane modules and ΔR_cs array do not match."
    for module_idx in range(1, length(vessel.modules_array))
        foul!(vessel.modules_array[module_idx], ΔCake_thicknesses_array[module_idx], ΔR_cs_array[module_idx])
    end
    return nothing
end