"""
    ReverseOsmosisProcesses

Module containing high-level reverse osmosis process models.
Currently implements single-pass RO systems.

# Exported Types and Functions
- `SinglePassRO`: Single-pass RO process structure
- `process_singlepass_RO!`: Process feed water through single-pass RO system
"""
module ReverseOsmosisProcesses

export
    SinglePassRO, process_singlepass_RO!

include("SinglePass.jl")

using .SinglePass

end # module ReverseOsmosisProcesses