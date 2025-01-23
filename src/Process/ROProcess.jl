module ReverseOsmosisProcesses

export
    SinglePassRO, process_singlepass_RO!,
    SemiBatchRO, process_semibatch_RO!

# TODO: Create CLI interface to construct new process.

include("SinglePass.jl")
include("SemiBatch.jl")
# include("SemiBatchODE.jl")

using .SinglePass
using .SemiBatch
# using .SemiBatchODE

end # module ReverseOsmosisProcesses