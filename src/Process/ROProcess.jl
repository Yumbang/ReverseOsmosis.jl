module ReverseOsmosisProcesses

export
    SinglePassRO, process_singlepass_RO!,
    CirculationPipe, CirculationPipe2, empty_pipe, empty_pipe2, fill_pipe!, update_pipe!,
    SemiBatchRO, process_semibatch_RO!

# TODO: Create CLI interface to construct new process.

include("SinglePass.jl")
include("SemiBatch.jl")
# include("SemiBatchODE.jl")

using .SinglePass
using .SemiBatch

end # module ReverseOsmosisProcesses