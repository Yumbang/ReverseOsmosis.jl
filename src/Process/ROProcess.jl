module ReverseOsmosisProcesses

export
    SinglePassRO, process_singlepass_RO!,
    CirculationPipe, SemiBatchRO, process_semibatch_RO!

# TODO: Create CLI interface to construct new process.

include("SinglePass.jl")
# include("SemiBatchODE.jl")

using .SinglePass
# using .SemiBatch

end # module ReverseOsmosisProcesses