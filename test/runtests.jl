using MriResearchTools
using Test
using TestItemRunner

@testset "MriResearchTools.jl" begin
    include("masking.jl")
    include("niftihandling.jl")
    include("unwrapping.jl")
    include("utility.jl")
    include("intensitycorrection.jl")
    include("methods.jl")
    include("mcpc3ds.jl")
    include("VSMbasedunwarping.jl")
    include("smoothing.jl")
    include("qsm.jl")
end
