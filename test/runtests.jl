using MriResearchTools
using Test

@testset "MriResearchTools.jl" begin
    include("unwrapping.jl")
    include("utility.jl")
    include("intensitycorrection.jl")
end
