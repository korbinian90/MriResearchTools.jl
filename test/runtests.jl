using MriResearchTools
using Test

@testset "MriResearchTools.jl" begin
    include("unwrapping_test.jl")
    include("utility_test.jl")
    include("intensitycorrection_test.jl")
end
