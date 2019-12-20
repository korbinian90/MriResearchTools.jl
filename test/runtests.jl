using MRI
using Test

@testset "MRI.jl" begin
    include("utility_test.jl")
    include("intensitycorrection_test.jl")
end
