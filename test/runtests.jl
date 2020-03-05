using MRI
using Test

@testset "MRI.jl" begin
    #include("unwrapping_test.jl")
    include("utility_test.jl")
    include("intensitycorrection_test.jl")
end
