@testitem "parse_array" begin
const parse_array = MriResearchTools.parse_array
const parse_weight_flags = MriResearchTools.parse_weight_flags


# Every form the command line accepts. `nargs = '+'` is why the same value
# arrives split across several strings when separated by spaces.
cases = [
    (["[1.5,3.0]"],     [1.5, 3.0]),
    (["[1.5,", "3.0]"], [1.5, 3.0]),
    (["1.5", "3.0"],    [1.5, 3.0]),
    (["[1.5", "3.0]"],  [1.5, 3.0]),
    (["1:3"],           [1, 2, 3]),
    (["3.5:3.5:14"],    [3.5, 7.0, 10.5, 14.0]),
    ([":"],             Colon()),
    (["3"],             3),
    (["[4,4,0]"],       [4, 4, 0]),
    (["1", "2", "3"],   [1, 2, 3]),
    (["-1.5,2"],        [-1.5, 2.0]),
    (["1e-3", "2e-3"],  [0.001, 0.002]),
]
for (input, want) in cases
    @test parse_array(input) == want
    @test typeof(parse_array(input)) == typeof(want)
end

# All-integer input stays integer, because echo selections index with it.
@test parse_array(["1:3"]) isa Vector{Int}
@test parse_array(["1,2.5"]) isa Vector{Float64}
@test parse_array(["3"]) isa Int

# What the user types is data, not code: anything that is not a number, range or
# colon is refused, with the offending text named.
for bad in (["run(`whoami`)"], ["exit(1)"], ["[1,foo]"], ["1:2:3:4"], [""])
    @test_throws ArgumentError parse_array(bad)
end

# Ranges and the colon are only ever used to index with.
@test (1:3)[parse_array([":"])] == [1, 2, 3]
@test (1:8)[parse_array(["1:3"])] == [1, 2, 3]

@test parse_weight_flags("1010") == Bool[1, 0, 1, 0, 0, 0]
@test parse_weight_flags("111111") == trues(6)
@test parse_weight_flags("romeo3") === nothing  # a named weighting, not flags
@test parse_weight_flags("") === nothing
@test parse_weight_flags("0101010") === nothing # seven bits: not flags either

end
