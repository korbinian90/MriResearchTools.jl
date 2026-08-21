@testitem "provenance" begin
# ROMEO is a dependency of the package, not a direct test dependency, so reach it
# through the parent rather than adding one just for this.
const ROMEO = MriResearchTools.ROMEO

settings = Dict{String,Any}(
    "phase" => "data/small/Phase.nii",
    "TEs" => [4.0, 8.0, 12.0],          # an array: must be recorded, not skipped
    "weights" => "romeo3",
    "verbose" => false,
    "header" => "should not be written", # deliberately excluded
)

dir = mktempdir()
write_provenance(dir, "toolname"; version="9.9.9", args=["-p", "p.nii", "-t", "[4,8,12]"],
                 settings, cite=[:romeo], optional=[:julia],
                 inputs=["phase" => "data/small/Phase.nii"],
                 packages=[ROMEO, MriResearchTools])

s = read(joinpath(dir, "settings_toolname.txt"), String)
c = read(joinpath(dir, "citations_toolname.txt"), String)

@test occursin("# toolname 9.9.9", s)
@test occursin("-p p.nii -t [4,8,12]", s)               # the command is recoverable
@test occursin("julia: $VERSION", s)
@test occursin("ROMEO: $(pkgversion(ROMEO))", s)
# Array settings used to be dropped, which lost the echo times - the one setting
# a result can least afford to be missing.
@test occursin("TEs: [4.0, 8.0, 12.0]", s)
@test occursin("(51, 51, 41, 3)", s)                    # inputs recorded with dimensions
@test !occursin("should not be written", s)             # header excluded

# Citations must cover what ran and nothing else: citing a method the user did
# not use is as wrong as omitting one they did.
@test occursin("Phase Unwrapping with a Rapid Opensource", c)
@test !occursin("ASPIRE", c)
@test occursin("# Optional citations:", c)
@test occursin("Julia: A fresh approach", c)
# Reference text is stored indented for readability; the file must not be.
@test !occursin("\n   Magnetic Resonance in Medicine", c)

# A method with a notice carries it into the record, exactly once.
dir2 = mktempdir()
write_provenance(dir2, "t2"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:aspire, :aspire])
c2 = read(joinpath(dir2, "citations_t2.txt"), String)
@test count("Computationally Efficient Combination", c2) == 1
@test occursin("US10605885B2", c2)
@test occursin("may not be used for diagnosis in humans", c2)

# ... and a run without that method says nothing about it.
dir3 = mktempdir()
write_provenance(dir3, "t3"; version="1", args=String[], settings=Dict{String,Any}(), cite=[:romeo])
c3 = read(joinpath(dir3, "citations_t3.txt"), String)
@test !occursin("US10605885B2", c3)

# Every citation key must resolve, and every notice must belong to a citation.
@test all(haskey(CITATIONS, k) for k in keys(NOTICES))
@test all(!isempty(v) for v in values(CITATIONS))
end
