@testitem "provenance" begin
# The writer lives in ROMEO; what this package owns is the references for the
# methods it implements, registered at load.
const ROMEO = MriResearchTools.ROMEO

for k in (:aspire, :homogeneity, :laplacian, :rts, :phase_based_masking, :qsmxt)
    @test haskey(ROMEO.CITATIONS, k)
end
# ... and it must not claim methods it does not implement.
@test !haskey(ROMEO.CITATIONS, :clearswi)

# Loading this package makes ROMEO's registry richer; it is one shared registry,
# not a copy per package.
@test ROMEO.CITATIONS === MriResearchTools.ROMEO.CITATIONS

settings = Dict{String,Any}("TEs" => [4.0, 8.0, 12.0], "weights" => "romeo3")

dir = mktempdir()
write_provenance(dir, "toolname"; version="9.9.9", args=["-p", "p.nii"],
                 settings, cite=[:romeo, :aspire], optional=[:julia],
                 inputs=["phase" => "data/small/Phase.nii"],
                 packages=[ROMEO, MriResearchTools], describe=describe_input)

s = read(joinpath(dir, "settings_toolname.txt"), String)
c = read(joinpath(dir, "citations_toolname.txt"), String)

@test occursin("MriResearchTools: $(pkgversion(MriResearchTools))", s)
# describe_input adds the dimensions, which is why it lives here rather than in
# the dependency-free writer.
@test occursin("(51, 51, 41, 3)", s)
@test occursin("TEs: [4.0, 8.0, 12.0]", s)

# MCPC-3D-S carries a patent, and the citations file is where someone looks
# before publishing or before shipping. It must appear when the method ran...
@test occursin("Computationally Efficient Combination", c)
@test count("Computationally Efficient Combination", c) == 1
@test occursin("US10605885B2", c)
@test occursin("may not be used for diagnosis in humans", c)

# ... and not otherwise.
dir2 = mktempdir()
write_provenance(dir2, "t2"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:romeo])
@test !occursin("US10605885B2", read(joinpath(dir2, "citations_t2.txt"), String))

# Every notice must belong to a citation.
@test all(haskey(ROMEO.CITATIONS, k) for k in keys(ROMEO.NOTICES))

@test occursin("(not found)", describe_input("does/not/exist.nii"))
end
