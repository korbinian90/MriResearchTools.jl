@testitem "provenance" begin
const ROMEO = MriResearchTools.ROMEO
# The registry holds the references for the methods this package and ROMEO
# implement, registered at load...
for k in (:romeo, :bestpath, :julia, :mcpc3ds, :aspire, :homogeneity, :laplacian, :rts, :phase_based_masking, :qsmxt)
    @test haskey(MriResearchTools.CITATIONS, k)
end
# ... and it must not claim methods it does not implement.
@test !haskey(MriResearchTools.CITATIONS, :clearswi)
@test !haskey(MriResearchTools.CITATIONS, :tgv)

@test package_version(ROMEO) == string(pkgversion(ROMEO))
@test package_version(MriResearchTools) == string(pkgversion(MriResearchTools))

settings = Dict{String,Any}(
    "phase" => "p.nii",
    "TEs" => [4.0, 8.0, 12.0],           # an array: must be recorded, not skipped
    "weights" => "romeo3",
    "header" => "should not be written", # deliberately excluded
)

dir = mktempdir()
write_provenance(dir, "toolname"; version="9.9.9", args=["-p", "p.nii", "-t", "[4,8,12]"],
                 settings, cite=[:romeo, :aspire], optional=[:julia],
                 inputs=["phase" => "data/small/Phase.nii"],
                 packages=[ROMEO, MriResearchTools], describe=describe_input)

s = read(joinpath(dir, "settings_toolname.txt"), String)
c = read(joinpath(dir, "citations_toolname.txt"), String)

@test occursin("# toolname 9.9.9", s)
@test occursin("-p p.nii -t [4,8,12]", s)
@test occursin("julia: $VERSION", s)
@test occursin("ROMEO: $(pkgversion(ROMEO))", s)
@test occursin("MriResearchTools: $(pkgversion(MriResearchTools))", s)
# describe_input adds the dimensions.
@test occursin("(51, 51, 41, 3)", s)
# Array settings used to be dropped, which lost the echo times - the one setting
# a result can least afford to be missing.
@test occursin("TEs: [4.0, 8.0, 12.0]", s)
@test !occursin("should not be written", s)

# Each reference is headed by the method that pulled it in, so a reader can tell
# which step of the run it belongs to without recognising the paper.
@test occursin("## ROMEO Unwrapping\nDymerska", c)
@test occursin("## Julia Scientific Programming Language\nBezanson", c)
@test occursin("# Optional citations:", c)
# Reference text is stored indented for readability; the file must not be.
@test !occursin("\n   Magnetic Resonance in Medicine", c)

# ASPIRE carries a patent, and the citations file is where someone looks before
# publishing or before shipping. It must appear when that method ran...
@test occursin("Computationally Efficient Combination", c)
@test count("Computationally Efficient Combination", c) == 1
@test occursin("US10605885B2", c)
@test occursin("may not be used for diagnosis in humans", c)

# ... and not otherwise.
dir2 = mktempdir()
write_provenance(dir2, "t2"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:romeo])
c2 = read(joinpath(dir2, "citations_t2.txt"), String)
@test occursin("Phase Unwrapping with a Rapid Opensource", c2)
@test !occursin("ASPIRE", c2)

# In particular not for MCPC-3D-S, which shares the ASPIRE paper as its reference
# but is a different method and is not patented.
dir3 = mktempdir()
write_provenance(dir3, "t3"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:mcpc3ds])
c3 = read(joinpath(dir3, "citations_t3.txt"), String)
@test occursin("Computationally Efficient Combination", c3)
@test !occursin("US10605885B2", c3)
@test !occursin("PATENT", c3)

# Every notice must belong to a citation.
@test all(haskey(MriResearchTools.CITATIONS, k) for k in keys(MriResearchTools.NOTICES))

@test occursin("(not found)", describe_input("does/not/exist.nii"))

# Registering is idempotent for identical text, and complains when two packages
# disagree about a reference rather than silently keeping one.
n = length(MriResearchTools.CITATIONS)
register_citation!(:romeo, MriResearchTools.CITATIONS[:romeo])
@test length(MriResearchTools.CITATIONS) == n
@test_logs (:warn,) register_citation!(:romeo, "something else entirely")
@test occursin("Rapid Opensource", MriResearchTools.CITATIONS[:romeo]) # first registration kept

# A method whose owning package is not loaded must be reported, not silently
# dropped: a missing citation is the failure this file exists to prevent.
dir4 = mktempdir()
@test_logs (:warn,) write_provenance(dir4, "t4"; version="1", args=String[],
                                     settings=Dict{String,Any}(), cite=[:not_a_real_method])

# A notice registered with a citation travels with it, exactly once.
register_citation!(:test_method, "Some Reference."; notice="A NOTICE.")
dir5 = mktempdir()
write_provenance(dir5, "t5"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:test_method, :test_method])
c5 = read(joinpath(dir5, "citations_t5.txt"), String)
@test count("Some Reference.", c5) == 1
@test occursin("A NOTICE.", c5)
delete!(MriResearchTools.CITATIONS, :test_method); delete!(MriResearchTools.NOTICES, :test_method)

# One heading per method, not per reference: a method with two references (TGV's
# paper and the abstract that first presented it) registers both under one label
# and they appear together, rather than as two apparent processing steps.
register_citation!(:m_a, "Reference A."; label="Shared Method")
register_citation!(:m_b, "Reference B."; label="Shared Method")
register_citation!(:m_c, "Reference C."; label="Other Method")
dir6 = mktempdir()
write_provenance(dir6, "t6"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:m_a, :m_b, :m_c])
c6 = read(joinpath(dir6, "citations_t6.txt"), String)
@test count("## Shared Method", c6) == 1
@test occursin("## Shared Method\nReference A.\n\nReference B.", c6)
@test occursin("## Other Method\nReference C.", c6)
# An unlabelled key still gets a heading, so the file never mixes labelled and
# bare blocks.
register_citation!(:m_d, "Reference D.")
dir7 = mktempdir()
write_provenance(dir7, "t7"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:m_d])
@test occursin("## m_d\nReference D.", read(joinpath(dir7, "citations_t7.txt"), String))
for k in (:m_a, :m_b, :m_c, :m_d)
    delete!(MriResearchTools.CITATIONS, k); delete!(MriResearchTools.LABELS, k)
end
end

@testitem "provenance: TGV registered by its extension" begin
using QuantitativeSusceptibilityMappingTGV
@test haskey(MriResearchTools.CITATIONS, :tgv)
@test haskey(MriResearchTools.CITATIONS, :tgv_original)
@test MriResearchTools.LABELS[:tgv] == MriResearchTools.LABELS[:tgv_original]
@test package_version(QuantitativeSusceptibilityMappingTGV) == string(pkgversion(QuantitativeSusceptibilityMappingTGV))
end
