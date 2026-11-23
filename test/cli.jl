@testitem "CLI" begin
CLI = MriResearchTools.CLI
spec = CLI.Spec("tool", "1.2.3", [
    CLI.Option("--phase", "-p", "the phase"),
    CLI.Option("--echo-times", "-t", "echo times"; nargs=:many),
    CLI.Option("--compute-B0", "-B", "B0"; nargs=:optional),
    CLI.Option("--verbose", "-v", "verbose"; nargs=:none),
    CLI.Option("--quality", "-q", "quality"; nargs=:none),
    CLI.Option("--threshold", "", "threshold"),
])

v = CLI.parse(spec, ["-p", "a.nii", "-t", "[2", "4", "6]", "-vq", "-B", "--threshold=4"])
@test v["phase"] == ["a.nii"]
@test v["echo-times"] == ["[2", "4", "6]"]
@test v["compute-B0"] == String[] # given without a value
@test haskey(v, "verbose") && haskey(v, "quality")
@test v["threshold"] == ["4"]
@test !haskey(v, "magnitude")
@test CLI.parse(spec, ["-B", "name", "-pa.nii"]) == Dict("compute-B0" => ["name"], "phase" => ["a.nii"])
@test CLI.parse(spec, ["--threshold", "-1"])["threshold"] == ["-1"] # a negative number is a value
@test CLI.parse(spec, ["-B", "-t", "2:2:6"]) == Dict("compute-B0" => String[], "echo-times" => ["2:2:6"])
# a unique prefix of a long name is accepted, as ArgParse does
@test CLI.parse(spec, ["--thresh", "4"]) == Dict("threshold" => ["4"])
@test CLI.parse(spec, ["--verb"]) == Dict("verbose" => String[])
two = CLI.Spec("tool", "1.2.3", [CLI.Option("--mask", "", "mask"), CLI.Option("--mag", "", "mag")])
@test CLI.parse(two, ["--mag", "a.nii"]) == Dict("mag" => ["a.nii"]) # an exact name wins
@test_throws ArgumentError CLI.parse(two, ["--ma", "a.nii"])
@test_throws ArgumentError CLI.parse(spec, ["--ver"]) # --verbose or --version
@test redirect_stdout(() -> CLI.parse(spec, ["--vers"]), devnull) === nothing
@test_throws ArgumentError CLI.parse(spec, ["--version=1"])
@test_throws ArgumentError CLI.parse(spec, ["--he="])
for bad in (["-x"], ["--threshold"], ["-t"], ["extra"], ["--verbose=1"], ["--bogus=1"])
    @test_throws ArgumentError CLI.parse(spec, bad)
end
@test redirect_stdout(() -> CLI.parse(spec, ["--help"]), devnull) === nothing
@test redirect_stdout(() -> CLI.parse(spec, ["--version"]), devnull) === nothing

h = CLI.help(spec)
@test occursin("usage: tool [-p PHASE] [-t ECHO-TIMES...] [-B [COMPUTE-B0]] [-v] [-q]", h)
@test occursin("  -t, --echo-times ECHO-TIMES...\n", h)
@test occursin("  --threshold THRESHOLD threshold\n", h)
@test all(length(line) <= 80 for line in split(h, '\n'))

@test CLI.format([2, 4, 6]) == "[2, 4, 6]"
@test CLI.format([2.5]) == "[2.5]"
@test CLI.format(["a", "b"]) == "a b"
@test CLI.format(String[]) == "(not set)"
@test CLI.format("") == "(not set)"
@test CLI.format(Inf) == "Inf"
@test CLI.format(true) == "true"
@test CLI.format(nothing) == "(not set)"
@test !CLI.static_binary() # this is a Julia session, not a compiled program

end
