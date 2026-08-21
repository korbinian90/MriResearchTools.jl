# Provenance records for command line tools.
#
# Every tool writes two files next to its output: what it was told to do
# (`settings_<tool>.txt`) and what has to be cited for the methods it actually
# ran (`citations_<tool>.txt`). Both used to be written by five near-identical
# `saveconfiguration` functions across three repositories, which is how the
# reference text drifted apart and how citations came to be printed for methods
# that were never used. The text lives here once instead.

# Reference text is written indented for readability in the source; strip that
# back out once, here, so every consumer gets left-aligned text.
function _dedent(text)
    lines = split(text, '\n')
    join([first(lines); lstrip.(lines[2:end])], '\n')
end
_dedent_all(d::Dict) = Dict(k => _dedent(v) for (k, v) in d)

"""
    CITATIONS

Reference text for every method in the toolbox, keyed by method. One source of
truth, so a reference cannot say different things in different tools.

See also [`write_provenance`](@ref).
"""
const CITATIONS = _dedent_all(Dict{Symbol,String}(
    :romeo => """Dymerska, B., Eckstein, K., Bachrata, B., Siow, B., Trattnig, S., Shmueli, K., Robinson, S.D., 2020.
                 Phase Unwrapping with a Rapid Opensource Minimum Spanning TreE AlgOrithm (ROMEO).
                 Magnetic Resonance in Medicine.
                 https://doi.org/10.1002/mrm.28563""",
    :aspire => """Eckstein, K., Dymerska, B., Bachrata, B., Bogner, W., Poljanc, K., Trattnig, S., Robinson, S.D., 2018.
                  Computationally Efficient Combination of Multi-channel Phase Data From Multi-echo Acquisitions (ASPIRE).
                  Magnetic Resonance in Medicine 79, 2996-3006.
                  https://doi.org/10.1002/mrm.26963""",
    :clearswi => """Eckstein, K., Bachrata, B., Hangel, G., Widhalm, G., Enzinger, C., Barth, M., Trattnig, S., Robinson, S.D., 2021.
                    Improved susceptibility weighted imaging at ultra-high field using bipolar multi-echo acquisition and optimized image processing: CLEAR-SWI.
                    NeuroImage 237, 118175.
                    https://doi.org/10.1016/j.neuroimage.2021.118175""",
    :homogeneity => """Eckstein, K., Trattnig, S., Robinson, S.D., 2019.
                       A Simple Homogeneity Correction for Neuroimaging at 7T.
                       Proceedings of the 27th Annual Meeting ISMRM. Presented at the ISMRM, Montreal, Quebec, Canada.
                       https://index.mirasmart.com/ISMRM2019/PDFfiles/2716.html""",
    :laplacian => """Schofield, M.A., Zhu, Y., 2003.
                     Fast phase unwrapping algorithm for interferometric applications.
                     Optics Letters 28, 1194-1196.
                     https://doi.org/10.1364/OL.28.001194""",
    :bestpath => """Abdul-Rahman, H.S., Gdeisat, M.A., Burton, D.R., Lalor, M.J., Lilley, F., Moore, C.J., 2007.
                    Fast and robust three-dimensional best path phase unwrapping algorithm.
                    Applied Optics 46, 6623-6635.
                    https://doi.org/10.1364/AO.46.006623""",
    :tgv => """Langkammer, C., Bredies, K., Poser, B.A., Barth, M., Reishofer, G., Fan, A.P., Bilgic, B., Fazekas, F., Mainero, C., Ropele, S., 2015.
               Fast quantitative susceptibility mapping using 3D EPI and total generalized variation.
               NeuroImage 111, 622-630.
               https://doi.org/10.1016/j.neuroimage.2015.02.041""",
    :tgv_original => """Bredies, K., Ropele, S., Poser, B.A., Barth, M., Langkammer, C., 2014.
                        Single-step quantitative susceptibility mapping using total generalized variation and 3D EPI.
                        Proceedings of the 22nd Annual Meeting ISMRM, p. 604.""",
    :rts => """Kames, C., Wiggermann, V., Rauscher, A., 2018.
               Rapid two-step dipole inversion for susceptibility mapping with sparsity priors.
               NeuroImage 167, 276-286.
               https://doi.org/10.1016/j.neuroimage.2017.11.018""",
    :phase_based_masking => """Hagberg, G.E., Eckstein, K., Tuzzi, E., Zhou, J., Robinson, S.D., Scheffler, K., 2022.
                               Phase-based masking for quantitative susceptibility mapping of the human brain at 9.4T.
                               Magnetic Resonance in Medicine.
                               https://doi.org/10.1002/mrm.29368""",
    :qsmxt => """Stewart, A.W., Robinson, S.D., O'Brien, K., Jin, J., Widhalm, G., Hangel, G., Walls, A., Goodwin, J., Eckstein, K., Tourell, M., Morgan, C., Narayanan, A., Barth, M., Bollmann, S., 2022.
                 QSMxT: Robust masking and artifact reduction for quantitative susceptibility mapping.
                 Magnetic Resonance in Medicine.
                 https://doi.org/10.1002/mrm.29048""",
    :julia => """Bezanson, J., Edelman, A., Karpinski, S., Shah, V.B., 2017.
                 Julia: A fresh approach to numerical computing.
                 SIAM Review 59, 65-98.
                 https://doi.org/10.1137/141000671""",
))

"""
    NOTICES

Non-citation notices that a method carries, keyed the same way as
[`CITATIONS`](@ref). Written into the citations file when the method was
actually used, because that file is where someone looks before publishing or
before shipping a product.
"""
const NOTICES = _dedent_all(Dict{Symbol,String}(
    :aspire => """PATENT: MCPC-3D-S / ASPIRE is covered by US10605885B2
                  (https://patents.google.com/patent/US10605885B2/en). Per the upstream ASPIRE
                  repository, no licence is required for scientific use and the method can be
                  applied free of charge, but a licence IS required for commercial use, and the
                  method is not a medical product, so it may not be used for diagnosis in humans.
                  Note that an MIT licence grants copyright permissions only, not patent rights.""",
))

_fmt(v::AbstractArray) = string(collect(v))
_fmt(v) = string(v)

function _describe_input(path)
    isnothing(path) && return nothing
    p = try abspath(String(path)) catch; String(path) end
    isfile(p) || return "$p (not found)"
    dims = try
        string(size(niread(p)))
    catch
        "unreadable"
    end
    return "$p  $dims"
end

"""
    write_provenance(dir, tool; version, args, settings, cite, optional=Symbol[],
                     inputs=(), packages=())

Write `settings_<tool>.txt` and `citations_<tool>.txt` into `dir`, recording what
was run and what has to be cited for it.

* `version` the application version string.
* `args` the raw command line arguments.
* `settings` any key-value collection of resolved settings. Written sorted, and
  array values are written out rather than skipped, so the echo times actually
  used are recoverable from the record.
* `cite` the methods that were actually used, as keys of [`CITATIONS`](@ref).
  Pass only what ran: a citation for a method the user did not use is as wrong as
  a missing one. Any [`NOTICES`](@ref) entry for those methods is included too.
* `optional` further methods to list under "Optional citations".
* `inputs` `name => path` pairs for the input files, recorded with their
  dimensions.
* `packages` modules whose versions did the work, recorded alongside the Julia
  version so a result can be traced to the code that produced it.

# Examples
```julia-repl
julia> write_provenance("out", "romeo"; version="4.7.1", args=ARGS, settings,
                        cite=[:romeo, :aspire], inputs=["phase" => fn_phase],
                        packages=[ROMEO, MriResearchTools])
```
"""
function write_provenance(dir, tool; version, args, settings, cite,
                          optional=Symbol[], inputs=(), packages=())
    dir = abspath(dir)
    mkpath(dir)
    _write_settings(dir, tool; version, args, settings, inputs, packages)
    _write_citations(dir, tool; cite, optional)
    return dir
end

function _write_settings(dir, tool; version, args, settings, inputs, packages)
    open(joinpath(dir, "settings_$tool.txt"), "w") do io
        println(io, "# $tool $version")
        println(io, "# written: ", _timestamp())
        println(io, "# command: ", join(args, " "))
        println(io)

        println(io, "[versions]")
        println(io, "julia: ", VERSION)
        for m in packages
            v = try string(pkgversion(m)) catch; "unknown" end
            println(io, nameof(m), ": ", v)
        end
        println(io)

        if !isempty(inputs)
            println(io, "[inputs]")
            for (name, path) in inputs
                d = _describe_input(path)
                isnothing(d) || println(io, name, ": ", d)
            end
            println(io)
        end

        println(io, "[settings]")
        for key in sort(collect(keys(settings)); by=string)
            key == "header" && continue
            println(io, key, ": ", _fmt(settings[key]))
        end
    end
end

function _timestamp()
    # Deliberately not using Dates: MriResearchTools does not depend on it and a
    # provenance record is not worth a new dependency for.
    t = round(Int, time())
    return string(Libc.strftime("%Y-%m-%dT%H:%M:%S", t), " (local), unix ", t)
end

function _write_citations(dir, tool; cite, optional)
    known = filter(k -> haskey(CITATIONS, k), unique(cite))
    open(joinpath(dir, "citations_$tool.txt"), "w") do io
        println(io, "# Citations for the methods this run actually used.")
        println(io, "# Methods that were available but not used are deliberately absent.")
        println(io)
        for k in known
            println(io, CITATIONS[k])
            println(io)
        end

        notices = [k for k in known if haskey(NOTICES, k)]
        if !isempty(notices)
            println(io, "# Notices for the methods used:")
            println(io)
            for k in notices
                println(io, NOTICES[k])
                println(io)
            end
        end

        opt = filter(k -> haskey(CITATIONS, k) && k ∉ known, unique(optional))
        if !isempty(opt)
            println(io, "# Optional citations:")
            println(io)
            for k in opt
                println(io, CITATIONS[k])
                println(io)
            end
        end
    end
end
