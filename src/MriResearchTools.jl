module MriResearchTools

using FFTW
using Interpolations
using NIfTI
using CodecZlib: GzipCompressorStream
using ROMEO
using Statistics
using DataStructures
using LocalFilters
using PaddedViews
using OffsetArrays

# Baked in at precompile time; include_dependency so a version bump invalidates
# the cache. Read from Project.toml rather than through pkgversion, which
# returns nothing when the package is precompiled with --strip-metadata, as
# juliac does.
const PKG_VERSION = let toml = joinpath(@__DIR__, "..", "Project.toml")
    include_dependency(toml)
    m = match(r"^version\s*=\s*\"([^\"]+)\""m, read(toml, String))
    m === nothing && error("no version field in $toml")
    VersionNumber(m.captures[1])
end

include("cli.jl")
include("parse.jl")
include("utility.jl")
include("smoothing.jl")
include("intensitycorrection.jl")
include("VSMbasedunwarping.jl")
include("methods.jl")
include("niftihandling.jl")
include("nifti_static.jl")
include("mcpc3ds.jl")
include("romeofunctions.jl")
include("ice2nii.jl")
include("laplacianunwrapping.jl")
include("masking.jl")
include("qsm_common.jl")
include("provenance.jl")
include("citations.jl")

# Placeholders until a backend extension adds the methods. They throw rather than
# warn: a caller that goes on with the `nothing` a warning returns fails later
# and less clearly, and a compiled program can leave out whatever follows a call
# that always throws.
_no_qsm(f) = error("No QSM implementation is loaded for `$f`. Type `using QuantitativeSusceptibilityMappingTGV` or `using QSM` to load the desired implementation.\n If already loaded, check the expected arguments via `?$f`")
qsm(args...; kwargs...) = _no_qsm("qsm")
qsm_average(args...; kwargs...) = _no_qsm("qsm_average")
qsm_B0(args...; kwargs...) = _no_qsm("qsm_B0")
qsm_laplacian_combine(args...; kwargs...) = _no_qsm("qsm_laplacian_combine")
qsm_romeo_B0(args...; kwargs...) = _no_qsm("qsm_romeo_B0")
phase_based_mask(args...; kwargs...) = error("Load ImageFiltering.jl to use this method: `using ImageFiltering`\n If already loaded, check the expected arguments via `?phase_based_mask`")
if !isdefined(Base, :get_extension)
    include("../ext/QSMExt.jl")
    include("../ext/PhaseBasedMaskingExt.jl")
end

export  readphase, readmag, niread, write_emptynii,
        loadnii, loadphase, loadmag, loadheader,
        header,
        savenii,
        estimatenoise,
        robustmask, robustmask!,
        phase_based_mask,
        mask_from_voxelquality,
        brain_mask,
        robustrescale,
        #combine_echoes,
        calculateB0_unwrapped, get_B0_snr,
        romeovoxelquality,
        getHIP,
        laplacianunwrap, laplacianunwrap!, laplacianunwrap_fft,
        getVSM,
        unwarp,
        thresholdforward,
        gaussiansmooth3d!, gaussiansmooth3d,
        gaussiansmooth3d_phase,
        makehomogeneous!, makehomogeneous,
        getsensitivity,
        getscaledimage,
        estimatequantile,
        RSS,
        mcpc3ds, mcpc3ds_meepi,
        unwrap, unwrap!, romeo, romeo!,
        unwrap_individual, unwrap_individual!,
        homodyne, homodyne!,
        to_dim,
        Ice_output_config, read_volume,
        write_provenance, write_citations, register_citation!, register_version!, describe_input, package_version,
        NumART2star, r2s_from_t2s,
        qsm_average, qsm_B0, qsm_laplacian_combine, qsm_romeo_B0, qsm_mask_filled

end # module
