module MriResearchTools

using FFTW
using Interpolations
using NIfTI
using ROMEO
using Statistics
using DataStructures
using ImageMorphology
using LocalFilters
using PaddedViews
using OffsetArrays
import StatsBase: countmap

include("utility.jl")
include("smoothing.jl")
include("intensitycorrection.jl")
include("VSMbasedunwarping.jl")
include("methods.jl")
include("niftihandling.jl")
include("mcpc3ds.jl")
include("romeofunctions.jl")
include("ice2nii.jl")
include("laplacianunwrapping.jl")
include("masking.jl")

qsm(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm`")
qsm_average(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm_average`")
qsm_B0(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm_B0`")
qsm_laplacian_combine(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm_laplacian_combine`")
qsm_romeo_B0(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm_romeo_B0`")
qsm_mask_filled(args...; kwargs...) = @warn("Type `using TGV_QSM` or `using QSM` to load the desired implementation \n If already loadad, check expected arguments via `?qsm_mask_filled`")
phase_based_mask(args...; kwargs...) = @warn("Load ImageFiltering.jl to use this method: `using ImageFiltering`\n If already loadad, check expected arguments via `?phase_based_masking`")
if !isdefined(Base, :get_extension)
    include("../ext/QSMExt.jl")
    include("../ext/PhaseBasedMaskingExt.jl")
end

export  readphase, readmag, niread, write_emptynii,
        header, affine,
        affine_transformation,
        savenii,
        estimatenoise,
        robustmask, robustmask!,
        phase_based_mask,
        mask_from_voxelquality,
        brain_mask,
        robustrescale,
        #combine_echoes,
        calculateB0_unwrapped,
        romeovoxelquality,
        getHIP,
        laplacianunwrap, laplacianunwrap!,
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
        NumART2star, r2s_from_t2s,
        qsm_average, qsm_B0, qsm_laplacian_combine, qsm_romeo_B0, qsm_mask_filled

end # module
