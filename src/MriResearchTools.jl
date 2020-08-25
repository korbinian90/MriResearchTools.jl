module MriResearchTools

using FFTW
using Interpolations
using NIfTI
using ROMEO
using Statistics
using NaNMath
using Plots

include("utility.jl")
include("smoothing.jl")
include("laplacianunwrapping.jl")
include("intensitycorrection.jl")
include("VSMbasedunwarping.jl")
include("methods.jl")
include("snr.jl")
include("tissue.jl")
include("scanner.jl")
include("snr_plots.jl")

romeo = unwrap # access unwrap function via alias romeo

export Data,
        readphase, readmag, niread,
        header,
        savenii,
        estimatenoise,
        robustmask, robustmask!,
        robustrescale,
        #combine_echoes,
        getHIP,
        laplacianunwrap, laplacianunwrap!,
        getVSM,
        unwarp,
        thresholdforward,
        gaussiansmooth3d!, gaussiansmooth3d,
        makehomogeneous!, makehomogeneous,
        getsensitivity,
        getscaledimage,
        estimatequantile,
        RSS,
        unwrap, unwrap!, romeo,
        unwrap_individual, unwrap_individual!,
        homodyne, homodyne!,
        to_dim,
        snr_graph,
        dutycycle_table

end # module
