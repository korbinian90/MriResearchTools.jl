module MriResearchTools

using FFTW
using Interpolations
using NIfTI
using ROMEO
using Statistics
using NaNMath

include("utility.jl")
include("smoothing.jl")
include("laplacianunwrapping.jl")
include("intensitycorrection.jl")
include("VSMbasedunwarping.jl")
include("romeo.jl")

export Data,
        readphase, readmag, niread,
        header,
        savenii,
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
        unwrap_individual, unwrap_individual!

end # module
