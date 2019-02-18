module MRI

using Statistics

include("NIfTI_mod.jl")
include("utility.jl")
include("smoothing.jl")
include("hussein3d.jl")
include("laplacianunwrapping.jl")
include("intensitycorrection.jl")
include("SWI.jl")
include("VSMbasedunwarping.jl")

export Data,
        readphase, readmag,
        savenii,
        createniiforwriting,
        getrobustmask, robustmask!,
        NIVolume,
        niread, niwrite,
        write_emptynii,
        getHIP,
        hussein3d, hussein3d!,
        laplacianunwrap, laplacianunwrap!,
        calculateSWI,
        getVSM,
        unwarp,
        thresholdforward,
        gaussiansmooth3d!, gaussiansmooth3d,
        makehomogeneous!, makehomogeneous,
        getscaledimage

end # module
