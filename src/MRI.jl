module MRI

using Statistics

include("NIfTI_mod.jl")
include("Utility.jl")
include("hussein3d.jl")
include("laplacianunwrapping.jl")

export Data,
        readphase,
        savenii,
        createniiforwriting,
        getrobustmask, robustmask!,
        NIVolume,
        niread,
        niwrite,
        write_emptynii,
        hussein3d, hussein3d!,
        laplacianunwrap, laplacianunwrap!


end # module
