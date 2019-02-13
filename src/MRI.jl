module MRI

using Statistics

include("NIfTI_mod.jl")
include("Utility.jl")

export Data,
        readphase,
        savenii,
        createniiforwriting,
        getrobustmask, robustmask!,
        NIVolume,
        niread,
        niwrite,
        write_emptynii


end # module
