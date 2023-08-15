module QSMExt

using MriResearchTools
using Statistics
import QSM: ismv, lbv, pdf, sharp, vsharp, nltv, rts, tikh, tkd, tsvd, tv

include("QSM_common.jl")

const γ = 267.52

function qsm(phase::AbstractArray, mask, TE, vsz; bfc_mask=mask, B0=3, bfc_algo=vsharp, qsm_algo=rts, unwrapping=laplacianunwrap, bdir=(0,0,1), kw...)
    vsz = Tuple(vsz)
    uphas = unwrapping(phase)
    uphas .*= inv(B0 * γ * TE) # convert units
    fl = bfc_algo(uphas, bfc_mask, vsz) # remove non-harmonic background fields

    # some background field correction methods require a mask update
    if fl isa Tuple
        fl, mask2 = fl
        mask = mask .& mask2
    end

    x = qsm_algo(fl, mask, vsz; bdir, kw...)
    return x
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask, vsz; bfc_mask=mask, B0=3, bfc_algo=vsharp, qsm_algo=rts, bdir=(0,0,1), kw...)
    scaled = B0_map .* (2π / (B0 * γ))
    fl = bfc_algo(scaled, bfc_mask, vsz)
    if fl isa Tuple
        fl, mask2 = fl
        mask = mask .& mask2
    end
    x = qsm_algo(fl, mask, vsz; bdir, kw...)
    return x
end

end
