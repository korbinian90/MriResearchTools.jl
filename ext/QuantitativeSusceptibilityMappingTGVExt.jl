module QuantitativeSusceptibilityMappingTGVExt

using MriResearchTools
using Statistics
using QuantitativeSusceptibilityMappingTGV

include("QSM_common.jl")

# TGV sits below MriResearchTools in the dependency graph, so it cannot register
# its own references; this extension does it on its behalf. Both references are
# for the one method, so they share a heading.
const TGV = QuantitativeSusceptibilityMappingTGV
const TGV_VERSION = let v = isdefined(TGV, :PKG_VERSION) ? TGV.PKG_VERSION : pkgversion(TGV)
    v === nothing ? "unknown" : string(v)
end

function __init__()
    register_version!(TGV, TGV_VERSION)
    register_citation!(:tgv,
        """Langkammer, C., Bredies, K., Poser, B.A., Barth, M., Reishofer, G., Fan, A.P., Bilgic, B., Fazekas, F., Mainero, C., Ropele, S., 2015.
           Fast quantitative susceptibility mapping using 3D EPI and total generalized variation.
           NeuroImage 111, 622-630.
           https://doi.org/10.1016/j.neuroimage.2015.02.041""";
        label = "TGV QSM")
    register_citation!(:tgv_original,
        """Bredies, K., Ropele, S., Poser, B.A., Barth, M., Langkammer, C., 2014.
           Single-step quantitative susceptibility mapping using total generalized variation and 3D EPI.
           Proceedings of the 22nd Annual Meeting ISMRM, p. 604.""";
        label = "TGV QSM")
end

function qsm(phase::AbstractArray, mask, TE, res; B0, kw...)
    qsm_tgv(phase, mask, res; TE, fieldstrength=B0, kw...)
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask::AbstractArray, res; B0, save=nothing, kw...)
    TE = 40 # ms
    map_scaled = B0_map .* 2pi .* TE # needs to be scaled to a reasonable TE for TGV    
    qsm_tgv(map_scaled, mask, res; TE=TE*1e-3, fieldstrength=B0, kw...)
end

end
