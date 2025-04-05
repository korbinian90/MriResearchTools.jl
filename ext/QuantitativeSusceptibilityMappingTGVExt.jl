module QuantitativeSusceptibilityMappingTGVExt

using MriResearchTools
using Statistics
using QuantitativeSusceptibilityMappingTGV

include("QSM_common.jl")

function qsm(phase::AbstractArray, mask, TE, res; B0, kw...)
    qsm_tgv(phase, mask, res; TE, fieldstrength=B0, kw...)
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask::AbstractArray, res; B0, save=nothing, kw...)
    TE = 40 # ms
    map_scaled = B0_map .* 2pi .* TE # needs to be scaled to a reasonable TE for TGV    
    qsm_tgv(map_scaled, mask, res; TE=TE*1e-3, fieldstrength=B0, kw...)
end

end
