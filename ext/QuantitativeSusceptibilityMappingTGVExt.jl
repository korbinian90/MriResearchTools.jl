module QuantitativeSusceptibilityMappingTGVExt

using MriResearchTools
using Statistics
using QuantitativeSusceptibilityMappingTGV

include("QSM_common.jl")

function qsm(phase::AbstractArray, mask, TE, res; B0, kw...)
    qsm_tgv(phase, mask, res; TE, fieldstrength=B0, kw...)
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask, res; B0, kw...)
    qsm_tgv(B0_map, mask, res; TE=1e-3, fieldstrength=B0, kw...)
end

end
