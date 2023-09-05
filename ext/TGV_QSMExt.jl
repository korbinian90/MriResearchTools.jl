module TGV_QSMExt

using MriResearchTools
using Statistics
using TGV_QSM

include("QSM_common.jl")

function qsm(phase::AbstractArray, mask, TE, res; kw...)
    qsm_tgv(phase, mask, res; TE, kw...)
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask, res; kw...)
    qsm_tgv(B0_map, mask, res; TE=1e-3, fieldstrength, iterations, erosions)
end

function MriResearchTools.qsm_average(phase::AbstractArray, mag, mask, TEs, res; kw...)
    weighted_average((qsm_tgv(phase[:,:,:,i], mask, res; TE=TEs[i], kw...) for i in axes(phase, 4)), mag, TEs)
end

end
