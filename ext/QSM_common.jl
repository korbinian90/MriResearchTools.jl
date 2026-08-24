# Backend-dependent QSM entry points.
#
# This file is `include`d into BOTH `QSMExt` and
# `QuantitativeSusceptibilityMappingTGVExt`, because every function here needs the
# module-local `qsm` / `qsm_B0` of the backend it is compiled into. That is also
# why loading both backends at once is not supported: the methods below have
# identical signatures in both extensions, so the one that loads last wins.
# Backend-independent helpers live in `src/qsm_common.jl` instead, so they are
# defined exactly once.

function MriResearchTools.qsm_average(phase::AbstractArray, mag::AbstractArray, mask::AbstractArray, TEs, vsz; kw...)
    MriResearchTools.weighted_average((qsm(phase[:,:,:,i], mask, TEs[i], vsz; kw...) for i in axes(phase, 4)), mag, TEs)
end

function MriResearchTools.qsm_romeo_B0(phase::AbstractArray, mag, mask, TEs, res; kw...)
    if size(phase, 4) > 1
        phase, _ = mcpc3ds(phase, mag; TEs)
    end
    unwrapped = romeo(phase; TEs, mag, individual=true, correctglobal=true)
    B0_map = calculateB0_unwrapped(unwrapped, mag, TEs .* 1e3, :average)
    if isnothing(mask)
        mask = qsm_mask_filled(B0_map .* 2pi .* 40, res) # needs to be scaled to a reasonable TE=40 for masking
    end
    haskey(kw, :save) && kw[:save](mask, "qsm_mask")
    return qsm_B0(B0_map, mask, res; kw...)
end

function MriResearchTools.qsm_laplacian_combine(phase::AbstractArray, mag, mask, TEs, res; laplacian_combine_type=:weighted_average, kw...)
    local_B0 = MriResearchTools.laplacian_combine(phase, mag, TEs; type=laplacian_combine_type)
    return qsm_B0(local_B0, mask, res; kw...)
end
