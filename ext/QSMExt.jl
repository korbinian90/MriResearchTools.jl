module QSMExt

using MriResearchTools
using Statistics
import QSM: ismv, lbv, pdf, sharp, vsharp, nltv, rts, tikh, tkd, tsvd, tv

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

function MriResearchTools.qsm_average(phase::AbstractArray, mag, mask, TEs, vsz; kw...)
    weighted_average((qsm(phase[:,:,:,i], mask, TEs[i], vsz; kw...) for i in axes(phase, 4)), mag, TEs)
end

function MriResearchTools.qsm_romeo_B0(phase::AbstractArray, mag, mask, TEs, vsz; kw...)
    phasecorr, _ = mcpc3ds(phase, mag; TEs)
    unwrapped = romeo(phasecorr; TEs, mag)
    B0_map = calculateB0_unwrapped(unwrapped, mag, TEs .* 1000)
    return qsm_B0(B0_map, mask, vsz; kw...)
end

function MriResearchTools.qsm_laplacian_combine(phase::AbstractArray, mag, mask, TEs, vsz; laplacian_combine_type=:weighted_average, kw...)
    local_B0 = laplacian_combine(phase, mag, TEs; type=laplacian_combine_type)
    return qsm_B0(local_B0, mask, vsz; kw...)
end

# output in [Hz]
function laplacian_combine(phase::AbstractArray, mag, TEs; type=:weighted_average)
    if type == :average
        return mean(laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)) ./ 2π
    elseif type == :weighted_average
        return weighted_average((laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)), mag, TEs) ./ 2π
    end
end

function MriResearchTools.qsm_mask_filled(phase::AbstractArray; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma=[5,5,5])
    mask_small = romeovoxelquality(phase) .> quality_thresh
    mask_filled = gaussiansmooth3d(mask_small, smooth_sigma) .> smooth_thresh
    return mask_filled
end

function weighted_average(image::AbstractArray{<:Number,4}, mag, TEs)
    w_sum = sum(image[:,:,:,i] .* mag[:,:,:,i].^2 .* TEs[i]^2 for i in axes(image,4)) ./ sum(mag[:,:,:,i].^2 .* TEs[i]^2 for i in axes(image,4))
    nans = isnan.(w_sum)
    w_sum[nans] .= mean(image[nans, i] for i in axes(image, 4))
    return w_sum
end

function weighted_average(images, mag, TEs)
    weighted_average(cat(images...; dims=4), mag, TEs)
end

end
