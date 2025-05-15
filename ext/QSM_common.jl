function MriResearchTools.qsm_average(phase::AbstractArray, mag::AbstractArray, mask::AbstractArray, TEs, vsz; kw...)
    weighted_average((qsm(phase[:,:,:,i], mask, TEs[i], vsz; kw...) for i in axes(phase, 4)), mag, TEs)
end

function MriResearchTools.qsm_romeo_B0(phase::AbstractArray, mag, mask, TEs, res; kw...)
    if size(phase, 4) > 1
        phase, _ = mcpc3ds(phase, mag; TEs)
    end
    unwrapped = romeo(phase; TEs, mag, individual=true, correct_global=true)
    B0_map = calculateB0_unwrapped(unwrapped, mag, TEs .* 1e3, :average)
    if isnothing(mask)
        mask = qsm_mask_filled(B0_map .* 2pi .* 40, res) # needs to be scaled to a reasonable TE=40 for masking
    end
    haskey(kw, :save) && kw[:save](mask, "qsm_mask")
    return qsm_B0(B0_map, mask, res; kw...)
end

function MriResearchTools.qsm_laplacian_combine(phase::AbstractArray, mag, mask, TEs, res; laplacian_combine_type=:weighted_average, kw...)
    local_B0 = laplacian_combine(phase, mag, TEs; type=laplacian_combine_type)
    return qsm_B0(local_B0, mask, res; kw...)
end

# output in [Hz]
function laplacian_combine(phase::AbstractArray, mag, TEs; type=:weighted_average)
    if type == :average
        return mean(laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)) ./ 2π
    elseif type == :weighted_average
        return weighted_average((laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)), mag, TEs) ./ 2π
    end
end

function MriResearchTools.qsm_mask_filled(phase::AbstractArray, vsz; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma_in_mm=[5,5,5])
    smooth_sigma = smooth_sigma_in_mm ./ vsz
    qsm_mask_filled(phase; quality_thresh, smooth_thresh, smooth_sigma)
end

function MriResearchTools.qsm_mask_filled(phase::AbstractArray; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma=[5,5,5])
    mask_small = (romeovoxelquality(phase) .> quality_thresh) .& (phase .!= 0)
    mask_filled = gaussiansmooth3d(mask_small, smooth_sigma; padding=true) .> smooth_thresh
    if sum(mask_filled) == 0
        @warn("QSM Mask cannot be created (maybe too low quality data?)")
        mask_filled .= true
    end
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
