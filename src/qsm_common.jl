# Backend-independent QSM helpers.
#
# These used to live in `ext/QSM_common.jl`, which is `include`d into *both*
# `QSMExt` and `QuantitativeSusceptibilityMappingTGVExt`. Anything defined there
# that extends a `MriResearchTools` function is therefore defined twice with the
# same signature, and the surviving method depends on extension load order.
# Nothing below needs a QSM backend, so it belongs in the package proper - that
# removes the duplicate definitions and makes the functions usable without
# loading a backend at all.

function qsm_mask_filled(phase::AbstractArray, vsz; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma_in_mm=[5,5,5])
    smooth_sigma = smooth_sigma_in_mm ./ vsz
    qsm_mask_filled(phase; quality_thresh, smooth_thresh, smooth_sigma)
end

function qsm_mask_filled(phase::AbstractArray; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma=[5,5,5])
    mask_small = (romeovoxelquality(phase) .> quality_thresh) .& (phase .!= 0)
    mask_filled = gaussiansmooth3d(mask_small, smooth_sigma; padding=true) .> smooth_thresh
    if sum(mask_filled) == 0
        @warn("QSM Mask cannot be created (maybe too low quality data?)")
        mask_filled .= true
    end
    return mask_filled
end

"""
    qsm_mask_filled(phase; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma=[5,5,5])

    qsm_mask_filled(phase, vsz; quality_thresh=0.5, smooth_thresh=0.5, smooth_sigma_in_mm=[5,5,5])

Creates a filled mask suitable for QSM from the ROMEO voxel quality of `phase`.
Voxels below `quality_thresh` are excluded, then the mask is smoothed and
re-thresholded at `smooth_thresh` to close holes.

Given the voxel size `vsz`, `smooth_sigma_in_mm` is used instead and converted to
voxels, which makes the result independent of the acquisition resolution.

# Examples
```julia-repl
julia> mask = qsm_mask_filled(phase[:,:,:,1]);
```

See also [`robustmask`](@ref), [`brain_mask`](@ref)
"""
qsm_mask_filled

# Echo-combination weighted by mag^2 * TE^2 (the phase-SNR weighting).
function weighted_average(image::AbstractArray{<:Number,4}, mag, TEs)
    w_sum = sum(image[:,:,:,i] .* mag[:,:,:,i].^2 .* TEs[i]^2 for i in axes(image,4)) ./ sum(mag[:,:,:,i].^2 .* TEs[i]^2 for i in axes(image,4))
    nans = isnan.(w_sum)
    w_sum[nans] .= mean(image[nans, i] for i in axes(image, 4))
    return w_sum
end

function weighted_average(images, mag, TEs)
    weighted_average(cat(images...; dims=4), mag, TEs)
end

# Laplacian-unwrapped echo combination, output in [Hz].
function laplacian_combine(phase::AbstractArray, mag, TEs; type=:weighted_average)
    if type == :average
        return mean(laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)) ./ 2π
    elseif type == :weighted_average
        return weighted_average((laplacianunwrap(phase[:,:,:,i]) ./ TEs[i] for i in axes(phase, 4)), mag, TEs) ./ 2π
    end
    error("laplacian_combine type '$type' is not defined!")
end
