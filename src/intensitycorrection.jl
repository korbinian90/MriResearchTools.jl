const DEBUG_PATH = "F:/MRI/Analysis/debug_hom"

"""
    makehomogeneous(mag::NIVolume; σ_mm=7, nbox=15)

Homogeneity correction for NIVolume from NIfTI files.

###  Keyword arguments:

- `σ_mm`: σ size for smoothing to obtain bias field. Takes NIfTI voxel size into account
- `nbox`: Number of boxes in each dimension for the box-segmentation step.

"""
function makehomogeneous(mag::NIVolume, datatype=eltype(mag); σ_mm=7, nbox=15)
    return makehomogeneous!(datatype.(mag); σ=mm_to_vox(σ_mm, mag), nbox=nbox)
end

"""
    makehomogeneous(mag; σ, nbox=15)

Homogeneity correction of 3D arrays. 4D volumes are corrected using the first 3D volume to
obtain the bias field.

###  Keyword arguments:

- `σ`: σ size in voxel for smoothing to obtain bias field. (mandatory)
- `nbox`: Number of boxes in each dimension for the box-segmentation step.

Larger σ-values make the bias field smoother, but might not be able to catch the
inhomogeneity. Smaller values can catch fast varying inhomogeneities but new inhomogeneities
might be created. The stronger the bias field, the more boxes are required for segmentation.
With too many boxes, it can happen that big darker structures are captured and appear
overbrightened.

"""
function makehomogeneous(mag, datatype=eltype(mag); σ, nbox=15)
    return makehomogeneous!(datatype.(mag); σ=σ, nbox=nbox)
end
function makehomogeneous!(mag; σ, nbox=15)
    lowpass = getsensitivity(mag; σ=σ, nbox=nbox)
    if eltype(mag) <: AbstractFloat
        mag ./= lowpass
    else # Integer doesn't support NaN
        lowpass[isnan.(lowpass) .| (lowpass .<= 0)] .= typemax(eltype(lowpass))
        mag .= div.(mag, lowpass ./ 2048) .|> x -> min(x, typemax(eltype(mag)))
    end
    mag
end

function getpixdim(nii::NIVolume)
    pixdim = nii.header.pixdim[2:(1+ndims(nii))]
    if all(pixdim .== 1)
        println("Warning! All voxel dimensions are 1 in NIfTI header, maybe they are wrong.")
    end
    return pixdim
end

mm_to_vox(mm, nii::NIVolume) = mm_to_vox(mm, getpixdim(nii))
mm_to_vox(mm, pixdim) = mm ./ pixdim


function getsensitivity(mag::NIVolume, datatype=eltype(mag); kw...)
    return getsensitivity(datatype.(mag), getpixdim(mag); kw...)
end
function getsensitivity(mag, pixdim; σ_mm=7, nbox=15)
    return getsensitivity(mag; σ=mm_to_vox(σ_mm, pixdim), nbox=nbox)
end
function getsensitivity(mag; σ, nbox=15)
    # segmentation
    firstecho = view(mag,:,:,:,1)
    @debug savenii(firstecho, "mag", DEBUG_PATH)
    mask = robustmask(firstecho)
    @debug savenii(mask, "mask", DEBUG_PATH)
    segmentation = boxsegment(firstecho, mask, nbox)
    @debug savenii(segmentation, "segmentation", DEBUG_PATH)
    # smoothing
    σ1, σ2 = getsigma(σ)
    lowpass = gaussiansmooth3d(firstecho, σ1; mask=segmentation, nbox=8)
    @debug savenii(lowpass, "lowpass_after_it", DEBUG_PATH)
    fillandsmooth!(lowpass, mean(firstecho[mask]), σ2)
    @debug savenii(lowpass, "lowpass_after_fillsmooth", DEBUG_PATH)

    return lowpass
end

# split sigma in two parts
function getsigma(σ)
    factorfinalsmoothing = 0.7
    σ1 = sqrt(1 - factorfinalsmoothing^2) .* σ
    σ2 = factorfinalsmoothing .* σ
    return σ1, σ2
end

function fillandsmooth!(lowpass, stablemean, σ2)
    stablethresh = stablemean / 4
    lowpassmask = (lowpass .< stablethresh) .| isnan.(lowpass) .| (lowpass .> 10 * stablemean)
    lowpass[lowpassmask] .= 3 * stablemean
    lowpassweight = 1.2 .- lowpassmask
    gaussiansmooth3d!(lowpass, σ2; weight = lowpassweight)
end

#threshold(image) = threshold(image, robustmask(image))
function threshold(image, mask; width=0.1)
    m = try quantile(skipmissing(image[mask]), 0.9) catch; 0 end
    return ((1 - width) * m .< image .< (1 + width) * m) .& mask
end

function boxsegment!(image::AbstractArray{<:AbstractFloat}, mask, nbox)
    image[boxsegment(image, mask, nbox)] .= NaN
    return image
end

function boxsegment(image, mask, nbox)
    N = size(image)
    dim = ndims(image)
    boxshift = ceil.(Int, N ./ nbox)

    segmented = zeros(UInt8, size(mask))
    for center in Iterators.product([1:boxshift[i]:N[i] for i in 1:dim]...)
        boxside(d) = max(1, center[d] - boxshift[d]):min(center[d] + boxshift[d], N[d])
        I = CartesianIndices(ntuple(boxside, dim))
        segmented[I] .+= threshold(image[I], mask[I])
    end
    return segmented .* mask .>= 2
end
