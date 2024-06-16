"""
    makehomogeneous(mag::NIVolume; sigma_mm=7, nbox=15)

Homogeneity correction for NIVolume from NIfTI files.

###  Keyword arguments:

- `sigma_mm`: sigma size for smoothing to obtain bias field. Takes NIfTI voxel size into account
- `nbox`: Number of boxes in each dimension for the box-segmentation step.

"""
function makehomogeneous(mag::NIVolume, datatype=eltype(mag); sigma_mm=7, nbox=15)
    return makehomogeneous!(datatype.(mag); sigma=mm_to_vox(sigma_mm, mag), nbox)
end

"""
    makehomogeneous(mag; sigma, nbox=15)

Homogeneity correction of 3D arrays. 4D volumes are corrected using the first 3D volume to
obtain the bias field.

###  Keyword arguments:

- `sigma`: sigma size in voxel for each dimension for smoothing to obtain bias field. (mandatory)
- `nbox`: Number of boxes in each dimension for the box-segmentation step.

Larger sigma-values make the bias field smoother, but might not be able to catch the
inhomogeneity. Smaller values can catch fast varying inhomogeneities but new inhomogeneities
might be created. The stronger the bias field, the more boxes are required for segmentation.
With too many boxes, it can happen that big darker structures are captured and appear
overbrightened.

Calculates the bias field using the `boxsegment` approach.
It assumes that there is a "main tissue" that is present in most areas of the object.
Published in [CLEAR-SWI](https://doi.org/10.1016/j.neuroimage.2021.118175).

See also [`getsensitivity`](@ref)
"""
makehomogeneous, makehomogeneous!

function makehomogeneous(mag, datatype=eltype(mag); sigma, nbox=15)
    return makehomogeneous!(datatype.(mag); sigma, nbox)
end
function makehomogeneous!(mag; sigma, nbox=15)
    lowpass = getsensitivity(mag; sigma, nbox)
    if eltype(mag) <: AbstractFloat
        mag ./= lowpass
    else # Integer doesn't support NaN
        lowpass[isnan.(lowpass).|(lowpass.<=0)] .= typemax(eltype(lowpass))
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


"""
    getsensitivity(mag; sigma, nbox=15)

    getsensitivity(mag, pixdim; sigma_mm=7, nbox=15)

    getsensitivity(mag::NIVolume, datatype=eltype(mag); sigma_mm=7, nbox=15)

Calculates the bias field using the `boxsegment` approach.
It assumes that there is a "main tissue" that is present in most areas of the object.
If not set, sigma_mm defaults to 7mm, with a maximum of 10% FoV. The sigma_mm value should 
correspond to the bias field, for a faster changing bias field this needs to be smaller.
Published in [CLEAR-SWI](https://doi.org/10.1016/j.neuroimage.2021.118175).

See also [`makehomogeneous`](@ref)
"""
function getsensitivity(mag::NIVolume, datatype=eltype(mag); kw...)
    return getsensitivity(datatype.(mag), getpixdim(mag); kw...)
end
function getsensitivity(mag, pixdim; sigma_mm=get_default_sigma_mm(mag, pixdim), nbox=15)
    return getsensitivity(mag; sigma=mm_to_vox(sigma_mm, pixdim), nbox)
end
function getsensitivity(mag; sigma, nbox=15)
    # segmentation
    firstecho = view(mag, :, :, :, 1)
    mask = robustmask(firstecho)
    segmentation = boxsegment(firstecho, mask, nbox)
    # smoothing
    sigma1, sigma2 = getsigma(sigma)
    lowpass = gaussiansmooth3d(firstecho, sigma1; mask=segmentation, nbox=8)
    fillandsmooth!(lowpass, mean(firstecho[mask]), sigma2)

    return lowpass
end

# Default is 7mm, but a maximum of 10% FoV
function get_default_sigma_mm(mag, pixdim)
    sigma_mm = zeros(min(ndims(mag), length(pixdim)))
    for i in eachindex(sigma_mm)
        sigma_mm[i] = pixdim[i] * size(mag, i)
    end
    sigma_mm = median(sigma_mm)
    sigma_mm = min(sigma_mm, 7)
    return sigma_mm
end

# split sigma in two parts
function getsigma(sigma)
    factorfinalsmoothing = 0.7
    sigma1 = sqrt(1 - factorfinalsmoothing^2) .* sigma
    sigma2 = factorfinalsmoothing .* sigma
    return sigma1, sigma2
end

function fillandsmooth!(lowpass, stablemean, sigma2)
    lowpassmask = (lowpass .< stablemean / 4) .| isnan.(lowpass) .| (lowpass .> 10 * stablemean)
    lowpass[lowpassmask] .= 3 * stablemean
    lowpassweight = 1.2 .- lowpassmask
    gaussiansmooth3d!(lowpass, sigma2; weight=lowpassweight)
end

#threshold(image) = threshold(image, robustmask(image))
function threshold(image, mask; width=0.1)
    m = try
        quantile(skipmissing(image[mask]), 0.9)
    catch
        0
    end
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
