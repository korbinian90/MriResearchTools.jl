const DEBUG_PATH = "F:/MRI/Analysis/debug_hom"

function makehomogeneous(mag::NIVolume, datatype=eltype(mag); σ_mm=7, kw...)
    return makehomogeneous!(datatype.(mag); σ=mm_to_vox(σ_mm, mag), kw...)
end
function makehomogeneous(mag, datatype=eltype(mag); σ, kw...)
    return makehomogeneous!(datatype.(mag); σ=σ, kw...)
end
function makehomogeneous!(mag; σ, nbox=15)
    lowpass = getsensitivity(mag; σ=σ)
    if eltype(mag) <: AbstractFloat
        mag ./= lowpass
    else # Integer doesn't support NaN
        lowpass[isnan.(lowpass) .| (lowpass .<= 0)] .= typemax(eltype(lowpass))
        mag .= div.(mag, lowpass ./ 2048) .|> x -> min(x, typemax(eltype(mag)))
    end
    mag
end

function mm_to_vox(mm, nii::NIVolume)
    pixdim = nii.header.pixdim[2:(1+ndims(nii))]
    if all(pixdim .== 1)
        println("Warning! All voxel dimensions are 1 in NIfTI header, maybe they are wrong.")
    end
    return mm_to_vox(mm, pixdim)
end
mm_to_vox(mm, pixdim) = mm ./ pixdim

function getsensitivity(mag::NIVolume, datatype=eltype(mag); σ_mm=7, kw...)
    return getsensitivity(datatype.(mag); σ=mm_to_vox(σ_mm, mag), kw...)
end
function getsensitivity(mag; σ, kw...)
    σ1, σ2 = getsigma(σ)
    firstecho = view(mag,:,:,:,1)
    @debug savenii(firstecho, "mag", DEBUG_PATH)

    mask = robustmask(firstecho)
    @debug savenii(mask, "mask", DEBUG_PATH)
    segmentation = boxsegment(firstecho, mask; kw...)
    @debug savenii(segmentation, "segmentation", DEBUG_PATH)
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

threshold(image) = threshold(image, robustmask(image))
function threshold(image, mask; width=0.1)
    m = try quantile(skipmissing(image[mask]), 0.9) catch; 0 end
    return ((1 - width) * m .< image .< (1 + width) * m) .& mask
end

function boxsegment!(image::AbstractArray{<:AbstractFloat}, mask; kw...)
    image[boxsegment(image, mask; kw...)] .= NaN
    return image
end

function boxsegment(image, mask; nboxes=15)
    N = size(image)
    dim = ndims(image)
    boxshift = ceil.(Int, N ./ nboxes)

    segmented = zeros(UInt8, size(mask))
    for center in Iterators.product([1:boxshift[i]:N[i] for i in 1:dim]...)
        boxside(d) = max(1, center[d] - boxshift[d]):min(center[d] + boxshift[d], N[d])
        I = CartesianIndices(ntuple(boxside, dim))
        segmented[I] .+= threshold(image[I], mask[I])
    end
    return segmented .* mask .>= 2
end
