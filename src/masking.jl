
function robustmask!(image; maskedvalue=if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!robustmask(image)] .= maskedvalue
    image
end

function robustmask(weight::AbstractArray; factor=1, threshold=nothing)
    if threshold isa Nothing
        μ, σ = estimatenoise(weight)
        m = mean(filter(isfinite, weight[weight .> 5σ]))
        threshold = minimum((5σ, m/2))
    end
    mask = weight .> (threshold * factor)
    # remove small holes and minimally grow
    boxsizes=[[5] for i in 1:ndims(weight)]
    mask = gaussiansmooth3d(mask; nbox=1, boxsizes) .> 0.4
    mask = fill_holes(mask)
    boxsizes=[[3,3] for i in 1:ndims(weight)]
    mask = gaussiansmooth3d(mask; nbox=2, boxsizes) .> 0.6
    return mask
end

"""
    robustmask(weight::AbstractArray; factor=1, threshold=nothing)

Creates a mask from an intensity/weight images by estimating a threshold and hole filling.
It assumes that at least one corner is without signal and only contains noise.
The automatic threshold is multiplied with `factor`.

# Examples
```julia-repl
julia> mask1 = robustmask(mag); # Using magnitude
julia> mask2 = phase_based_mask(phase); # Using phase
julia> mask3 = robustmask(romeovoxelquality(phase; mag)); # Using magnitude and phase
julia> brain = brain_mask(robustmask(romeovoxelquality(phase; mag); threshold=0.9));
```

See also [`ROMEO.romeovoxelquality`](@ref), [`phase_based_mask`](@ref), [`brain_mask`](@ref)
"""
robustmask, robustmask!

function sphere(radius, dim=3)
    len = 2radius + 1
    arr = OffsetArray(falses(repeat([len], dim)...), repeat([-radius:radius], dim)...)
    for I in CartesianIndices(arr)
        arr[I] = sqrt(sum(Tuple(I).^2)) < radius
    end
    return arr
end

"""
    phase_based_mask(phase; filter=true, threshold=1.0)

Creates a mask from a phase image.
Morphological filtering is activated by default.
To return the mask before thresholding pass `threshold=nothing`.

# Examples
```julia-repl
julia> phase_mask = phase_based_mask(phase);
```

See also [`ROMEO.romeovoxelquality`](@ref), [`ROMEO.romeo`](@ref), [`robustmask`](@ref), [`brain_mask`](@ref)

Original MATLAB algorithm:
    se=strel('sphere',6);
    L=del2(sign(wr));
    test=convn(abs(L),se.Neighborhood,'same');
    PB=mask{1}.*(test<500);
    PB=imclose(PB,se);
    mask{2}=round(imopen(PB,se));
"""
function phase_based_mask(phase; filter=true, threshold=1.0)
    strel = sphere(6, ndims(phase))
    laplacian = imfilter(sign.(phase), Kernel.Laplacian(1:ndims(phase), ndims(phase)))
    test = imfilter(abs.(laplacian), strel)
    if isnothing(threshold)
        return test * 500 * 6
    end
    PB = test .< (500 * 6 * threshold)
    if filter
        PB = LocalFilters.closing(PB, strel)
        PB = LocalFilters.opening(PB, strel)
    end
    return PB
end

"""
    mask_from_voxelquality(qmap::AbstractArray, threshold=:auto)

Creates a mask from a quality map. Another option is to use `robustmask(qmap)`

# Examples
```julia-repl
julia> qmap = romeovoxelquality(phase_3echo; TEs=[1,2,3]);
julia> mask = mask_from_voxelquality(qmap);
```

See also [`ROMEO.romeovoxelquality`](@ref), [`ROMEO.romeo`](@ref), [`robustmask`](@ref), [`brain_mask`](@ref)
"""
mask_from_voxelquality = robustmask

function fill_holes(mask; max_hole_size=length(mask) / 20)
    return .!imfill(.!mask, (1, max_hole_size)) # fills all holes up to max_hole_size (uses 6 connectivity as default for 3D)
end

function get_largest_connected_region(mask)
    labels = label_components(mask)
    return labels .== argmax(countmap(labels[labels .!= 0]))
end

"""
    brain_mask(mask)

Tries to extract the brain from a mask with skull and a gap between brain and skull.

# Examples
```julia-repl
julia> mask = robustmask(mag)
julia> brain = brain_mask(mask)
```

See also [`robustmask`](@ref), [`ROMEO.romeovoxelquality`](@ref), [`phase_based_mask`](@ref)
"""
function brain_mask(mask, strength=7)
    # set border to false
    shrink_mask = copy(mask)
    if ndims(shrink_mask) == 3 && all(size(shrink_mask) .> 5)
        shrink_mask[:,:,[1,end]] .= false
        shrink_mask[[1,end],:,:] .= false
        shrink_mask[:,[1,end],:] .= false
    end

    boxsizes=[[strength] for i in 1:ndims(shrink_mask)]
    smoothed = gaussiansmooth3d(shrink_mask; nbox=1, boxsizes)
    shrink_mask2 = smoothed .> 0.7

    brain_mask = get_largest_connected_region(shrink_mask2)

    # grow brain mask
    boxsizes=[[strength,strength] for i in 1:ndims(shrink_mask2)]
    smoothed = gaussiansmooth3d(brain_mask; nbox=2, boxsizes)
    brain_mask = smoothed .> 0.2
    return brain_mask .& mask
end
