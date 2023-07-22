
function robustmask!(image; maskedvalue=if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!robustmask(image)] .= maskedvalue
    image
end

function robustmask(weight::AbstractArray; factor=1, threshold=nothing)
    if threshold isa Nothing
        w = sample(weight)
        q05, q15, q8, q99 = quantile.(Ref(w), (0.05, 0.15, 0.8, 0.99))
        high_intensity = mean(filter(isfinite, w[q8 .<= w .<= q99]))
        noise = mean(filter(isfinite, w[w .<= q15]))
        if noise > high_intensity/10
            noise = mean(filter(isfinite, w[w .<= q05]))
            if noise > high_intensity/10
                noise = 0 # no noise detected
            end
        end
        threshold = max(5noise, high_intensity/5)
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
