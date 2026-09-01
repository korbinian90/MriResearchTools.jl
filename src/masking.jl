
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
julia> # Multi-echo phase additionally needs the echo times:
julia> mask4 = robustmask(romeovoxelquality(phase_3echo; mag=mag_3echo, TEs=[1,2,3]));
julia> brain = brain_mask(robustmask(romeovoxelquality(phase; mag); threshold=0.9));
```

See also [`brain_mask`](@ref)
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

See also [`robustmask`](@ref), [`brain_mask`](@ref)
"""
const mask_from_voxelquality = robustmask

# Connected components, 6-connectivity in 3D (4 in 2D, the N-dimensional
# diamond in general). This used to be ImageMorphology's `label_components` and
# `imfill`, which were the only two things this package took from that
# dependency; ImageMorphology in turn is the sole reason the whole
# LoopVectorization stack (15 packages) was in the manifest, which is a large
# part of the compiled binaries and keeps the stack off Julia 1.12. The
# partition produced here is identical to ImageMorphology's.
#
# Union-find over linear indices, unioning each true voxel with its already
# visited neighbour in each dimension, so one pass suffices.
function connected_components(mask::AbstractArray{Bool,N}) where N
    parent = collect(1:length(mask))
    function find(x)
        root = x
        while parent[root] != root
            root = parent[root]
        end
        while parent[x] != root # path compression
            parent[x], x = root, parent[x]
        end
        return root
    end
    lin = LinearIndices(mask)
    for I in CartesianIndices(mask)
        mask[I] || continue
        for d in 1:N
            J = I - CartesianIndex(ntuple(k -> k == d ? 1 : 0, N))
            (checkbounds(Bool, mask, J) && mask[J]) || continue
            a, b = find(lin[I]), find(lin[J])
            a != b && (parent[a] = b)
        end
    end
    labels = zeros(Int, size(mask))
    sizes = Int[]
    label_of_root = Dict{Int,Int}()
    for I in CartesianIndices(mask)
        mask[I] || continue
        root = find(lin[I])
        label = get(label_of_root, root, 0)
        if label == 0
            push!(sizes, 0)
            label = label_of_root[root] = length(sizes)
        end
        labels[I] = label
        sizes[label] += 1
    end
    return labels, sizes
end

function fill_holes(mask; max_hole_size=length(mask) / 20)
    labels, sizes = connected_components(.!mask) # holes are background components
    filled = copy(mask)
    for I in CartesianIndices(mask)
        label = labels[I]
        label != 0 && sizes[label] <= max_hole_size && (filled[I] = true)
    end
    return filled
end

function get_largest_connected_region(mask)
    labels, sizes = connected_components(mask)
    isempty(sizes) && return falses(size(mask))
    return labels .== argmax(sizes)
end

"""
    brain_mask(mask)

Tries to extract the brain from a mask with skull and a gap between brain and skull.

# Examples
```julia-repl
julia> mask = robustmask(mag)
julia> brain = brain_mask(mask)
```

See also [`robustmask`](@ref)
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
