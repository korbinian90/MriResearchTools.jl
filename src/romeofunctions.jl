const romeo = unwrap # access unwrap function via alias romeo
const romeo! = unwrap!

"""
    mask_from_voxelquality(qmap::AbstractArray, threshold=0.5)

Creates a mask from a quality map. Another option is to use `robustmask(qmap)`

# Examples
```julia-repl
julia> qmap = romeovoxelquality(phase_3echo; TEs=[1,2,3]);
julia> mask = mask_from_voxelquality(qmap);
```

See also [`romeovoxelquality`](@ref), [`romeo`](@ref), [`robustmask`](@ref)
"""
function mask_from_voxelquality(qmap::AbstractArray, threshold=0.2, th2=0.8)
    qmap_bin = qmap .> threshold # NaN defaults to false (0)
    qmap_bin = fill_holes(qmap_bin)
    return gaussiansmooth3d(qmap_bin, (1,1,1)) .> th2
end

function fill_holes(mask; max_hole_size=length(mask) / 20)
    return .!imfill(.!mask, (1, max_hole_size)) # fills all holes up to max_hole_size (uses 6 connectivity as default for 3D)
end

function get_largest_connected_region(mask)
    labels = label_components(mask)
    return labels .== argmax(countmap(labels[labels .!= 0]))
end

function brain_mask(mask)
    bm = get_largest_connected_region(mask)
    return bm
end

function ROMEO.calculateweights(phase::AbstractArray{T,4}; TEs, template=2, p2ref=1, keyargs...) where T
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = phase[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag) && size(args[:mag], 4) > 1
        args[:mag] = args[:mag][:,:,:,template]
    end
    return ROMEO.calculateweights(view(phase,:,:,:,template); args...)
end

"""
    calculateB0_unwrapped(unwrapped_phase, mag, TEs)

Calculates B0 in [Hz] from unwrapped phase.
TEs in [ms].
The phase offsets have to be removed prior.

See also [`mcpc3ds`](@ref) and [`romeo`](@ref)
"""
function calculateB0_unwrapped(unwrapped_phase, mag, TEs)
    dims = 4
    TEs = to_dim(TEs, 4)
    return (1000 / 2π) * sum(unwrapped_phase .* mag .* mag .* TEs; dims) ./ sum(mag .* mag .* TEs.^2; dims) |> I -> dropdims(I; dims)
end

"""
    romeovoxelquality(phase::AbstractArray; keyargs...)

Calculates a quality for each voxel. The voxel quality can be used to create a mask.

# Examples
```julia-repl
julia> qmap = romeovoxelquality(phase_3echo; TEs=[1,2,3]);
julia> mask1 = mask_from_voxelquality(qmap);
julia> mask2 = robustmask(qmap);
```
     
Takes the same inputs as `romeo`/`unwrap`:
$(@doc unwrap)

See also [`mask_from_voxelquality`](@ref), [`romeo`](@ref), [`robustmask`](@ref)
""" 
function romeovoxelquality(phase; keyargs...)
    weights = ROMEO.calculateweights(phase; type=Float32, rescale=x->x, keyargs...) # [0;1]
    qmap = dropdims(sum(weights; dims=1); dims=1)
    qmap[2:end,:,:] .+= weights[1,1:end-1,:,:]
    qmap[:,2:end,:] .+= weights[2,:,1:end-1,:]
    qmap[:,:,2:end] .+= weights[3,:,:,1:end-1]
    return qmap ./ 6 # [0;1]
end
