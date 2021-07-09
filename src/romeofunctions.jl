romeo = unwrap # access unwrap function via alias romeo
romeo! = unwrap!

function mask_from_voxelquality(qmap, threshold)
    qmap_bin = qmap .> threshold # NaN defaults to false (0)
    max_hole_size = length(qmap) / 20
    qmap_bin = .!imfill(.!qmap_bin, (1, max_hole_size)) # fills all holes up to max_hole_size (uses 6 connectivity as default)
    return gaussiansmooth3d(qmap_bin) .> 0.8 # hardcoded final threshold
end

function ROMEO.calculateweights(phase::AbstractArray{T,4}; TEs, template=2, p2ref=1, keyargs...) where T
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = phase[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag)
        args[:mag] = args[:mag][:,:,:,template]
    end
    return ROMEO.calculateweights(view(phase,:,:,:,template); args...)
end

# calculates B0 in [Hz]
function calculateB0_unwrapped(unwrapped_phase, mag, TEs)
    dims = 4
    TEs = to_dim(TEs, 4)
    return (1000 / 2π) * sum(unwrapped_phase .* mag .* mag .* TEs; dims=dims) ./ sum(mag .* mag .* TEs.^2; dims=dims) |> I -> dropdims(I; dims=dims)
end

"""
Calculates a quality for each voxel. Takes the same inputs as ROMEO
"""
function romeovoxelquality(phase; keyargs...)
    weights = ROMEO.calculateweights(phase; type=Float32, rescale=x->x, keyargs...)
    return dropdims(sum(weights; dims=1); dims=1) ./ 3 # [0;1]
end
