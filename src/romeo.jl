# multi echo unwrapping
function ROMEO.unwrap!(
    wrapped::AbstractArray{T, 4}; TEs=1:size(wrapped, 4), template=2, p2ref=1, keyargs...
) where {T <: AbstractFloat}
    args = Dict{Symbol, Any}(keyargs)
    args[:phase2] = wrapped[:,:,:,p2ref]
    args[:TEs] = TEs[[template, p2ref]]
    if haskey(args, :mag)
        args[:mag] = args[:mag][:,:,:,template]
    end
    unwrap!(view(wrapped,:,:,:,template); args...)
    for iEco in [(template-1):-1:1; (template+1):length(TEs)]
        iRef = if (iEco < template) iEco + 1 else iEco - 1 end
        wrapped[:,:,:,iEco] .= ROMEO.unwrapvoxel.(wrapped[:,:,:,iEco], wrapped[:,:,:,iRef] .* (TEs[iEco] / TEs[iRef]))
        #magf = if !haskey(keyargs, :mag) nothing else view(mag,:,:,:,iEco) end
        #unwrapfilter!(view(wrapped,:,:,:,iEco), magf) # TODO
    end
    return wrapped
end

function unwrapfilter!(phase, mag)
    smoothedphase = if mag == nothing
        gaussiansmooth3d(phase; boxsizes=[3,3,3], nbox=1)
    else
        # corresponds to weighted meanfilter with nbox = 1
        gaussiansmooth3d(phase; weight=Float32.(mag), boxsizes=[3,3,3], nbox=1)
    end
    phase .= unwrapvoxel.(phase, smoothedphase)
    return
end

unwrap_individual(wrapped; keyargs...) = unwrap_individual!(Float32.(wrapped); keyargs...)
function unwrap_individual!(
    wrapped::AbstractArray{T, 4}; TEs=1:size(wrapped, 4), keyargs...
) where {T <: AbstractFloat}
    for i in 1:length(TEs)
        e2 = (i == 1) ? 2 : (i - 1)
        echoes = [i, e2]
        args = Dict()
        if haskey(keyargs, :mag) args[:mag] = keyargs[:mag][:,:,:,i] end
        unwrap!(view(wrapped,:,:,:,i); phase2=wrapped[:,:,:,e2], TEs=TEs[echoes], args...)
    end
    return wrapped
end
