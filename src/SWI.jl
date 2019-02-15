function calculateSWI(mag, phase, TEs; σ = [2, 2, 1], power = 1, mode = :positive)
    getswimag(mag) .* getswiphase(phase, mag, TEs, σ, mode).^power
end

getswimag(mag; keyargs...) = makehomogeneous!(combine_echoes(mag); keyargs...)

function getswiphase(phase, mag, TEs, σ, mode; keyargs...)
    mask = getrobustmask(view(mag,:,:,:,1))
    if haskey(keyargs, :writedir)
        savenii(mask, "maskforphase", keyargs[:writedir])
    end

    combined = getcombinedphase(phase, mag, σ, TEs, mask; keyargs...)
    scaleandthreshold!(combined, mask; mode = mode)
end

function getcombinedphase(phase, mag, σ, TEs, mask; writedir = nothing, unwrapping = :kunwrap)
    echodim = ndims(phase)
    # reshape TEs to the highest dimension, eg. "size(TEs) = (1, 1, 1, 3)"
    TEs = reshape(TEs, ones(Int, echodim - 1)..., length(TEs))

    if unwrapping == :laplacian
        combined = sum(mag[:,:,:,iEco] .* unwrap(rescale(view(phase,:,:,:,iEco), -π, π; datatype = Float32)) for iEco in 1:size(phase, echodim))
    elseif unwrapping == :kunwrap
        combined = dropdims(sum(mag .* kunwrap(phase, mag, TEs; mask = mask); dims = echodim); dims = echodim)
    else
        error("Unwrapping $unwrapping not defined!")
    end

    #end
    if writedir != nothing savenii(combined, "$writedir/combinedphase.nii") end

    magdiv = dropdims(sum(mag .* Float32.(TEs); dims = echodim); dims = echodim)
    if writedir != nothing savenii(magdiv, "$writedir/magdiv.nii") end

    combined ./= magdiv
    if writedir != nothing savenii(combined, "$writedir/scaledphase.nii") end

    combined .-= gaussiansmooth3d(combined, σ; mask = mask, dims = 1:2)
    if writedir != nothing savenii(combined, "$writedir/filteredphase.nii") end

    #TODO change order
    combined
end


function highpassunwrap(phase, σ, mask; name)
    # TODO Float64 or Float32?
    unwrapped = unwrap(rescale(phase, -π, π; datatype = Float32))
    #savenii(unwrapped, "/data/korbi/data/julia/SWI_temp/unwrapped_$iEco.nii")
    #TODO check if mask is neccessary
    filtered = unwrapped .- gaussiansmooth3d(unwrapped, σ; mask = mask, dims = 1:3)
    if name != nothing
        savenii(filtered, name)
    end
    filtered
end

function scaleandthreshold!(swiphase, mask; mode = :tanh)
    swiphase[.!mask] .= 0

    pos = swiphase .> 0
    if mode == :tanh
        f(x) = 0.2 + 0.4(1 + tanh(100(-x + 0.015)))
        swiphase .= f.(swiphase)
    elseif mode == :positive
        swiphase[.!pos] .= 1
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0)
    elseif mode == :negative
        swiphase[pos] .= 1
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1)
    elseif mode == :triangular
        swiphase[pos] .= robustrescale(swiphase[pos], 1, 0)
        swiphase[.!pos] .= robustrescale(swiphase[.!pos], 0, 1)
    end
    swiphase[swiphase .< 0] .= 0
    swiphase
end
