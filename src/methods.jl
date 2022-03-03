"""
    homodyne(mag, phase)
    
    homodyne(mag, phase; dims, σ)

    homodyne(I; kwargs...)

Performs homodyne filtering via division of complex complex smoothing.
"""
homodyne, homodyne!

homodyne(mag, phase; kwargs...) = homodyne(mag .* exp.(1im .* phase); kwargs...)
homodyne(I; kwargs...) = homodyne!(copy(I); kwargs...)
homodyne!(I; dims=1:2, σ=8 .* ones(length(dims))) = I ./= gaussiansmooth3d(I, σ; dims)

"""
    NumART2star(image::AbstractArray{T,4}, TEs) where T

Performs T2* calculation on 4D-multi-echo magnitude data.
https://doi.org/10.1002/mrm.10283
"""
function NumART2star(image::AbstractArray{T,4}, TEs) where T
    neco = length(TEs)
    t2star(m) = (TEs[end]-TEs[1]) / 2(neco-1) * (m[1]+m[end]+2sum(m[2:end-1])) / (m[1]-m[end])
    return [t2star(image[I,:]) for I in CartesianIndices(size(image)[1:3])]
end

"""
    r2s_from_t2s(t2s) = 1000 ./ t2s

Converts from T2* [ms] to R2* [1/s] values.
"""
r2s_from_t2s(t2s) = 1000 ./ t2s
