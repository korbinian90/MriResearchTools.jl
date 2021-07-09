homodyne(mag, phase; kwargs...) = homodyne(mag .* exp.(1im .* phase); kwargs...)
homodyne(I; kwargs...) = homodyne!(copy(I); kwargs...)
homodyne!(I; dims=1:2, σ=8 .* ones(length(dims))) = I ./= gaussiansmooth3d(I, σ; dims=dims)

function NumART2star(image, TE)
    neco = length(TE)
    t2star(m) = (TE[end]-TE[1]) / 2(neco-1) * (m[1]+m[end]+2sum(m[2:end-1])) / (m[1]-m[end])
    return [t2star(image[I,:]) for I in CartesianIndices(size(image)[1:3])]
end

r2s_from_t2s(t2s) = 1000 ./ t2s
