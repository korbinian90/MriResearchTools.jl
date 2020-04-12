homodyne(mag, phase; kwargs...) = homodyne(mag .* exp.(1im .* phase); kwargs...)
homodyne(I; kwargs...) = homodyne!(copy(I); kwargs...)
homodyne!(I; σ=8, dims=1:2) = I ./= gaussiansmooth3d(I, σ; dims=dims)
