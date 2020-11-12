homodyne(mag, phase; kwargs...) = homodyne(mag .* exp.(1im .* phase); kwargs...)
homodyne(I; kwargs...) = homodyne!(copy(I); kwargs...)
homodyne!(I; dims=1:2, σ=8 .* ones(length(dims))) = I ./= gaussiansmooth3d(I, σ; dims=dims)

# calculates B0 in [Hz]
function calculateB0_unwrapped(unwrapped_phase, mag, TEs)
    dims=4
    return (1000 / 2π) * sum(unwrapped_phase .* mag; dims=dims) ./ sum(mag .* to_dim(TEs,4); dims=dims) |> I -> dropdims(I; dims=dims)
end
