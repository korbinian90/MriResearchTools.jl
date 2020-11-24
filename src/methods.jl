homodyne(mag, phase; kwargs...) = homodyne(mag .* exp.(1im .* phase); kwargs...)
homodyne(I; kwargs...) = homodyne!(copy(I); kwargs...)
homodyne!(I; dims=1:2, σ=8 .* ones(length(dims))) = I ./= gaussiansmooth3d(I, σ; dims=dims)

# calculates B0 in [Hz]
function calculateB0_unwrapped(unwrapped_phase, mag, TEs)
    dims = 4
    TEs = to_dim(TEs, 4)
    return (1000 / 2π) * sum(unwrapped_phase .* mag .* mag .* TEs; dims=dims) ./ sum(mag .* mag .* TEs.^2; dims=dims) |> I -> dropdims(I; dims=dims)
end
