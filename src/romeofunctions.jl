const romeo = unwrap # access unwrap function via alias romeo
const romeo! = unwrap!

"""
    calculateB0_unwrapped(unwrapped_phase, mag, TEs)

Calculates B0 in [Hz] from unwrapped phase.
TEs in [ms].
The phase offsets have to be removed prior.

See also [`mcpc3ds`](@ref)
"""
function calculateB0_unwrapped(unwrapped_phase, mag, TEs)
    dims = 4
    TEs = to_dim(TEs, 4)
    B0 = (1000 / 2π) * sum(unwrapped_phase .* mag .* mag .* TEs; dims) ./ sum(mag .* mag .* TEs.^2; dims) |> I -> dropdims(I; dims)
    B0[.!isfinite.(B0)] .= 0
    return B0
end

romeovoxelquality = voxelquality
