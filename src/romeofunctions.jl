const romeo = unwrap # access unwrap function via alias romeo
const romeo! = unwrap!

"""
    calculateB0_unwrapped(unwrapped_phase, mag, TEs)

Calculates B0 in [Hz] from unwrapped phase.
TEs in [ms].
The phase offsets have to be removed prior.

See also [`mcpc3ds`](@ref)
"""
function calculateB0_unwrapped(unwrapped_phase, mag, TEs, type=:phase_snr)
    dims = 4
    TEs = to_dim(TEs, 4)
    weight = get_B0_phase_weighting(mag, TEs, type)
    B0 = (1000 / 2π) * sum(unwrapped_phase ./ TEs .* weight; dims) ./ sum(weight; dims) |> I -> dropdims(I; dims)
    B0[.!isfinite.(B0)] .= 0
    return B0
end

function get_B0_phase_weighting(mag, TEs, type)
    if type == :phase_snr
        mag .* TEs
    elseif type == :average
        to_dim(ones(length(TEs)), 4)
    elseif type == :TEs
        TEs
    elseif type == :mag
        mag
    elseif type == :simulated_mag
        mag = to_dim(exp.(-TEs / 20), 4)
        mag .* TEs
    else
        error("The phase weighting option '$type' is not defined!")
    end
end

function get_B0_snr(mag, TEs, type=:phase_snr)
    weight = get_B0_phase_weighting(mag, to_dim(TEs, 4), type)
    sum(mag .* weight; dims=4) ./ sum(weight; dims=4)
end

romeovoxelquality = voxelquality
