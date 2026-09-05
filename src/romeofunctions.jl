const romeo = unwrap # access unwrap function via alias romeo
const romeo! = unwrap!

"""
    calculateB0_unwrapped(unwrapped_phase, mag, TEs)

Calculates B0 in [Hz] from unwrapped phase.
TEs in [ms].
The phase offsets have to be removed prior.

See also [`mcpc3ds`](@ref)
"""
calculateB0_unwrapped(unwrapped_phase, mag, TEs, type::Symbol=:phase_snr) = calculateB0_unwrapped(unwrapped_phase, mag, TEs, Val(type))
# with the weighting as a Val every array type is static, which compiled programs need
function calculateB0_unwrapped(unwrapped_phase, mag, TEs, ::Val{type}) where type
    TEs = to_dim(TEs, Val(4))
    weight = get_B0_phase_weighting(mag, TEs, Val(type))
    B0 = (1000 / 2π) * sum(unwrapped_phase ./ TEs .* weight; dims=4) ./ sum(weight; dims=4)
    B0 = reshape(B0, size(B0, 1), size(B0, 2), size(B0, 3))
    B0[.!isfinite.(B0)] .= 0
    return B0
end

get_B0_phase_weighting(mag, TEs, type::Symbol) = get_B0_phase_weighting(mag, TEs, Val(type))
function get_B0_phase_weighting(mag, TEs, ::Val{type}) where type
    if type == :phase_snr
        mag .* TEs
    elseif type == :phase_var
        mag .* mag .* TEs .* TEs
    elseif type == :average
        to_dim(ones(length(TEs)), Val(4))
    elseif type == :TEs
        TEs
    elseif type == :mag
        mag
    elseif type == :simulated_mag
        mag = to_dim(exp.(-TEs / 20), Val(4))
        mag .* TEs
    else
        error("The phase weighting option '$type' is not defined!")
    end
end

get_B0_snr(mag, TEs, type::Symbol=:phase_snr) = get_B0_snr(mag, TEs, Val(type))
function get_B0_snr(mag, TEs, ::Val{type}) where type
    weight = get_B0_phase_weighting(mag, to_dim(TEs, Val(4)), Val(type))
    sum(mag .* weight; dims=4) ./ sum(weight; dims=4)
end

const romeovoxelquality = voxelquality
