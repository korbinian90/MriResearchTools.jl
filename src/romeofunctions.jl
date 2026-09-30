const romeo = unwrap # access unwrap function via alias romeo
const romeo! = unwrap!

"""
    calculateB0_unwrapped(unwrapped_phase, mag, TEs)

Calculates B0 in [Hz] from unwrapped phase.
TEs in [ms].
The phase offsets have to be removed prior.

See also [`mcpc3ds`](@ref)
"""
calculateB0_unwrapped(unwrapped_phase, mag, TEs) = calculateB0_unwrapped(unwrapped_phase, mag, TEs, Val(:phase_snr))
calculateB0_unwrapped(unwrapped_phase, mag, TEs, type::Symbol) = calculateB0_unwrapped(unwrapped_phase, mag, TEs, Val(type))
# with the weighting as a Val every array type is static, which compiled programs need
function calculateB0_unwrapped(unwrapped_phase, mag, TEs, ::Val{type}) where type
    TEs = to_dim(TEs, Val(4))
    weight = get_B0_phase_weighting(mag, TEs, Val(type))
    B0 = (1000 / 2π) * sum(unwrapped_phase ./ TEs .* weight; dims=4) ./ sum(weight; dims=4)
    B0 = _drop_echo_dim(B0)
    B0[.!isfinite.(B0)] .= 0
    return B0
end

# drops the summed echo dimension 4, keeping a channel dimension 5
_drop_echo_dim(B0::AbstractArray{<:Any,3}) = B0
_drop_echo_dim(B0::AbstractArray{<:Any,N}) where N = reshape(B0, ntuple(i -> size(B0, i < 4 ? i : i + 1), Val(N - 1)))

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
