# TODO σ_mm=5, bipolar_correction

# Functions for convenient calling
mcpc3ds(phase; keyargs...) = angle.(mcpc3ds(exp.(1im .* phase); keyargs...))
mcpc3ds(phase, mag; keyargs...) = mcpc3ds(mag .* exp.(1im .* phase); keyargs...) |> I -> (angle.(I), abs.(I))
# MCPC3Ds in complex
function mcpc3ds(compl::AbstractArray{<:Complex}; TEs, eco=[1,2], σ=[10,10,5],
        po=zeros(eltype(compl),(size(compl)[1:3]..., size(compl,5)))
    )
    ΔTE = TEs[eco[2]] - TEs[eco[1]]
    hip = getHIP(compl; echoes=eco) # complex
    weight = sqrt.(abs.(hip))
    mask = robustmask(weight)
    phaseevolution = (TEs[eco[1]] / ΔTE) .* romeo(angle.(hip); mag=weight, mask=mask) # different from ASPIRE
    po .= compl[:,:,:,eco[1],:] ./ exp.(1im .* phaseevolution)
    gaussiansmooth3d!(po, σ; mask=mask) # weighted by mag1 and masked
    po ./= abs.(po)
    return combinewithPO(compl, po)
end

function combinewithPO(compl, po)
    combined = zeros(eltype(compl), size(compl)[1:4])
    for iCha = 1:size(po, 5)
        combined .+= abs.(compl[:,:,:,:,iCha]) .* compl[:,:,:,:,iCha] ./ po[:,:,:,iCha]
    end
    return combined ./ sqrt.(abs.(combined))
end
