# TODO σ_mm=5, bipolar_correction

# memory efficient inplace only for complex nifti
mcpc3ds(mag, phase; keyargs...) = mcpc3ds(mag .* exp.(1im .* phase); keyargs...)
function mcpc3ds(compl; TEs, eco=[1,2], σ=[10,10,5],
        po=zeros(eltype(compl),(size(compl)[1:3]..., size(compl,5)))
    )
    ΔTE = TEs[eco[2]] - TEs[eco[1]]
    hip = getHIP(compl; echoes=eco) # complex
    weight = sqrt.(abs.(hip))
    mask = robustmask(weight)
    @time phaseevolution = (TEs[eco[1]] / ΔTE) .* romeo(angle.(hip); mag=weight, mask=mask) # only difference to ASPIRE
    @time po .= compl[:,:,:,eco[1],:] ./ exp.(1im .* phaseevolution)
    @time gaussiansmooth3d!(po, σ; mask=mask) # weighted by mag1 and masked
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
