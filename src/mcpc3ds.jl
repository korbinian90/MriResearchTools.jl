# TODO σ_mm=5, bipolar_correction
struct PhaseMag
    phase
    mag
end
Base.iterate(t::PhaseMag, i) = if i == 1 (t.phase, 2) else (t.mag, 3) end
Base.eltype(t::PhaseMag) = promote_type(eltype(t.mag), eltype(t.phase))
Base.size(t::PhaseMag, args...) = size(t.phase, args...)
MrData = Union{PhaseMag, AbstractArray{<:Complex}}
getHIP(data::PhaseMag; keyargs...) = getHIP(data.mag, data.phase; keyargs...)
getangle(c, echo) = angle.(view(c,:,:,:,echo,:))
getangle(d::PhaseMag, echo) = view(d.phase,:,:,:,echo,:)
getmag(c, echo) = abs.(view(c,:,:,:,echo,:))
getmag(d::PhaseMag, echo) = view(d.mag,:,:,:,echo,:)

# MCPC-3D-S is implemented for both complex and phase, magnitude data to allow
# working with large mmaped files from disk

# Functions for convenient calling
mcpc3ds(phase; keyargs...) = angle.(mcpc3ds(exp.(1im .* phase); keyargs...))
mcpc3ds(phase, mag; keyargs...) = mcpc3ds(PhaseMag(phase, mag); keyargs...)
# MCPC3Ds in complex
function mcpc3ds(image::MrData; TEs, eco=[1,2], σ=[10,10,5],
        bipolar_correction=false,
        po=zeros(Float64,(size(image)[1:3]..., size(image,5)))
    )
    ΔTE = TEs[eco[2]] - TEs[eco[1]]
    @time hip = getHIP(image; echoes=eco) # complex
    @time weight = sqrt.(abs.(hip))
    @time mask = robustmask(weight)
    @time phaseevolution = (TEs[eco[1]] / ΔTE) .* romeo(angle.(hip); mag=weight, mask=mask) # different from ASPIRE
    @time po .= getangle(image, eco[1]) .- phaseevolution
    @time for icha in 1:size(po, 4)
        po[:,:,:,icha] .= angle.(gaussiansmooth3d(exp.(1im .* po[:,:,:,icha]), σ; mask=mask))
    end
    @time combined = combinewithPO(image, po)
    if bipolar_correction
        bipolar_correction!(combined; TEs=TEs, σ=σ, mask=mask)
    end
    return combined
end

function combinewithPO(compl, po)
    combined = zeros(eltype(compl), size(compl)[1:4])
    for icha = 1:size(po, 4)
        @views combined .+= abs.(compl[:,:,:,:,icha]) .* compl[:,:,:,:,icha] ./ exp.(1im .* po[:,:,:,icha])
    end
    return combined ./ sqrt.(abs.(combined))
end

function combinewithPO(image::PhaseMag, po)
    combined = zeros(Complex{eltype(image)}, size(image)[1:4])
    for icha = 1:size(po, 4)
        @views combined .+= image.mag[:,:,:,:,icha] .* image.mag[:,:,:,:,icha] .* exp.(1im .* (image.phase[:,:,:,:,icha] .- po[:,:,:,icha]))
    end
    return PhaseMag(angle.(combined), sqrt.(abs.(combined)))
end

function bipolar_correction!(image; TEs, σ, mask, unwrapping=TEs[3] % TEs[1] > 0.1)
    G = artefact(image, TEs)
    artefact = angle.(gaussiansmooth3d!(exp.(1im .* artefact)), σ; mask=mask)
    if unwrapping
        romeo!(G; mag=getmag(image, 1))
    end
    return remove_artefact!(image, G, TEs)
end

artefact(I, TEs) = (TEs[3] / TEs[1]) .* getangle(I, 1) .- getangle(I, 3)

function remove_artefact!(image, G, TEs)
    m = TEs[1] / (TEs[2] - TEs[1])
    for ieco in 1:size(image, 4)
        f = if iseven(ieco) (m + 1) else m end
        subtract_angle!(image, echo, f * G)
    end
    return image
end

function subtract_angle!(I::PhaseMag, echo, sub)
    I.phase[:,:,:,echo] .-= sub
end
function subtract_angle!(I, echo, sub)
    I[:,:,:,echo] ./= exp.(1im .* sub)
end
