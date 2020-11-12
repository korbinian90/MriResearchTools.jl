# TODO σ_mm=5


# MCPC-3D-S is implemented for both complex and phase + magnitude data to allow
# working with large mmaped files from disk

# Functions for convenient calling
mcpc3ds(phase::AbstractArray{<:Real}; keyargs...) = angle.(mcpc3ds(exp.(1im .* phase); keyargs...))
mcpc3ds(phase, mag; keyargs...) = mcpc3ds(PhaseMag(phase, mag); keyargs...)
# MCPC3Ds in complex (or PhaseMag)
function mcpc3ds(image; TEs, eco=[1,2], σ=[10,10,5],
        bipolar_correction=false,
        po=zeros(Float64,(size(image)[1:3]..., size(image,5)))
    )
    ΔTE = TEs[eco[2]] - TEs[eco[1]]
    hip = getHIP(image; echoes=eco) # complex
    weight = sqrt.(abs.(hip))
    mask = robustmask(weight)
    phaseevolution = (TEs[eco[1]] / ΔTE) .* romeo(angle.(hip); mag=weight, mask=mask) # different from ASPIRE
    po .= getangle(image, eco[1]) .- phaseevolution
    for icha in 1:size(po, 4)
        po[:,:,:,icha] .= gaussiansmooth3d_phase(po[:,:,:,icha], σ; mask=mask)
    end
    combined = combinewithPO(image, po)
    if bipolar_correction
        G = bipolar_correction!(combined; TEs=TEs, σ=σ, mask=mask)
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

function bipolar_correction!(image; TEs, σ, mask)
    fG = artefact(image, TEs)
    fG .= gaussiansmooth3d_phase(fG, σ; mask=mask)
    romeo!(fG; mag=getmag(image, 1)) # can be replaced by gradient-subtraction-unwrapping
    remove_artefact!(image, fG, TEs)
    return fG
end

getm(TEs) = TEs[1] / (TEs[2] - TEs[1])
getk(TEs) = (TEs[1] + TEs[3]) / TEs[2]

function artefact(I, TEs)
    k = getk(TEs)
    p2 = if abs(k - round(k)) < 0.01 # no unwrapping for integer k required
            getangle(I, 2)
        else
            romeo(getangle(I, 2); mag=getmag(I, 2))
        end
    return  getangle(I, 1) .+ getangle(I, 3) .- k .* p2
end

function remove_artefact!(image, fG, TEs)
    m = getm(TEs)
    k = getk(TEs)
    f = (2 - k) * m - k
    for ieco in 1:size(image, 4)
        t = ifelse(iseven(ieco), m + 1, m) / f
        subtract_angle!(image, ieco, t .* fG)
    end
end

function subtract_angle!(I, echo, sub)
    I[:,:,:,echo] ./= exp.(1im .* sub)
end

## PhaseMag functions
struct PhaseMag
    phase
    mag
end

function combinewithPO(image::PhaseMag, po)
    combined = zeros(Complex{eltype(image)}, size(image)[1:4])
    for icha = 1:size(po, 4)
        @views combined .+= image.mag[:,:,:,:,icha] .* image.mag[:,:,:,:,icha] .* exp.(1im .* (image.phase[:,:,:,:,icha] .- po[:,:,:,icha]))
    end
    return PhaseMag(angle.(combined), sqrt.(abs.(combined)))
end

function subtract_angle!(I::PhaseMag, echo, sub)
    I.phase[:,:,:,echo] .-= sub
    I.phase .= rem2pi.(I.phase, RoundNearest)
end

## Utility
Base.iterate(t::PhaseMag, i) = if i == 1 (t.phase, 2) else (t.mag, 3) end
Base.iterate(t::PhaseMag) = t.phase, t.mag
Base.eltype(t::PhaseMag) = promote_type(eltype(t.mag), eltype(t.phase))
Base.size(t::PhaseMag, args...) = size(t.phase, args...)
getHIP(data::PhaseMag; keyargs...) = getHIP(data.mag, data.phase; keyargs...)
getangle(c, echo) = angle.(ecoview(c, echo))
getangle(d::PhaseMag, echo) = ecoview(d.phase, echo)
getmag(c, echo) = abs.(ecoview(c, echo))
getmag(d::PhaseMag, echo) = ecoview(d.mag, echo)
ecoview(a, echo) = dimview(a, 4, echo)
dimview(a, dim, i) = view(a, ntuple(x -> if x == dim i else (:) end, ndims(a))...)
