# TODO sigma_mm=5


# MCPC-3D-S is implemented for both complex and phase + magnitude data to allow

"""
    mcpc3ds(phase, mag; TEs, keyargs...)

    mcpc3ds(compl; TEs, keyargs...)

    mcpc3ds(phase; TEs, keyargs...)

Perform MCPC-3D-S coil combination and phase offset removal on 4D (multi-echo) and 5D (multi-echo, uncombined) input.

## Optional Keyword Arguments
- `echoes`: only use the defined echoes. default: `echoes=[1,2]`
- `sigma`: smoothing parameter for phase offsets. default: `sigma=[10,10,5]`
- `bipolar_correction`: removes linear phase artefact. default: `bipolar_correction=false`
- `po`: phase offsets are stored in this array. Can be used to retrieve phase offsets or work with memory mapping.

# Examples
```julia-repl
julia> phase = readphase("phase5D.nii")
julia> mag = readmag("mag5D.nii")
julia> combined = mcpc3ds(phase, mag; TEs=[4,8,12])
```
For very large files that don't fit into memory, the uncombined data can be processed with memory mapped to disk:
```julia-repl
julia> phase = readphase("phase5D.nii"; mmap=true)
julia> mag = readmag("mag5D.nii"; mmap=true)
julia> po_size = (size(phase)[1:3]..., size(phase,5))
julia> po = write_emptynii(po_size, "po.nii")
julia> combined = mcpc3ds(phase, mag; TEs=[4,8,12], po)
```
"""
mcpc3ds(phase::AbstractArray{<:Real}; keyargs...) = angle.(mcpc3ds(exp.(1im .* phase); keyargs...))
mcpc3ds(phase, mag; keyargs...) = mcpc3ds(PhaseMag(phase, mag); keyargs...)
# MCPC3Ds in complex (or PhaseMag)
function mcpc3ds(image; TEs, echoes=[1,2], sigma=[10,10,5],
        bipolar_correction=false,
        po=zeros(getdatatype(image),(size(image)[1:3]..., size(image,5)))
    )
    ΔTE = TEs[echoes[2]] - TEs[echoes[1]]
    hip = getHIP(image; echoes) # complex
    weight = sqrt.(abs.(hip))
    mask = robustmask(weight)
    # TODO try to include additional second-phase information in the case of 3+ echoes for ROMEO, maybe phase2=phase[3]-phase[2], TEs=[dTE21, dTE32]
    phaseevolution = (TEs[echoes[1]] / ΔTE) .* romeo(angle.(hip); mag=weight, mask) # different from ASPIRE
    po .= getangle(image, echoes[1]) .- phaseevolution
    for icha in axes(po, 4)
        po[:,:,:,icha] .= gaussiansmooth3d_phase(view(po,:,:,:,icha), sigma; mask)
    end
    combined = combinewithPO(image, po)
    if bipolar_correction
        fG = bipolar_correction!(combined; TEs, sigma, mask)
    end
    return combined
end

"""
    mcpc3ds_meepi(phase, mag; TEs, keyargs...)

    mcpc3ds_meepi(compl; TEs, keyargs...)

    mcpc3ds_meepi(phase; TEs, keyargs...)

Perform MCPC-3D-S phase offset removal on 5D MEEPI (multi-echo, multi-timepoint) input.
The phase offsets are calculated for the template timepoint and removed from all volumes.

## Optional Keyword Arguments
- `echoes`: only use the defined echoes. default: `echoes=[1,2]`
- `sigma`: smoothing parameter for phase offsets. default: `sigma=[10,10,5]`
- `po`: phase offsets are stored in this array. Can be used to retrieve phase offsets or work with memory mapping.
- `template_tp`: timepoint for the template calculation. default: `template_tp=1`
"""
mcpc3ds_meepi(phase::AbstractArray{<:Real}; keyargs...) = mcpc3ds_meepi(exp.(1im .* phase); keyargs...)
mcpc3ds_meepi(phase, mag; keyargs...) = mcpc3ds_meepi(PhaseMag(phase, mag); keyargs...)
function mcpc3ds_meepi(image; template_tp=1, po=zeros(getdatatype(image),size(image)[1:3]), kwargs...)
    template = selectdim(image, 5, template_tp)
    mcpc3ds(template; po, bipolar_correction=false, kwargs...) # calculates and sets po

    corrected_phase = similar(po, size(image))
    for tp in axes(image, 5)
        corrected_phase[:,:,:,:,tp] = getangle(combinewithPO(selectdim(image, 5, tp), po))
    end
    return corrected_phase
end

function combinewithPO(compl, po)
    combined = zeros(eltype(compl), size(compl)[1:4])
    @sync for iecho in axes(combined, 4)
        Threads.@spawn @views combined[:,:,:,iecho] = sum(abs.(compl[:,:,:,iecho,icha]) .* compl[:,:,:,iecho,icha] ./ exp.(1im .* po[:,:,:,icha]) for icha in axes(po,4))
    end
    return combined ./ sqrt.(abs.(combined))
end

## Bipolar correction
# see https://doi.org/10.34726/hss.2021.43447, page 53, 3.1.3 Bipolar Corrections
function bipolar_correction!(image; TEs, sigma, mask)
    fG = artefact(image, TEs)
    fG .= gaussiansmooth3d_phase(fG, sigma; mask)
    romeo!(fG; mag=getmag(image, 1), correctglobal=true) # can be replaced by gradient-subtraction-unwrapping
    remove_artefact!(image, fG, TEs)
    return fG
end

getm(TEs) = TEs[1] / (TEs[2] - TEs[1])
getk(TEs) = (TEs[1] + TEs[3]) / TEs[2]

function artefact(I, TEs)
    k = getk(TEs)
    ϕ1 = getangle(I, 1)
    ϕ2 = getangle(I, 2)
    ϕ3 = getangle(I, 3)
    if abs(k - round(k)) < 0.01 # no unwrapping for integer k required
        ϕ2 = romeo(ϕ2; mag=getmag(I, 2))
    end
    return  ϕ1 .+ ϕ3 .- k .* ϕ2
end

function remove_artefact!(image, fG, TEs)
    m = getm(TEs)
    k = getk(TEs)
    f = (2 - k) * m - k
    for ieco in axes(image, 4)
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
    for icha in axes(po, 4)
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
Base.axes(t::PhaseMag, dim) = 1:size(t.phase, dim)
getHIP(data::PhaseMag; keyargs...) = getHIP(data.mag, data.phase; keyargs...)
getangle(c, echo=:) = angle.(ecoview(c, echo))
getangle(d::PhaseMag, echo=:) = ecoview(d.phase, echo)
getmag(c, echo=:) = abs.(ecoview(c, echo))
getmag(d::PhaseMag, echo=:) = ecoview(d.mag, echo)
Base.selectdim(A::PhaseMag, d, i) = PhaseMag(selectdim(A.mag, d, i), selectdim(A.phase, d, i))
ecoview(a, echo) = dimview(a, 4, echo)
dimview(a, dim, i) = view(a, ntuple(x -> if x == dim i else (:) end, ndims(a))...)
getdatatype(cx::AbstractArray{<:Complex{T}}) where T = T
getdatatype(other) = eltype(other)
