function readphase(fn; rescale=true, keyargs...)
    phase = niread(fn; keyargs...)
    if phase.header.scl_slope == 0 || rescale
        minp, maxp = Float32.(approxextrema(phase))
        if isapprox(maxp - minp, 2π; atol=0.1) # no rescaling required
            return phase
        end
        minp, maxp = Float32.(approxextrema(phase.raw))
        if isapprox(maxp - minp, 2π; atol=0.1) # no rescaling required, but header wrong
            phase.header.scl_slope = 1
            phase.header.scl_inter = 0
        else # rescaling
            phase.header.scl_slope = 2pi / (maxp - minp)
            phase.header.scl_inter = -pi - minp * phase.header.scl_slope
        end
    end
    return phase
end

function readmag(fn; rescale=false, keyargs...)
    mag = niread(fn; keyargs...)
    if mag.header.scl_slope == 0 || rescale
        mini, maxi = Float32.(approxextrema(mag.raw))
        mag.header.scl_slope = 1 / (maxi - mini)
        mag.header.scl_inter = - mini * mag.header.scl_slope
    end
    return mag
end

Base.copy(x::NIfTI.NIfTI1Header) = NIfTI.NIfTI1Header([getfield(x, k) for k ∈ fieldnames(NIfTI.NIfTI1Header)]...)

function Base.similar(header::NIfTI.NIfTI1Header)
    hdr = copy(header)
    hdr.scl_inter = 0
    hdr.scl_slope = 1
    return hdr
end

header(v::NIVolume) = similar(v.header)

function savenii(image, name, writedir, header=nothing)
    if isnothing(writedir) return end
    if splitext(name)[2] != ".nii"
        name = name * ".nii"
    end
    savenii(image, joinpath(writedir, name); header)
end
"""
    savenii(image, filepath; header=nothing)
save the image at the path
Warning: MRIcro can only open images with types Int32, Int64, Float32, Float64
"""
function savenii(image::AbstractArray, filepath; header=nothing)
    vol = NIVolume([h for h in [header] if h !== nothing]..., image)
    niwrite(filepath, vol)
    return filepath
end

ConvertTypes = Union{BitArray, AbstractArray{UInt8}} #TODO debug NIfTI
MriResearchTools.savenii(image::ConvertTypes, args...;kwargs...) = savenii(Float32.(image), args...;kwargs...)

function write_emptynii(sz, path; datatype=Float64, header=NIVolume(zeros(datatype, 1)).header)
    header = copy(header)
    header.dim = Int16.((length(sz), sz..., ones(8-1-length(sz))...))
    header.datatype = NIfTI.nidatatype(datatype)
    header.bitpix = NIfTI.nibitpix(datatype)

    if isfile(path) rm(path) end
    file = open(path, "w")
    write(file, header)
    close(file)
    GC.gc()
    return niread(path; mmap=true, mode="r+")
end

mmtovoxel(sizemm, nii::NIVolume) = mmtovoxel(sizemm, nii.header)
mmtovoxel(sizemm, header::NIfTI.NIfTI1Header) = mmtovoxel(sizemm, header.pixdim)
mmtovoxel(sizemm, pixdim) = sizemm ./ pixdim

getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)