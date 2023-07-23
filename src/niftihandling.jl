"""
    readphase(filename; rescale=true, keyargs...)

Reads the NIfTI phase with sanity checking and optional rescaling to [-π;π].

# Examples
```julia-repl
julia> phase = readphase("phase.nii")
```

### Optional keyargs are forwarded to `niread`:
$(@doc niread)
"""
function readphase(filename; rescale=true, keyargs...)
    phase = niread(filename; keyargs...)
    if phase.header.scl_slope == 0 # slope of 0 is always wrong
        phase.header.scl_slope = 1
    end
    if rescale
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

"""
    readmag(filename; rescale=false, keyargs...)

Reads the NIfTI magnitude with sanity checking and optional rescaling to [0;1].

# Examples
```julia-repl
julia> mag = readmag("mag.nii")
```

### Optional keyargs are forwarded to `niread`:
$(@doc niread)
"""
function readmag(fn; rescale=false, keyargs...)
    mag = niread(fn; keyargs...)
    if mag.header.scl_slope == 0
        mag.header.scl_slope = 1
    end
    if rescale
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

"""
    header(v::NIVolume)

Returns a copy of the header with the orientation information.

# Examples
```julia-repl
julia> vol = readmag("image.nii")
julia> hdr = header(vol)
julia> savenii(vol .+ 10, "vol10.nii"; header=hdr)
```
"""
header(v::NIVolume) = similar(v.header)

function savenii(image, name, writedir, header=nothing)
    if isnothing(writedir) return end
    if !(last(splitext(name)) in [".nii", ".gz"])
        name = "$name.nii"
    end
    savenii(image, joinpath(writedir, name); header)
end
"""
    savenii(image::AbstractArray, filepath; header=nothing)

    savenii(image::AbstractArray, name, writedir, header=nothing)

Warning: MRIcro can only open images with types Int32, Int64, Float32, Float64

# Examples
```julia-repl
julia> savenii(ones(64,64,5), "image.nii")

julia> savenii(ones(64,64,5), "image2", "folder")
```
"""
function savenii(image::AbstractArray, filepath; header=nothing)
    vol = NIVolume([h for h in [header] if h !== nothing]..., image)
    niwrite(filepath, vol)
    return filepath
end

ConvertTypes = Union{BitArray, AbstractArray{UInt8}} #TODO debug NIfTI
MriResearchTools.savenii(image::ConvertTypes, args...;kwargs...) = savenii(Float32.(image), args...;kwargs...)

"""
    write_emptynii(size, path; datatype=Float32, header=NIVolume(zeros(datatype, 1)).header)

Writes an empty NIfTI image to disk that can be used for memory-mapped access.

# Examples
```julia-repl
julia> vol = write_emptynii((64,64,64), "empty.nii")
julia> vol.raw[:,:,1] .= ones(64,64) # synchronizes mmapped file on disk
```

Warning: MRIcro can only open images with types Int32, Int64, Float32, Float64
"""
function write_emptynii(sz, path; datatype=Float32, header=NIVolume(zeros(datatype, 1)).header)
    header = copy(header)
    header.dim = Int16.((length(sz), sz..., ones(8-1-length(sz))...))
    header.datatype = NIfTI.eltype_to_int16(datatype)
    header.bitpix = NIfTI.nibitpix(datatype)

    if isfile(path) rm(path) end
    open(path, "w") do file
        write(file, header)
        write(file, Int32(0)) # offset of 4 bytes
    end
    return niread(path; mmap=true, mode="r+")
end

mmtovoxel(sizemm, nii::NIVolume) = mmtovoxel(sizemm, nii.header)
mmtovoxel(sizemm, header::NIfTI.NIfTI1Header) = mmtovoxel(sizemm, header.pixdim)
mmtovoxel(sizemm, pixdim) = sizemm ./ pixdim

getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)

function Base.setindex!(vol::NIVolume{<:AbstractFloat}, v, i...)
    scaled = v / vol.header.scl_slope + vol.header.scl_inter
    setindex!(vol.raw, scaled, i...)
end

function affine(nii::NIVolume)
    fields = [:srow_x, :srow_y, :srow_z]
    rows = (transpose(collect(getfield(nii.header, f))) for f in fields)
    return vcat(rows..., [0 0 0 1])
end
