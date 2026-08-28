"""
    readphase(filename; rescale=true, fix_ge=false, keyargs...)

Reads the NIfTI phase with sanity checking and optional rescaling to [-π;π].
Warning for GE: `fix_ge=true` is probably required and will add pi to every second slice.

# Examples
```julia-repl
julia> phase = readphase("phase.nii")
```

### Optional keyargs are forwarded to `niread`:
$(@doc niread)
"""
function readphase(filename; rescale=true, fix_ge=false, keyargs...)
    phase = niread(filename; keyargs...)
    if phase.header.scl_slope == 0 # slope of 0 is always wrong
        phase.header.scl_slope = 1
    end
    if rescale
        minp, maxp = Float32.(approxextrema(phase))
        if !isapprox(maxp - minp, 2π; atol=0.1) # rescaling required
            minp, maxp = Float32.(approxextrema(phase.raw))
            if isapprox(maxp - minp, 2π; atol=0.1) # no rescaling required, but header wrong
                phase.header.scl_slope = 1
                phase.header.scl_inter = 0
            else # rescaling
                phase.header.scl_slope = 2pi / (maxp - minp)
                phase.header.scl_inter = -pi - minp * phase.header.scl_slope
            end
        end
    end
    if fix_ge
        fix_ge_phase!(phase.raw)
        phase.header.scl_slope = -phase.header.scl_slope # phase is inverted
    end
    return phase
end

# Add pi to every second slice
function fix_ge_phase!(phase::AbstractArray{T}) where T
    minp, maxp = approxextrema(phase)
    if T <: Integer
        pi = round(T, (maxp - minp) / 2)
    else
        pi = (maxp - minp) / 2
    end
    every_second_slice = selectdim(phase, 3, 2:2:size(phase, 3))
    every_second_slice .= rem.(every_second_slice .+ pi, 2pi, RoundNearest)
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

function savenii(image, name, writedir, header=nothing; kwargs...)
    if isnothing(writedir) return end
    if !(last(splitext(name)) in [".nii", ".gz"])
        name = "$name.nii"
    end
    savenii(image, joinpath(writedir, name); header, kwargs...)
end
"""
    savenii(image::AbstractArray, filepath; header=nothing, datatype=Float32, kwargs...)

    savenii(image::AbstractArray, name, writedir, header=nothing, kwargs...)

The image is written as `datatype`, which defaults to `Float32` (`ComplexF32`
for complex data) whatever the array holds - a Float64 result is half the file
for no less information, and an algorithm that needs Float64 internally should
narrow here rather than compute in Float32. Pass `datatype=nothing` to write the
array's own element type unchanged, or any type to force it.

Note that integer element types are promoted to `Float32` by that default, on
purpose: NIfTI.jl 0.6.2 writes a wrong `bitpix` for every integer width other
than 32 bits, so a `UInt8` mask would get a header contradicting itself. Pass
`datatype=UInt8` to write it anyway.

Warning: MRIcro can only open images with types Int32, Int64, Float32, Float64

# Examples
```julia-repl
julia> savenii(ones(64,64,5), "image.nii")

julia> savenii(ones(64,64,5), "image2", "folder")

julia> savenii(ones(64,64,5), "image2", "folder"; voxel_size=(0.54,0.54,2.0))
```
"""
function savenii(image::AbstractArray, filepath; header=nothing, datatype=default_output_type(eltype(image)), kwargs...)
    image = to_output_type(image, datatype)
    vol = NIVolume([h for h in [header] if h !== nothing]..., image; kwargs...)
    dir = dirname(filepath)
    if !isdir(dir)
        mkpath(dir)
    end
    niwrite(filepath, vol)
    return filepath
end

# The on-disk NIfTI datatype comes from the eltype of the array handed to
# niwrite and from nothing else - a header passed as `header=` has its datatype
# and bitpix discarded by NIfTI.jl's niupdate - so converting the array is the
# only way to control it. This is the one function every package in the family
# writes through, which is why the rule lives here.
#
# Integers are promoted to Float32 as a workaround, not a preference: NIfTI.jl
# 0.6.2 derives bitpix from typeof(one(T)*1.0f0+1.0f0) rather than the element
# type, so it writes a self-contradicting header for every integer width other
# than 32 bits (UInt8 -> datatype 2 with bitpix 32, Int16 -> 4 with 32, Int64 ->
# 1024 with 32). A UInt8 mask, one byte per voxel instead of four, is exactly
# what cannot be written correctly today. Pass datatype=UInt8 to do it anyway.
"""
    default_output_type(T)

The element type `savenii` writes for an array of element type `T`: `Float32`,
or `ComplexF32` for complex data. See [`savenii`](@ref) to override it.
"""
default_output_type(::Type{<:Complex}) = ComplexF32
default_output_type(::Type) = Float32

to_output_type(image, ::Nothing) = image
to_output_type(image, ::Type{T}) where {T} = eltype(image) === T ? image : T.(image)

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
