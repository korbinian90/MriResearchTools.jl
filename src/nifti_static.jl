# NIfTI reading and writing with statically known types.
#
# NIfTI.jl decides the element type of a volume from the file header at run time
# and serialises headers by iterating over field names. Neither can be resolved
# statically, so a program that reads through `niread` cannot be compiled with
# `juliac --trim`. The readers here decide the element type in one explicit chain
# over the NIfTI datatypes and return Float32 arrays with five dimensions; the
# writer packs the header field by field. The values match what `readphase`,
# `readmag` and `savenii` produce through NIfTI.jl, and written files are
# byte-identical.

const HEADER1_SIZE = Int(NIfTI.SIZEOF_HDR1)

# offsets of the packed on-disk layout, which has no alignment padding
const HEADER1_OFFSETS = let offsets = Int[], offset = 0
    for T in fieldtypes(NIfTI.NIfTI1Header)
        push!(offsets, offset)
        offset += sizeof(T)
    end
    @assert offset == HEADER1_SIZE
    Tuple(offsets)
end

@generated function _copy_header(x::NIfTI.NIfTI1Header)
    H = NIfTI.NIfTI1Header
    return Expr(:new, H, (:(getfield(x, $i)) for i in 1:fieldcount(H))...)
end
Base.copy(x::NIfTI.NIfTI1Header) = _copy_header(x)

_swap(x::Union{Int8,UInt8}, swapped) = x
_swap(x::Number, swapped) = swapped ? ntoh(x) : x
_swap(x::Tuple, swapped) = map(y -> _swap(y, swapped), x)

@generated function _unpack_header(bytes::Vector{UInt8}, swapped::Bool)
    H = NIfTI.NIfTI1Header
    fields = [:(_swap(unsafe_load(Ptr{$(fieldtype(H, i))}(pointer(bytes) + $(HEADER1_OFFSETS[i]))), swapped))
              for i in 1:fieldcount(H)]
    return :(GC.@preserve bytes $(Expr(:new, H, fields...)))
end

@generated function _pack_header(hdr::NIfTI.NIfTI1Header)
    H = NIfTI.NIfTI1Header
    stores = [:(unsafe_store!(Ptr{$(fieldtype(H, i))}(pointer(bytes) + $(HEADER1_OFFSETS[i])), getfield(hdr, $i)))
              for i in 1:fieldcount(H)]
    return quote
        bytes = zeros(UInt8, HEADER1_SIZE)
        GC.@preserve bytes begin
            $(stores...)
        end
        bytes
    end
end

# Returns the header and whether the file is big-endian
function _read_header(io)
    bytes = read(io, HEADER1_SIZE)
    length(bytes) == HEADER1_SIZE || throw(ArgumentError("file is shorter than a NIfTI-1 header"))
    sizeof_hdr = GC.@preserve bytes unsafe_load(Ptr{Int32}(pointer(bytes)))
    if sizeof_hdr == NIfTI.SIZEOF_HDR1
        return _unpack_header(bytes, false), false
    elseif bswap(sizeof_hdr) == NIfTI.SIZEOF_HDR1
        return _unpack_header(bytes, true), true
    elseif sizeof_hdr == NIfTI.SIZEOF_HDR2 || bswap(sizeof_hdr) == NIfTI.SIZEOF_HDR2
        throw(ArgumentError("NIfTI-2 files are not supported"))
    end
    throw(ArgumentError("not a NIfTI file"))
end

_open(filename) = NIfTI.niopen(open(filename, "r"))

"""
    loadheader(filename)

The NIfTI-1 header of `filename`, read without the data. See [`loadnii`](@ref).
"""
function loadheader(filename)
    io = _open(filename)
    try
        return first(_read_header(io))
    finally
        close(io)
    end
end

# The .img file that belongs to a .hdr of dual file storage
function _img_path(filename)
    base = _strip_extension(endswith(filename, ".gz") ? filename[1:end-3] : filename)
    isfile(base * ".img") && return base * ".img"
    isfile(base * ".img.gz") && return base * ".img.gz"
    throw(ArgumentError("NIfTI file is dual file storage, but $base.img does not exist"))
end

function _strip_extension(name)
    i = lastindex(name)
    while i >= firstindex(name)
        name[i] == '.' && return name[1:prevind(name, i)]
        i = prevind(name, i)
    end
    return name
end

# Five dimensions, with trailing ones for lower rank
function _dims5(dim::NTuple{8,Int16})
    nd = Int(dim[1])
    1 <= nd <= 5 || throw(ArgumentError("$(nd)D NIfTI data is not supported"))
    return ntuple(i -> i <= nd ? Int(dim[i+1]) : 1, Val(5))
end

function _read_raw(io, ::Type{T}, dims::NTuple{5,Int}, swapped::Bool) where T
    raw = read!(io, Array{T,5}(undef, dims))
    if swapped
        raw .= ntoh.(raw)
    end
    return raw
end

# One statically resolvable call of `f(raw, header)` per NIfTI datatype
function _read_typed(f::F, io, hdr::NIfTI.NIfTI1Header, dims, swapped) where F
    dt = hdr.datatype
    dt == NIfTI.eltype_to_int16(Float32) && return f(_read_raw(io, Float32, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(Int16) && return f(_read_raw(io, Int16, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(Float64) && return f(_read_raw(io, Float64, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(UInt8) && return f(_read_raw(io, UInt8, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(Int32) && return f(_read_raw(io, Int32, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(UInt16) && return f(_read_raw(io, UInt16, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(Int8) && return f(_read_raw(io, Int8, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(UInt32) && return f(_read_raw(io, UInt32, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(Int64) && return f(_read_raw(io, Int64, dims, swapped), hdr)
    dt == NIfTI.eltype_to_int16(UInt64) && return f(_read_raw(io, UInt64, dims, swapped), hdr)
    throw(ArgumentError("NIfTI datatype $dt is not supported"))
end

# Calls `f(raw, header)` with the raw array in the element type of the file and
# returns its result together with the header. The header scaling is cleared
# afterwards: `f` applies it to the data, and a header that still carried it
# would scale the data a second time when written out with it.
function _read_nii(f::F, filename) where F
    io = _open(filename)
    try
        hdr, swapped = _read_header(io)
        dims = _dims5(hdr.dim)
        if hdr.magic == NIfTI.NP1_MAGIC
            read(io, Int(hdr.vox_offset) - HEADER1_SIZE) # extensions are not used
        else
            close(io)
            io = _open(_img_path(filename))
        end
        data = _read_typed(f, io, hdr, dims, swapped)
        hdr.scl_slope = 1
        hdr.scl_inter = 0
        return data, hdr
    finally
        close(io)
    end
end

function _slope_inter(hdr)
    slope = hdr.scl_slope == 0 ? 1f0 : hdr.scl_slope # slope of 0 is always wrong
    return slope, hdr.scl_inter
end

# raw * slope + inter as Float32, which is what indexing a NIVolume computes
function _scale(raw::Array{T,5}, slope::Float32, inter::Float32) where T
    if T === Float32 && slope == 1 && inter == 0
        return raw
    end
    out = Array{Float32,5}(undef, size(raw))
    @inbounds for i in eachindex(raw, out)
        out[i] = raw[i] * slope + inter
    end
    return out
end

"""
    loadnii(filename) -> (data, header)

Reads a NIfTI file into a `Float32` array with five dimensions, `(x, y, z, echo,
channel)` with size 1 for the dimensions the file does not have, and returns it
with the file header. The header scaling is applied to the data and cleared in
the returned header, so the header describes the data and can be passed to
[`savenii`](@ref) as it is. Unlike [`niread`](@ref), the element type of the
result does not depend on the file, which is what static compilation with
`juliac --trim` requires.

See also [`loadphase`](@ref), [`loadmag`](@ref), [`loadheader`](@ref).
"""
loadnii(filename) = _read_nii(filename) do raw, hdr
    _scale(raw, _slope_inter(hdr)...)
end

"""
    loadphase(filename; rescale=true, fix_ge=false) -> (phase, header)

[`readphase`](@ref) as a `Float32` array with five dimensions, see [`loadnii`](@ref).
"""
function loadphase(filename; rescale=true, fix_ge=false)
    return _read_nii(filename) do raw, hdr
        slope, inter = _slope_inter(hdr)
        if rescale
            minr, maxr = approxextrema(raw)
            minp, maxp = minmax(Float32(minr * slope + inter), Float32(maxr * slope + inter))
            if !isapprox(maxp - minp, 2π; atol=0.1) # rescaling required
                minr, maxr = Float32.((minr, maxr))
                if isapprox(maxr - minr, 2π; atol=0.1) # no rescaling required, but header wrong
                    slope, inter = 1f0, 0f0
                else
                    slope = Float32(2pi / (maxr - minr))
                    inter = Float32(-pi - minr * slope)
                end
            end
        end
        if fix_ge
            fix_ge_phase!(raw)
            slope = -slope # phase is inverted
        end
        _scale(raw, slope, inter)
    end
end

"""
    loadmag(filename; rescale=false) -> (mag, header)

[`readmag`](@ref) as a `Float32` array with five dimensions, see [`loadnii`](@ref).
"""
function loadmag(filename; rescale=false)
    return _read_nii(filename) do raw, hdr
        slope, inter = _slope_inter(hdr)
        if rescale
            mini, maxi = Float32.(approxextrema(raw))
            slope = 1 / (maxi - mini)
            inter = -mini * slope
        end
        _scale(raw, slope, inter)
    end
end

# Writes `image` with `header`, byte-identical to `niwrite` of a `NIVolume` built
# from them, without NIfTI.jl's field-name iteration.
function _niwrite(filepath, header::NIfTI.NIfTI1Header, image::AbstractArray{T}) where T
    hdr = copy(header)
    hdr.dim = NIfTI.to_dim_i16(size(image))
    hdr.datatype = NIfTI.eltype_to_int16(T)
    hdr.bitpix = NIfTI.nibitpix(T)
    hdr.vox_offset = hdr.sizeof_hdr + 4 # no extensions
    bytes = _pack_header(hdr)
    open(filepath, "w") do io
        if endswith(filepath, ".gz")
            stream = GzipCompressorStream(io)
            _write_volume(stream, bytes, image)
            close(stream)
        else
            _write_volume(io, bytes, image)
        end
    end
end
_niwrite(filepath, header, image) = niwrite(filepath, NIVolume(header, image))

# Converts to the output element type and writes. One statically resolvable call
# per NIfTI element type, so that a compiled program does not need to dispatch
# on a type that is only known at run time.
_write_as(filepath, header, image, ::Nothing) = _niwrite(filepath, header, image)
function _write_as(filepath, header, image, T::Type)
    T === eltype(image) && return _niwrite(filepath, header, image)
    T === Float32 && return _niwrite(filepath, header, Float32.(image))
    T === Float64 && return _niwrite(filepath, header, Float64.(image))
    T === UInt8 && return _niwrite(filepath, header, UInt8.(image))
    T === Int16 && return _niwrite(filepath, header, Int16.(image))
    T === Int32 && return _niwrite(filepath, header, Int32.(image))
    T === Int64 && return _niwrite(filepath, header, Int64.(image))
    T === Int8 && return _niwrite(filepath, header, Int8.(image))
    T === UInt16 && return _niwrite(filepath, header, UInt16.(image))
    T === UInt32 && return _niwrite(filepath, header, UInt32.(image))
    T === UInt64 && return _niwrite(filepath, header, UInt64.(image))
    T === ComplexF32 && return _niwrite(filepath, header, ComplexF32.(image))
    T === ComplexF64 && return _niwrite(filepath, header, ComplexF64.(image))
    T === Bool && return _niwrite(filepath, header, Bool.(image))
    throw(ArgumentError("NIfTI cannot store this element type"))
end

function _write_volume(io, bytes, image)
    write(io, bytes)
    write(io, Int32(0)) # extender bytes: no extensions
    write(io, eltype(image) === Bool ? BitArray(image) : image)
end
