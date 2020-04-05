function readphase(fn; keyargs...)
    phase = niread(fn; keyargs...)
    minp, maxp = approxextrema(phase)
    phase.header.scl_slope = 2pi / (maxp - minp)
    phase.header.scl_inter = -pi - minp * phase.header.scl_slope
    return phase
end

function readmag(fn; normalize=false, keyargs...)
    mag = niread(fn; keyargs...)
    if mag.header.scl_slope == 0 || normalize
        _, maxi = approxextrema(mag)
        mag.header.scl_slope = 1 / maxi
        mag.header.scl_inter = 0
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

header(v::NIfTI.NIVolume) = similar(v.header)

Base.minimum(I::Array{AbstractFloat}) = NaNMath.minimum(I)
Base.maximum(I::Array{AbstractFloat}) = NaNMath.maximum(I)

approxextrema(image::NIVolume) = image.raw[1:30:end] |> I -> (minimum(I), maximum(I)) # sample every 30th value

savenii(image, name, writedir::Nothing, header=nothing) = nothing
function savenii(image, name, writedir::String, header=nothing)
    if splitext(name)[2] != ".nii"
        name = name * ".nii"
    end
    savenii(image, joinpath(writedir, name); header=header)
end
"""
    savenii(image, filepath; header=nothing)
save the image at the path
Warning: MRIcro can only open images with types Int32, Int64, Float32, Float64
"""
function savenii(image::AbstractArray, filepath::AbstractString; header=nothing)
    vol = NIVolume([h for h in [header] if h != nothing]..., image)
    niwrite(filepath, vol)
end

function write_emptynii(sz, path; datatype=Float64, header=NIVolume(zeros(datatype, 1)).header)
    header = copy(header)
    header.dim = Int16.((length(sz), sz..., ones(8-1-length(sz))...))
    header.datatype = NIfTI.nidatatype(datatype)
    header.bitpix = NIfTI.nibitpix(datatype)

    if isfile(path) rm(path) end
    file = open(path, "w")
    write(file, header)
    close(file)
end

"""
    getHIP(mag, phase; echoes=[1,2])
return the hermitian inner product between the specified echoes.
"""
function getHIP(mag, phase; echoes=[1,2])
    e1, e2 = echoes
    compl = zeros(ComplexF64, size(mag)[1:3])
    for iCha in 1:size(mag, 5)
        compl .+= exp.(1.0im .* (phase[:,:,:,e2,iCha] .- phase[:,:,:,e1,iCha])) .* mag[:,:,:,e1,iCha] .* mag[:,:,:,e2,iCha]
    end
    compl
end


function estimatenoise(weight)
    # find corner with lowest intensity
    d = size(weight)
    n = min.(10, d .÷ 3) # take 10 voxel but maximum a third
    getrange(num, len) = [1:len, (len-num+1):len] # first and last voxels
    corners = Iterators.product(getrange.(n, d)...)
    lowestmean = Inf
    sigma = 0
    for I in corners
        m = mean(weight[I...])
        if m < lowestmean
            lowestmean = m
            sigma = std(weight[I...])
        end
    end
    return lowestmean, sigma
end

function robustmask!(image; maskedvalue=if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!robustmask(image)] .= maskedvalue
    image
end
function robustmask(weight)
    μ, σ = estimatenoise(weight)
    return weight .> μ + 3σ
end


getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)
getcomplex(fnmag::AbstractString, fnphase::AbstractString) = getcomplex(niread(fnmag), niread(fnphase))

function getcomplex(mag, phase)
    higherdims = ones(Int, length(size(phase)) - 2)
    minp = minimum(phase[:,:,higherdims...])
    maxp = maximum(phase[:,:,higherdims...])

    mag .* exp.((2im * pi / (maxp - minp)) .* phase)
end

function readfromtextheader(filename, searchstring)
    for line in readlines(open(filename, "r"))
        if occursin(searchstring, line)
            # regex to search for "= " or ": " and return the following non-whitespace characters
            return match(r"(?<=(= |: ))(\S+)", line).match
        end
    end
end

# root sum of squares combination
RSS(mag; dim = ndims(mag)) = dropdims(.√sum(mag.^Float32(2); dims = dim); dims = dim)

function getscaledimage(array, div::Number, offset = 0, type::Symbol = :trans)
    array = reshape(array, size(array)[1:2]) # drops trailing singleton dimensions
    scaled = if offset != 0
        (array .- offset) .* (1 / div) .+ 0.5
    else
        array .* (1 / div)
    end
    scaled[isnan.(scaled) .| (scaled .< 0)] .= 0
    scaled[scaled .> 1] .= 1
    if type == :trans
        scaled = reverse(permutedims(scaled, [2 1]); dims = 1)
    else
    end
    scaled
end

function getscaledimage(array, type::Symbol = :trans)
    scaled = robustrescale(array, 0, 1, threshold = true)
    getscaledimage(scaled, 1, 0, type)
end

robustrescale(array, newmin, newmax; threshold = false, mask = trues(size(array)), datatype = Float64) =
    robustrescale!(datatype.(array), newmin, newmax; threshold = threshold, mask = mask)

function robustrescale!(array, newmin, newmax; threshold = false, mask = trues(size(array)))
    array[isnan.(array)] .= minimum(array[.!isnan.(array)])
    q = [0.01, 0.99] # quantiles
    oldq = estimatequantile(array[mask], q)
    oldrange = (oldq[2] - oldq[1]) / (q[2] - q[1])
    oldmin = oldq[1] - q[1] * oldrange
    newrange = newmax - newmin

    array .= (array .- oldmin) .* (newrange / oldrange) .+ newmin

    if threshold
        array[array .< newmin] .= newmin
        array[array .> newmax] .= newmax
    end
    array
end

function estimatequantile(array, p)
    samples = 1e5
    step = ceil(Int, length(array) / samples)
    quantile(array[1:step:end], p)
end

function rescale(array, newmin, newmax; datatype = eltype(array))
    rescale!(datatype.(array), newmin, newmax)
end

function rescale!(array, newmin, newmax)
    oldmin, oldmax = approxextrema(array)
    factor = (newmax - newmin) / (oldmax - oldmin)
    array .= (array .- oldmin) .* factor .+ newmin
end

mmtovoxel(sizemm, nii::NIVolume) = mmtovoxel(sizemm, nii.header)
mmtovoxel(sizemm, header::NIfTI.NIfTI1Header) = mmtovoxel(sizemm, header.pixdim)
mmtovoxel(sizemm, pixdim) = sizemm ./ pixdim
