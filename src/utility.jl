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

function approxextrema(I)
    arr = sample(I)
    return (minimum(arr), maximum(arr))
end

function estimatequantile(array, p)
    try 
        return quantile(sample(array; n=1e5), p)
    catch
        @warn "quantile could not be estimated! (maybe only NaNs)"
        return 0
    end
end

function sample(I; n=10000)
    n = min(n, length(I))
    len = ceil(Int, √n) # take len blocks of len elements
    startindices = round.(Int, range(firstindex(I) - 1, lastindex(I) - len; length=len))
    indices = vcat((i .+ (1:len) for i in startindices)...)
    ret = filter(isfinite, I[indices])
    if isempty(ret)
        ret = filter(isfinite, I)
    end
    return ret
end

savenii(image, name, writedir::Nothing, header=nothing) = nothing
function savenii(image, name, writedir, header=nothing)
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
function savenii(image::AbstractArray, filepath; header=nothing)
    vol = NIVolume([h for h in [header] if h !== nothing]..., image)
    niwrite(filepath, vol)
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

function getHIP(compl; echoes=[1,2])
    e1, e2 = echoes
    c = zeros(eltype(compl), size(compl)[1:3])
    for iCha in 1:size(compl, 5)
        c .+=  compl[:,:,:,e2,iCha] .* conj.(compl[:,:,:,e1,iCha])
    end
    return c
end

function get_corner_indices(I; max_length=10)
    d = size(I)
    n = min.(max_length, ceil.(Int, d ./ 3)) # n voxels for each dim
    getrange(num, len) = [1:num, (len-num+1):len] # first and last voxels
    return collect(Iterators.product(getrange.(n, d)...))
end

# estimate noise parameters from corner without signal
function estimatenoise(image)
    corners = get_corner_indices(image)
    (lowestmean, ind) = findmin(mean.(filter(isfinite, image[I...]) for I in corners))
    sigma = std(filter(isfinite, image[corners[ind]...]))
    if isnan(sigma) # no corner available
        # estimation that is only true if have the image is signal and half noise
        sigma = 2estimatesigma_from_quantile(image, 1/4)
        lowestmean = sigma / 2
    end
    return lowestmean, sigma
end

# sigma is only calculated for quantile (biased)
function estimatesigma_from_quantile(image, quantile)
    q = estimatequantile(image, quantile)
    samples = filter(x -> x < q, sample(image))
    return std(samples)
end

function robustmask!(image; maskedvalue=if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!robustmask(image)] .= maskedvalue
    image
end
function robustmask(weight)
    μ, σ = estimatenoise(weight)
    m = mean(filter(isfinite, weight[weight .> μ + 4σ]))
    return weight .> maximum((μ + 3σ, m/5))
end


getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)
getcomplex(fnmag::AbstractString, fnphase::AbstractString) = getcomplex(niread(fnmag), niread(fnphase))

function getcomplex(mag, phase)
    minp, maxp = approxextrema(phase)
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
    array[isnan.(array)] .= minimum(filter(isfinite, array))
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

to_dim(V::AbstractVector, dim::Int) = reshape(V, ones(Int, dim-1)..., :)
