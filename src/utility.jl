function approxextrema(I)
    arr = sample(I)
    return (minimum(arr), maximum(arr))
end

"""
    estimatequantile(array, p)

Quickly estimates the quantile `p` of a possibly large array by using a subset of the data.
"""
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

function get_corner_indices(I; max_length=10)
    d = size(I)
    n = min.(max_length, ceil.(Int, d ./ 3)) # n voxels for each dim
    getrange(num, len) = [1:num, (len-num+1):len] # first and last voxels
    return collect(Iterators.product(getrange.(n, d)...))
end

# estimate noise parameters from corner without signal
"""
    estimatenoise(image::AbstractArray)

Estimates the noise from the corners of the image.
It assumes that at least one corner is without signal and only contains noise.
"""
function estimatenoise(image::AbstractArray)
    corners = get_corner_indices(image)
    (lowestmean, ind) = findmin(mean.(filter(isfinite, image[I...]) for I in corners))
    sigma = std(filter(isfinite, image[corners[ind]...]))
    if isnan(sigma) # no corner available
        # estimation that is correct if half the image is signal and half noise
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


"""
    robustmask!(image)
    robustmask!(image; maskedvalue)

Creates a mask and applies it inplace.
It assumes that at least one corner is without signal and only contains noise.
"""
function robustmask!(image; maskedvalue=if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!robustmask(image)] .= maskedvalue
    image
end
"""
    robustmask(weight::AbstractArray)

Creates a mask from a weights images.
It assumes that at least one corner is without signal and only contains noise.
"""
function robustmask(weight::AbstractArray)
    μ, σ = estimatenoise(weight)
    m = mean(filter(isfinite, weight[weight .> 5σ]))
    mask = weight .> maximum((5σ, m/5))
    # remove small holes and minimally grow
    boxsizes=[[3,3] for i in 1:ndims(weight)]
    return gaussiansmooth3d(mask; nbox=2, boxsizes) .> 0.55
end

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
"""
    RSS(mag; dims=ndims(mag))

Performs root-sum-of-squares combination along the last dimension of `mag`.
The dimension can be specificed via the `dims` keyword argument.
"""
RSS(mag; dims=ndims(mag)) = dropdims(.√sum(mag.^Float32(2); dims); dims)

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
    scaled = robustrescale(array, 0, 1, threshold=true)
    getscaledimage(scaled, 1, 0, type)
end

"""
    robustrescale(array, newmin, newmax; threshold=false, mask=trues(size(array)), datatype=Float64)

Rescales the image to the the new range, disregarding outliers.
Only values inside `mask` are used for estimating the rescaling option
"""
robustrescale(array, newmin, newmax; threshold=false, mask=trues(size(array)), datatype=Float64) =
    robustrescale!(datatype.(array), newmin, newmax; threshold, mask)

function robustrescale!(array, newmin, newmax; threshold=false, mask=trues(size(array)))
    mask[isnan.(array)] .= false
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


function rescale(array, newmin, newmax; datatype=eltype(array))
    rescale!(datatype.(array), newmin, newmax)
end

function rescale!(array, newmin, newmax)
    oldmin, oldmax = approxextrema(array)
    factor = (newmax - newmin) / (oldmax - oldmin)
    array .= (array .- oldmin) .* factor .+ newmin
end

"""
    to_dim(V::AbstractVector, dim::Int)

    to_dim(a::Real, dim::Int)

Converts a vector or number to a higher dimension.

# Examples
```julia-repl
julia> to_dim(5, 3)
1×1×1 Array{Int64, 3}:
[:, :, 1] =
 5
julia> to_dim([1,2], 2)
 1×2 Matrix{Int64}:
  1  2
```
"""
to_dim(a::Real, dim::Int) = to_dim([a], dim)
to_dim(V::AbstractVector, dim::Int) = reshape(V, ones(Int, dim-1)..., :)

"""
    getHIP(mag, phase; echoes=[1,2])

    getHIP(compl; echoes=[1,2])

Calculates the Hermitian Inner Product between the specified echoes.
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
