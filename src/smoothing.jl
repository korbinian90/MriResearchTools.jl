"""
    gaussiansmooth3d(image)

    gaussiansmooth3d(image, sigma=[5,5,5];
        mask=nothing,
        nbox=ifelse(isnothing(mask), 3, 6), 
        weight=nothing, dims=1:min(ndims(image),3), 
        boxsizes=getboxsizes.(sigma, nbox)
        )

Performs Gaussian smoothing on `image` with `sigma` as standard deviation of the Gaussian.
By application of `nbox` times running average filters in each dimension.
The length of `sigma` and the length of the `dims` that are smoothed have to match. (Default `3`)

Optional arguments:
- `mask`: Smoothing can be performed using a mask to inter-/extrapolate missing values.
- `nbox`: Number of box applications. Default is `3` for normal smoothing and `6` for masked smoothing.
- `weight`: Apply weighted smoothing. Either weighted or masked smoothing can be porformed.
- `dims`: Specify which dims should be smoothed. Corresponds to manually looping of the other dimensions.
- `boxizes`: Manually specify the boxsizes, not using the provided sigma. `length(boxsizes)==length(dims) && length(boxsizes[1])==nbox`
"""
gaussiansmooth3d, gaussiansmooth3d!

function gaussiansmooth3d(image, sigma=[5,5,5]; padding=false, kwargs...)
    if padding
        image = pad_image(image, sigma)
    end
    smoothed = gaussiansmooth3d!(0f0 .+ copy(image), sigma; kwargs...)
    if padding
        smoothed = remove_padding(smoothed, sigma)
    end
    return smoothed
end

pad_image(image, sigma) = PaddedView(0, image, Tuple(size(image) .+ 2sigma), Tuple(sigma .+ 1))
remove_padding(image, sigma) = image[[sigma[i]+1:size(image,i)-sigma[i] for i in 1:ndims(image)]...]

"""
    gaussiansmooth3d_phase(phase, sigma=[5,5,5]; weight=1, kwargs...)

Smoothes the phase via complex smoothing. A weighting image can be given.
The same keyword arguments are supported as in `gaussiansmooth3d`:
$(@doc gaussiansmooth3d)
"""
function gaussiansmooth3d_phase(phase, sigma=[5,5,5]; weight=1, kwargs...)
    clx = weight .* exp.(1im .* phase)
    phase_real = real.(clx)
    phase_imag = imag.(clx)
    gaussiansmooth3d!(phase_real, sigma; kwargs...)
    gaussiansmooth3d!(phase_imag, sigma; kwargs...)
    return angle.(complex.(phase_real, phase_imag))
end

function gaussiansmooth3d!(image, sigma=[5,5,5]; mask=nothing, nbox=ifelse(isnothing(mask), 3, 4), weight=nothing, dims=1:min(ndims(image),3), boxsizes=getboxsizes.(sigma, nbox))
    if length(sigma) < length(dims) @error "Length of sigma and dims does not match!" end
    if length(boxsizes) < length(dims) || length(boxsizes[1]) != nbox @error "boxsizes has wrong size!" end
    if typeof(mask) != Nothing
        image .*= ifelse.(mask .== 0, NaN, 1) # 0 in mask -> NaN in image
    end
    if typeof(weight) != Nothing
        weight = Float32.(weight)
        weight[weight .== 0] .= minimum(weight[weight .!= 0])
    end
    checkboxsizes!(boxsizes, size(image), dims)

    for ibox in 1:nbox, dim in dims
        bsize = boxsizes[dim][ibox]
        if size(image, dim) == 1 || bsize < 3
            continue
        end
        linefilter = getfilter(image, weight, mask, bsize, size(image, dim))
        K = ifelse(mask isa Nothing || isodd(ibox), :, size(image, dim):-1:1)

        # TODO parallel? -> Distributed arrays? -> use slices
        #loop = Iterators.product((size(image) |> sz -> (sz[1:(dim-1)], sz[(dim+1):end]) .|> CartesianIndices)...)
        #Threads.@threads for (I, J) in collect(loop)
            
        for J in CartesianIndices(size(image)[(dim+1):end])
            for I in CartesianIndices(size(image)[1:(dim-1)])
                w = if weight isa Nothing nothing else view(weight,I,:,J) end
                linefilter(view(image,I,K,J), w)
            end
        end
    end
    return image
end

## Calculate the filter sizes to achieve a given sigma

function getboxsizes(sigma, n)
    try
        wideal = √( (12sigma^2 / n) + 1 )
        wl::Int = round(wideal - (wideal + 1) % 2) # next lower odd integer
        wu::Int = wl + 2

        mideal = (12sigma^2 - n*wl.^2 - 4n*wl - 3n) / (-4wl - 4)
        m = round(mideal)

        [if i <= m wl else wu end for i in 1:n]
    catch
        zeros(n)
    end
end

function checkboxsizes!(boxsizes, sz, dims)
    for dim in dims
        bs = boxsizes[dim]
        for i in eachindex(bs)
            if iseven(bs[i])
                bs[i] += 1
            end
            if bs[i] > sz[dim] / 2
                val = sz[dim] ÷ 2
                if iseven(val) val += 1 end
                bs[i] = val
            end
        end
    end
end

## Function to initialize the filters

function getfilter(image, weight::Nothing, mask::Nothing, bsize, len)
    q = CircularBuffer{eltype(image)}(bsize)
    return (im, _) -> boxfilterline!(im, bsize, q)
end
function getfilter(image, weight, mask::Nothing, bsize, len)
    q = CircularBuffer{eltype(image)}(bsize)
    qw = CircularBuffer{eltype(weight)}(bsize)
    return (im, w) -> boxfilterline!(im, bsize, w, q, qw)
end
function getfilter(image, weight, mask, bsize, len)
    buffer = ones(eltype(image), len + bsize - 1) * NaN16
    return (im, _) -> nanboxfilterline!(im, bsize, buffer)
end

## Running Average Filters

function boxfilterline!(line::AbstractVector, boxsize::Int, q::CircularBuffer)
    r = div(boxsize, 2)
    initvals = view(line, 1:r)
    lsum = sum(initvals)
    append!(q, initvals)

    # start with edge effect
    @inbounds for i in 1:(r+1)
        lsum += line[i+r]
        push!(q, line[i+r])
        line[i] = lsum / (r + i)
    end

    # middle part
    @inbounds for i in (r+2):(length(line)-r)
        lsum += line[i+r] - popfirst!(q)
        push!(q, line[i+r])
        line[i] = lsum / boxsize
    end

    # end with edge effect
    @inbounds for i in (length(line)-r+1):length(line)
        lsum -= popfirst!(q)
        line[i] = lsum / (r + length(line) - i + 1)
    end
end

function boxfilterline!(line::AbstractVector, boxsize::Int, weight::AbstractVector, lq::CircularBuffer, wq::CircularBuffer)
    r = div(boxsize, 2)

    wsmooth = wsum = sum = eps() # slightly bigger than 0 to avoid division by 0
    @inbounds for i in 1:boxsize
        sum += line[i] * weight[i]
        wsum += weight[i]
        wsmooth += weight[i]^2
        push!(lq, line[i])
        push!(wq, weight[i])
    end

    @inbounds for i in (r+2):(length(line)-r)
        w = weight[i+r]
        l = line[i+r]
        wold = popfirst!(wq)
        lold = popfirst!(lq)
        push!(wq, w)
        push!(lq, l)

        sum += l * w - lold * wold
        wsum += w - wold
        line[i] = sum / wsum
        wsmooth += w^2 - wold^2
        weight[i] = wsmooth / wsum
    end
end

function nanboxfilterline!(line::AbstractVector, boxsize::Int, orig::AbstractVector)
    n = length(line)
    r = div(boxsize, 2)
    maxfills = r

    orig[r+1:r+n] .= line
    orig[r+n+1:end] .= NaN

    lsum = sum(view(orig,r+1:2r))
    if isnan(lsum) lsum = 0. end
    nfills = 0
    nvalids = 0

    mode = :nan

    @inbounds for i in eachindex(line)
        if isnan(lsum) @warn "lsum nan"; break end

        # check for mode change
        if mode == :normal
            if isnan(orig[i+2r])
                mode = :fill
            end

        elseif mode == :nan
            if isnan(orig[i+2r])
                nvalids = 0
            else
                nvalids += 1
            end
            if nvalids == boxsize
                mode = :normal
                lsum = sum(view(orig,i:(i+2r)))
                line[i] = lsum / boxsize
                continue # skip to next loop iteration
            end

        elseif mode == :fill
            if isnan(orig[i+2r])
                nfills += 1
                if nfills > maxfills
                    mode = :nan
                    nfills = 0
                    lsum = 0
                    nvalids = 0
                end
            else
                mode = :normal
                nfills = 0
            end
        end

        # perform operation
        if mode == :normal
            lsum += orig[i+2r] - orig[i-1]
            line[i] = lsum / boxsize
        elseif mode == :fill
            lsum -= orig[i-1]
            line[i] = (lsum - orig[i]) / (boxsize - 2)

            orig[i+2r] = 2line[i] - line[i-r] # TODO maybe clamp the value
            if (i+r < n) line[i+r] = orig[i+2r] end
            lsum += orig[i+2r]
        end

    end
end
