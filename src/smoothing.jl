# mask should be 3D
# image can have any higher dimension

function gaussiansmooth3d(image, σ=[5,5,5]; kwargs...)
    gaussiansmooth3d!(0f0 .+ copy(image), σ; kwargs...)
end

function gaussiansmooth3d!(image, σ=[5,5,5]; mask=nothing, nbox=4, weight=nothing, dims=1:ndims(image), boxsizes=nothing)
    if σ isa Number
        σ = σ * ones(ndims(image))
    end
    if typeof(mask) != Nothing
        nbox *= 2
        # TODO do we need small boxsize?
        #@show boxsizes = getboxsizes_small.(σ, nbox, 5)
        #boxsizes = getboxsizes.(σ, nbox)
        image[mask .== 0] .= NaN
    end
    if typeof(weight) != Nothing
          w = Float32.(weight)
          w[w .== 0] .= minimum(w[w .!= 0])
    end
    if boxsizes === nothing boxsizes = getboxsizes.(σ, nbox) end
    checkboxsizes!(boxsizes, size(image), dims)

    for ibox in 1:nbox, dim in dims
        bsize = boxsizes[dim][ibox]
        if size(image, dim) == 1 || bsize < 3
            continue
        end
        q = CircularBuffer{eltype(image)}(bsize)
        q2 = CircularBuffer{eltype(image)}(bsize)
        K = ifelse(isodd(ibox), :, size(image, dim):-1:1)
        # TODO parallel? -> Distributed arrays?
        #loop = Iterators.product((size(image) |> sz -> (sz[1:(dim-1)], sz[(dim+1):end]) .|> CartesianIndices)...)
        #Threads.@threads for (I, J) in collect(loop)
        for I in CartesianIndices(size(image)[1:(dim-1)])
            for J in CartesianIndices(size(image)[(dim+1):end])
                if typeof(mask) != Nothing
                    im = view(image, I, K, J)
                    nanboxfilterline!(im, boxsizes[dim][ibox])
                elseif typeof(weight) != Nothing
                    boxfilterline!(view(image,I,:,J), boxsizes[dim][ibox], view(w,I,:,J))
                else
                    boxfilterline!(view(image,I,:,J), boxsizes[dim][ibox], q)
                end
            end
        end
    end
    image
end

function getboxsizes(σ, n)
    try
        wideal = √( (12σ^2 / n) + 1 )
        wl::Int = round(wideal - (wideal + 1) % 2) # next lower odd integer
        wu::Int = wl + 2

        mideal = (12σ^2 - n*wl.^2 - 4n*wl - 3n) / (-4wl - 4)
        m = round(mideal)

        [if i <= m wl else wu end for i in 1:n]
    catch
        zeros(n)
    end
end

# TODO compare with MATLAB
function getboxsizes_small(σ, n::Int, smallsize::Int)
    smallat = [3; 4]
    nsmall = length(smallat)

    smallsize = 2ceil(smallsize / 2) - 1

    wideal = √( (12σ^2 - (smallsize^2 - 1)nsmall) / (n - nsmall) + 1 )

    wl::Int = wideal - (wideal + 1) % 2 # next lower odd integer
    wu::Int = wl + 2

    mideal = (12σ^2 - (n - nsmall) * (wu^2 - 1) - (smallsize^2 - 1)nsmall) / (wl^2 - wu^2)
    m = round(mideal)

    boxsizes = [if i <= m wl else wu end for i in 1:(n - nsmall)]
    for ismall in smallat
        insert!(boxsizes, ismall, smallsize)
    end
    boxsizes
end

function checkboxsizes!(boxsizes, sz, dims)
    for dim in dims
        bs = boxsizes[dim]
        for i in eachindex(bs)
        if iseven(bs[i])
            @warn "boxsize $i is even: $(bs[i]); it was changed to next bigger odd integer!"
            bs[i] += 1
        end
            if bs[i] > sz[dim] / 2
                @warn "boxsize $i is limited to half the image; it was changed from $(bs[i]) to $(sz[dim]÷2)!"
                bs[i] = sz[dim] ÷ 2
    end
end
    end
end

function boxfilterline!(line::AbstractVector, boxsize::Int, q::CircularBuffer)
    r = div(boxsize, 2)
    initvals = view(line, 1:boxsize)
    lsum = sum(initvals)
    append!(q, initvals)

    @inbounds for i in (r+2):(length(line)-r)
        lsum += line[i+r] - popfirst!(q)
        push!(q, line[i+r])
        line[i] = lsum / boxsize
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

function nanboxfilterline!(line::AbstractVector, boxsize::Int)
    n = length(line)
    r = div(boxsize, 2)
    maxfills = r

    orig = [repeat([NaN], r); line; repeat([NaN], r)]


    lsum = sum(orig[r+1:2r])
    if isnan(lsum) lsum = 0. end
    nfills = 0
    nvalids = 0

    mode = :nan

    @inbounds for i in 1:length(line)
        # TODO remove isnan check if it runs stable
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
                lsum = sum(orig[i:(i+2r)])
                line[i] = lsum / boxsize
                continue
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

            orig[i+2r] = 2line[i] - line[i-r]
            if (i+r < n) line[i+r] = orig[i+2r] end
            lsum += orig[i+2r]
        end

    end
end
