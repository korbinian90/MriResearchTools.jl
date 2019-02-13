#module Utility

#using Statistics, Images
#using NIfTI_mod
using Statistics

export Data, readphase, savenii, createniiforwriting, getHIP, getrobustmask, robustmask!, getcomplex, readfromtextheader, combine_echoes, getscaledimage, RSS, rescale!, rescale, robustrescale, robustrescale!

struct Data
    path
    dirphase
    fnphase
    dirmag
    fnmag
    name
    TEs
    rbw
    readoutdim
end
function Data(;path = "", dirphase = "", fnphase = joinpath(path, dirphase, "Image.nii"), dirmag = "", fnmag = joinpath(path, dirmag, "Image.nii"), name, TEs = [], rbw = 0, readoutdim = 0)
    Data(path, dirphase, fnphase, dirmag, fnmag, name, TEs, rbw, readoutdim)
end

function readphase(fn; keyargs...)
    phase = niread(fn; keyargs...)
    if eltype(phase.raw) <: Integer && phase.header.scl_slope == 0
        minp, maxp = minmaxfirstslice(phase.raw)
        phase.header.scl_slope = 2pi / (maxp - minp + 1)
        phase.header.scl_inter = -pi - minp * phase.header.scl_slope
    end
    return phase
end

function minmaxfirstslice(image)
    ones = repeat([1], ndims(image)-2)
    image[:,:,ones...] |> I -> (minimum(I), maximum(I))
end

savenii(image, name, writedir; kwargs...) = savenii(image, joinpath(writedir, name * ".nii"); kwargs...)
savenii(image, filepath; kwargs...) = niwrite(filepath, NIVolume([k[2] for k in kwargs]..., image))
savenii(image::BitArray, filepath; kwargs...) = niwrite(filepath, NIVolume([k[2] for k in kwargs]..., UInt8.(image)))
#savenii(image, filepath) = niwrite(filepath, NIVolume(image))

function createniiforwriting(im, name::AbstractString, writedir::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header)
    filepath = joinpath(writedir, name * ".nii")
    createniiforwriting(im, filepath; datatype = datatype, header = header)
end

function createniiforwriting(im::AbstractArray, filepath::AbstractString; datatype::DataType = eltype(im), header = NIVolume(zeros(datatype, 1)).header)
    nii = createniiforwriting(size(im), filepath; datatype = datatype, header = header)
    nii .= im
end

createniiforwriting(sz::Tuple, filepath::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header) = createniiforwriting([sz...], filepath; datatype = datatype, header = header)
function createniiforwriting(sz::AbstractVector{Int}, filepath::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header)
    write_emptynii(sz, filepath, datatype = datatype, header = header)
    niread(filepath, mmap = true, write = true).raw
end

function getHIP(mag, phase; echoes = [1, 3])
    minp = minimum(phase[:,:,1,1,1])
    maxp = maximum(phase[:,:,1,1,1])
    e1, e2 = echoes
    dropdims(sum(exp.((2im * pi / (maxp - minp)) .* (Float64.(phase[:,:,:,e2,:]) .- phase[:,:,:,e1,:])) .* Float64.(mag[:,:,:,e1,:]) .* mag[:,:,:,e2,:]; dims = 5); dims = 4)
end


function robustmask!(image; maskedvalue = if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!getrobustmask(image)] .= maskedvalue
    image
end

function getrobustmask(weight)
    noisemask = weight .<= mean(weight)
    noisemean = mean(weight[noisemask])

    signalmean = mean(weight[.!noisemask])
    noise_σ = std(weight[noisemask])

    weight .> noisemean + noise_σ
end

function getcomplex(mag, phase)
    higherdims = ones(Int, length(size(phase)) - 2)
    minp = minimum(phase[:,:,higherdims...])
    maxp = maximum(phase[:,:,higherdims...])

    mag .* exp.((2im * pi / (maxp - minp)) .* phase)
end

getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)

getcomplex(fnmag::AbstractString, fnphase::AbstractString) = getcomplex(niread(fnmag), niread(fnphase))

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

combine_echoes(mag) = RSS(mag)

function combine_echoes(phase, mag, TEs)
    TEs = reshape(TEs, ones(Int, ndims(phase)-1)..., length(TEs))
    dim = ndims(phase)
    combined = dropdims(sum(phase .* mag .* TEs; dims = dim); dims = dim)
    combined ./= dropdims(sum(mag .* TEs; dims = dim); dims = dim)
end


function getscaledimage(array, div::Number, offset = 0, type::Symbol = :trans)
    array = reshape(array, size(array)[1:2]) # drops singleton dimensions
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
    scaled .|> Gray
end

function getscaledimage(array, type::Symbol = :trans)
    scaled = robustrescale(array, 0, 1, threshold = true)
    getscaledimage(scaled, 1, 0, type)
end

function estimatequantile(array, p)
    # only use 5 samples of each higher dimension
    # TODO take 100 x 16 consecutive with linear indexing (bound checking!)
    inds = [unique(round.(Int, LinRange(1, sz, 5))) for sz in size(array)[2:end]]
    quantile(array[:, inds...][:], p)
end

function robustrescale(array, newmin, newmax; threshold = false, mask = trues(size(array)), datatype = Float64)
    robustrescale!(datatype.(array), newmin, newmax; threshold = threshold, mask = mask)
end
#=
function robustrescalefactors(array, newmin, newmax; mask = trues(size(array)))
    q = [0.01, 0.99] # quantiles
    oldq = estimatequantile(array[mask], q)
    oldrange = (oldq[2] - oldq[1]) / (q[2] - q[1])
    oldmin = oldq[1] - q[1] * oldrange
    newrange = newmax - newmin
end
=#
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

function rescale(array, newmin, newmax; datatype = eltype(array))
    rescale!(datatype.(array), newmin, newmax)
end

function rescale!(array, newmin, newmax)
    # TODO sample function to take only every.. elem block
    oldmin, oldmax = extrema(array)
    factor = (newmax - newmin) / (oldmax - oldmin)
    array .= (array .- oldmin) .* factor .+ newmin
end

mmtovoxel(sizemm, nii::NIVolume) = mmtovoxel(sizemm, nii.header)
mmtovoxel(sizemm, header::NIfTI1Header) = mmtovoxel(sizemm, header.pixdim)
mmtovoxel(sizemm, pixdim) = sizemm ./ pixdim

#end
