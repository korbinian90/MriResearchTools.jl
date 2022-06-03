struct Ice_output_config
    name::String
    path::String
    nslices::Int
    nfiles::Int
    nechoes::Int
    nchannels::Int
    dtype::Type
    size::Tuple{Int, Int}
end

function Ice_output_config(name, path, nslices, nfiles; nechoes=1, nchannels=1, dtype=Int16)
    return Ice_output_config(name, path, nslices, nfiles, nechoes, nchannels, dtype, getsize(path))
end

function get_setting(T, lines, setting; offset=3, default=0)
    for iLine in 1:length(lines)
        if occursin(setting, lines[iLine])
            try
                return parse(T, lines[iLine + offset])
            catch
                return default
            end
        end
    end
    return default
end

function read_volume(cfg)
    volume = create_volume(cfg)

    for i in 1:cfg.nfiles
        num = lpad(i, 5, '0')
        imahead = joinpath(cfg.path, "MiniHead_ima_$num.IceHead")
        file = joinpath(cfg.path, "WriteToFile_$num.ima")
        
        if occursin(cfg.name, read(imahead, String))
            vol = Array{cfg.dtype}(undef, cfg.size...)
            read!(file, vol)
            lines = readlines(imahead)
            eco = get_setting(Int, lines, "EchoNumber"; default=1)
            slc = getslice(lines)
            rescale_slope = get_setting(Float32, lines, "RescaleSlope"; offset=4, default=1)
            rescale_intercept = get_setting(Float32, lines, "RescaleIntercept"; offset=4, default=0)
            volume[:,:,slc,eco] .= (vol .+ rescale_intercept) .* rescale_slope
        end
    end

    return volume
end

function getslice(lines)
    slc = get_setting(Int, lines, "Actual3DImaPartNumber"; default=nothing)
    if (isnothing(slc))
        slc = get_setting(Int, lines, "AnatomicalSliceNo")
    end
    return slc + 1
end

function getsize(path)
    lines = readlines(joinpath(path, "MiniHead_ima_00001.IceHead"))
    return (get_setting(Int, lines, "NoOfCols"), get_setting(Int, lines, "NoOfRows"))
end

function create_volume(cfg)
    return zeros(Float32, cfg.size..., cfg.nslices, cfg.nechoes)
end