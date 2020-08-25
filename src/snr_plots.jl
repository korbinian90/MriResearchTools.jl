## INIT
YNORM = 1
plot_type(args...; keyargs...) = plot_type!(plot(), args...; keyargs...)
function plot_type!(h, TRs, Nx, FOV, snr_type; field=:B7T, equidistant=0, ME=:bipolar, NE=(3:2:16)', label=["$ME: $i echoes" for i in NE], kw...)
    t = getTimes.(Ref(TRs), Nx, FOV, NE, ME, field, equidistant)
    if NE != 1
        t = hcat(t...)
    end
    plot_type!(h, t, snr_type; field=field, label=label, kw...)
end
function plot_type!(h, t, snr_type; field=:B7T, y_norm=:off, kw...)
    global YNORM
    y = snr_ratio.(t, Ref(snr_type), field)
    #y[.!isfinite.(y)] .= 0

    if y_norm == :calculate
        YNORM = maximum(filter(!isnan, y))
        y_norm = :on
    end
    if y_norm == :on
        y ./= YNORM
    end

    plot!(h, getfield.(t, :TR), y; kw...) # Ref() to avoid broadcasting
end

## SNR Graph Plot
function snr_graph()
    f = :B7T
    NE = [3 5 7 9]
    Nx = 800
    FOV = 210
    TR = 0:0.1:70
    snrtype = :SNR => :wm
    config = (TR, Nx, FOV, snrtype)
    config_me = (equidistant=NE, y_norm=:on, field=f, NE=NE, linecolor=(2:5)')
    h = plot_type(config...;y_norm=:calculate, field=f, NE=1, label="single-echo", xaxis="TR", yaxis="nSNR", legend=:bottomright)
    h = plot_type!(h, config...; config_me..., ME=:monopolar, linestyle=:dash, label=nothing)
    h = plot_type!(h, config...; config_me..., ME=:bipolar, label=["$i echoes" for i in NE])
    #=
    figpath = raw"F:\MRI\paperfigures\simulation"
    mkpath(figpath)
    for type in ["png", "pdf", "svg"]
        savefig(h, joinpath(figpath, "snr_equidistant.$type"))
    end
    =#
    return h
end

## Dutycycle Table
function dutycycle_table()
    FOV = 210 # in [m]
    Nx = [64 128 256 512 768 960]
    ME = [:monopolar, :bipolar]
    NE = 8
    ΔTE = 4

    times = findtimesforΔTE.(Nx, NE, ME, FOV, ΔTE)
    return getdutycycle.(times)
end

#=
## TABLE SNR (abstract)
s = SCANNER[:B7T]

y = snr_ratio.(times, :SNR => :wm, :B7T)
times_SE = getTimes.(getfield.(times, :TR), Nx, FOV, 1)
se = snr_ratio.(times_SE, :SNR => :wm, :B7T)
@show y ./ se

TR_se = times[1].TR
Tacq_se = TR_se - s[:Texc] - s[:Tspoil]
TE_se = s[:Texc] + Tacq_se / 2
times_SE = Times(Tacq_se, [TE_se], TR_se)
snr_se = snr_ratio(times_SE, :SNR => :wm, :B7T)

times_perfect = Times(4, 4:4:32, 36)
snr_perfect = snr_ratio(times_perfect, :SNR => :wm, :B7T)

@show snr_table = sqrt.(dutycycle) .* (snr_perfect / snr_se)
=#

#=
## add REAL measurement point to plot
include("../datasets.jl")
d = loadconfig(:tumor_31)
snrtype = :SNR => :wm

Tacq = 2getTacq(d[:Tdwell], d[:Nx])
t = Times(Tacq, d[:TE], d[:TR])
y = ratio_me(t, snrtype, f)

#t = getTimes.(31, 736, 220, 6)

plot_type!(h, [t], snrtype; field=:B7T, y_norm=:on, marker=:dot, label=d[:tag])
=#
