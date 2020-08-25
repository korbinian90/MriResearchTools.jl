struct Times
    Tacq
    TE
    TR
end

## SNR Theory

function S(Tacq, TE, TR, T2s, T1, PD=1, α=ernst(TR))
    return PD * angleterm(TR, T1, α) * signalterm(T2s, Tacq) * t2s_decay(TE, T2s)
end

ernst(TR; T1=1900) = acos(exp(-TR/T1))
t2s_decay(TE, T2s) = exp(-TE/T2s)
signalterm(T2s, Tacq) = (2T2s * (1 - exp(-Tacq/2T2s))) # close to Tacq for Tacq << T2s
function angleterm(TR, T1, α=ernst(TR))
    E1 = exp(-TR/T1)
    angleterm = ((1 - E1) * sin(α)) / (1 - E1 * cos(α))
end

snr(args...) = S(args...) / noise(first(args)) # first(args) is Tacq
noise(Tacq) = √Tacq
combine_snr(snr) = √sum(snr.^2)

## Convenience functions

S(t::Times, tissue::Symbol, field) = S(t, gettissue(field, tissue)...)
S(t::Times, args...) = S(t.Tacq, t.TE, t.TR, args...)
noise(times::Times) = noise(times.Tacq)

## Create tables for plotting

function ratio_type(times, tissue, type, field)
    times = Times.(times.Tacq, times.TE, times.TR)
    if type == :CNR
        T2s, T1, PD, factor = gettissue_all(field)
        f = factor[tissue[1]] / factor[tissue[2]]
        easy_contrast.(times, f, T2s[tissue[1]], T2s[tissue[2]]) ./ noise.(times)
    elseif type == :SNR
        snr.(times, tissue, field)
    elseif type == :CNR_phase
        snr.(times, tissue, field) .* getfield.(times, :TE)
    end
end

function snr_ratio(times::Times, type, field; normalize=true)
    type, tissue = type
    r = combine_snr(ratio_type(times, tissue, type, field))
    if normalize
        r /= √times.TR
    end
    return r
end

## Simulation

function getTimes(TRret, Nx, FOV, NE=3, ME=:bipolar, field=:B7T, equidistant=true)
    TEs, TR, Tacq = buildTable(FOV, Nx, NE, ME, field, equidistant)
    t = Times.(0, Ref(zeros(NE)), TRret) # dummy times
    j = length(TR)
    for iret in 1:length(TRret)
        while TR[j-1] ≤ TRret[iret]
            if j == 1 continue end
            j -= 1
        end
        if TR[j] > TRret[iret] continue end
        ttemp = Times(Tacq[j], TEs[j], TRret[iret])
        if typeof(t) <: AbstractArray
            t[iret] = ttemp
        else
            return ttemp
        end
    end
    return t
end

function buildTable(FOV, Nx, NE, ME, field, equidistant; te_digits=2, points=100000)
    s = SCANNER[field]
    FOV /= 1e3 # in [m]
    GMsampling = Nx / s[:γ] / FOV
    # calculate for all possible Gradients
    pow = 4 # for more dense sampling for low Gradients
    times = Array{Times}(undef, points)
    for (i, Gs) in enumerate(range(0, s[:Glimit]^(1/pow); length=points))
        times[i] = getTEs_TR_Tacq(Gs^pow, GMsampling, s, NE, ME, field, te_digits, equidistant)
    end
    return getfield.(times, :TE), getfield.(times, :TR), getfield.(times, :Tacq) # matching arrays for all possible Gs
end

function getTEs_TR_Tacq(G, GMsampling, s, NE, ME, field, digits, equidistant)
    r(x) = round(x, RoundUp; digits=digits) # round to precision available on special card
    Tacq, Techo, Tprewind = 1e3 .* getTecho(G, s[:SlewR], GMsampling, s[:Glimit], ME)
    TE1 = Techo/2 + s[:Texc]
    ΔTE = r(Techo + Tprewind)
    TEs = if NE == 1
        [r(TE1)]
    elseif equidistant > 1
        TE2 = max(2TE1, 2ΔTE)
        TE1 = TE2 / 2
        TEs = r(TE1) .* collect(1:NE) # equidistant
        for i in (equidistant+1):NE
            TEs[i] = TEs[i-1] + ΔTE # replace with closest echo
        end
        TEs
    else # closest echoes
        TE1 .+ (0:NE-1) .* ΔTE
    end
    TR = TEs[NE] + Techo/2 + s[:Tspoil]
    return Times(Tacq, r.(TEs), r(TR))
end

# in [s]
function getTecho(G, SlewR, GMsampling, Glimit, ME)
    Tramp = G / SlewR
    Tsample = GMsampling / G
    GMtotal = Tramp * G + GMsampling
    return Tsample, 2Tramp + Tsample, getTprewind(GMtotal, SlewR, Glimit, ME)
end

# in [s]
function getTprewind(GM, SlewR, G_limit, ME)
    if ME == :bipolar return 0 end
    tPrewind = sqrt(4 * GM / SlewR)
    maxG = tPrewind * SlewR / 2

    if maxG > G_limit
        tRamp = G_limit / SlewR
        GM_Ramp = tRamp * G_limit # 2 ramps

        GM_box = GM - GM_Ramp
        tBox = GM_box / G_limit
        tPrewind = 2 * tRamp + tBox
    end
    if tPrewind < 0 throw("tPrewind negative") end
    return tPrewind
end

function easy_contrast(times, factor, T2s1, T2s2)
    T1 = 1900 # factor is used insteead
    return S(times, T2s1, T1) - factor * S(times, T2s2, T1)
end

## Alternative Simulation (duty cycle)

function findtimesforΔTE(Nx, NE, ME, FOV, ΔTE, field=:B7T)
    s = SCANNER[field]
    FOV /= 1e3 # in [m]
    GMsampling = Nx / s[:γ] / FOV
    config = (GMsampling, s, NE, ME, field, 2, 0)
    # find the smallest gradient, for which: (TE2 - TE1) < ΔTE
    f_mΔTE21(g) = getTEs_TR_Tacq(g, config...).TE |> TEs -> -(TEs[2] - TEs[1])
    G = binarysearch(f_mΔTE21, -ΔTE, 0, s[:Glimit])
    if -f_mΔTE21(G) ≈ ΔTE
        return getTEs_TR_Tacq(G, config...)
    else
        return Times(0,zeros(NE),0) # not possible for these settings
    end
end

# binarysearch needs monotone growing function (adjust sign of f and val)
function binarysearch(f, val, min, max)
    a, b = Float64.((min, max))
    for _ in 1:30
        mid = (a + b) / 2
        if f(mid) < val
            a = mid
        else
            b = mid
        end
    end
    b
end

function getdutycycle(t::Times)
    return t.Tacq / (t.TE[2] - t.TE[1])
end

getTacq(Tdwell, Nx) = Tdwell * Nx / 1e6