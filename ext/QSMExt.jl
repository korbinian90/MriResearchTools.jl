# Backend: QSM.jl (https://github.com/kamesy/QSM.jl).
#
# NOT the recommended backend. QuantitativeSusceptibilityMappingTGV is what the
# compiled mritools binaries use and what is supported here. QSM.jl is kept
# reachable for the dipole inversions it offers that TGV does not, but upstream
# has had no release since 2023-12-03. QSM 0.5.4 calls `FFTW.libfftw3[]`, which
# throws from FFTW 1.9 on, so this extension only works if the user pins FFTW to
# 1.8 in their own environment - the package no longer caps FFTW on everyone's
# behalf to accommodate it. The README says how; test/Project.toml pins it so
# this extension stays under test. See kamesy/QSM.jl#13, open since 2026-01-08.
#
# This extension and QuantitativeSusceptibilityMappingTGVExt define the same
# `qsm_B0` and friends with identical signatures, so loading both means the one
# that loads second wins. Load exactly one per session.
module QSMExt

using MriResearchTools
using Statistics
import QSM: ismv, lbv, pdf, sharp, vsharp, nltv, rts, tikh, tkd, tsvd, tv

include("QSM_common.jl")

const γ = 267.52

function qsm(phase::AbstractArray, mask::AbstractArray, TE, vsz; bfc_mask=mask, B0=3, bfc_algo=vsharp, qsm_algo=rts, unwrapping=laplacianunwrap, bdir=(0,0,1), kw...)
    vsz = Tuple(vsz)
    uphas = unwrapping(phase)
    uphas .*= inv(B0 * γ * TE) # convert units
    fl = bfc_algo(uphas, bfc_mask, vsz) # remove non-harmonic background fields

    # some background field correction methods require a mask update
    if fl isa Tuple
        fl, mask2 = fl
        mask = mask .& mask2
    end

    # remove unsupported keywords
    kw = filter(x -> x in [:pad, :Dkernel, :bdir, :lstol, :delta, :mu, :rho, :tol, :maxit, :verbose, :tau, :gamma], keys(kw))
    x = qsm_algo(fl, mask, vsz; bdir, kw...)
    return x
end

function MriResearchTools.qsm_B0(B0_map::AbstractArray, mask::AbstractArray, vsz; bfc_mask=mask, B0=3, bfc_algo=vsharp, qsm_algo=rts, bdir=(0,0,1), kw...)
    scaled = B0_map .* (2π / (B0 * γ))
    fl = bfc_algo(scaled, bfc_mask, vsz)
    if fl isa Tuple
        fl, mask2 = fl
        mask = mask .& mask2
    end
    
    # remove unsupported keywords
    kw = filter(x -> x in [:pad, :Dkernel, :bdir, :lstol, :delta, :mu, :rho, :tol, :maxit, :verbose, :tau, :gamma], keys(kw))
    x = qsm_algo(fl, mask, vsz; bdir, kw...)
    return x
end

end
