# References for the methods this package implements.
#
# The writer itself lives in ROMEO, which has no dependencies, so every package
# can reach it without a version constraint pointing back up the graph. What
# lives here is only what MriResearchTools owns: register the reference next to
# the code that implements the method, and the arrows all point one way.

function __init__()
    # MCPC-3D-S and ASPIRE were published together in the ASPIRE paper, so they
    # share a reference but are two different methods, and only ASPIRE is
    # patented. This package implements MCPC-3D-S: mcpc3ds unwraps the HIP with
    # ROMEO on every path and never takes the ASPIRE shortcut (see mcpc3ds.jl).
    # So :mcpc3ds is what the code cites, and it carries no patent notice.
    register_citation!(:mcpc3ds,
        """Eckstein, K., Dymerska, B., Bachrata, B., Bogner, W., Poljanc, K., Trattnig, S., Robinson, S.D., 2018.
           Computationally Efficient Combination of Multi-channel Phase Data From Multi-echo Acquisitions (ASPIRE).
           Magnetic Resonance in Medicine 79, 2996-3006.
           https://doi.org/10.1002/mrm.26963""";
        label = "MCPC-3D-S Coil Combination")

    # Registered for the method itself, which nothing here runs today: ASPIRE
    # skips unwrapping when the echo times satisfy TE2 = n*TE1. If that path is
    # ever implemented, cite :aspire and the notice comes with it.
    register_citation!(:aspire,
        """Eckstein, K., Dymerska, B., Bachrata, B., Bogner, W., Poljanc, K., Trattnig, S., Robinson, S.D., 2018.
           Computationally Efficient Combination of Multi-channel Phase Data From Multi-echo Acquisitions (ASPIRE).
           Magnetic Resonance in Medicine 79, 2996-3006.
           https://doi.org/10.1002/mrm.26963""";
        label = "ASPIRE Coil Combination",
        notice =
        """PATENT: ASPIRE is covered by US10605885B2
           (https://patents.google.com/patent/US10605885B2/en). Per the upstream ASPIRE
           repository, no licence is required for scientific use and the method can be
           applied free of charge, but a licence IS required for commercial use, and the
           method is not a medical product, so it may not be used for diagnosis in humans.
           Note that an MIT licence grants copyright permissions only, not patent rights.""")

    register_citation!(:homogeneity,
        """Eckstein, K., Trattnig, S., Robinson, S.D., 2019.
           A Simple Homogeneity Correction for Neuroimaging at 7T.
           Proceedings of the 27th Annual Meeting ISMRM. Presented at the ISMRM, Montreal, Quebec, Canada.
           https://index.mirasmart.com/ISMRM2019/PDFfiles/2716.html""";
        label = "Homogeneity Correction")

    register_citation!(:laplacian,
        """Schofield, M.A., Zhu, Y., 2003.
           Fast phase unwrapping algorithm for interferometric applications.
           Optics Letters 28, 1194-1196.
           https://doi.org/10.1364/OL.28.001194""";
        label = "Laplacian Unwrapping")

    register_citation!(:rts,
        """Kames, C., Wiggermann, V., Rauscher, A., 2018.
           Rapid two-step dipole inversion for susceptibility mapping with sparsity priors.
           NeuroImage 167, 276-286.
           https://doi.org/10.1016/j.neuroimage.2017.11.018""";
        label = "RTS Dipole Inversion")

    register_citation!(:phase_based_masking,
        """Hagberg, G.E., Eckstein, K., Tuzzi, E., Zhou, J., Robinson, S.D., Scheffler, K., 2022.
           Phase-based masking for quantitative susceptibility mapping of the human brain at 9.4T.
           Magnetic Resonance in Medicine.
           https://doi.org/10.1002/mrm.29368""";
        label = "Phase-based Masking")

    register_citation!(:qsmxt,
        """Stewart, A.W., Robinson, S.D., O'Brien, K., Jin, J., Widhalm, G., Hangel, G., Walls, A., Goodwin, J., Eckstein, K., Tourell, M., Morgan, C., Narayanan, A., Barth, M., Bollmann, S., 2022.
           QSMxT: Robust masking and artifact reduction for quantitative susceptibility mapping.
           Magnetic Resonance in Medicine.
           https://doi.org/10.1002/mrm.29048""";
        label = "QSMxT Masking")
end

"""
    describe_input(path)

Render an input file for a provenance record: the absolute path plus its
dimensions. Pass as the `describe` keyword of `write_provenance`.
"""
function describe_input(path)
    p = try abspath(String(path)) catch; String(path) end
    isfile(p) || return "$p (not found)"
    dims = try string(size(niread(p))) catch; "unreadable" end
    return "$p  $dims"
end
