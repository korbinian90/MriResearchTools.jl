@testitem "Test QSM integration" begin
using Statistics

# using QSM
using QSM

cd(@__DIR__)
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
TEs = 4:4:12

vsz = header(phase_nii).pixdim[2:4] .* 10 # reduces required iterations for testing
phase = Float32.(phase_nii)
mag = Float32.(mag_nii)
mask = qsm_mask_filled(phase[:,:,:,1])
B0 = 3

args = (phase, mag, mask, TEs, vsz)

# QSM single-echo
qsm_single = qsm_romeo_B0(phase[:,:,:,1], mag[:,:,:,1], mask, TEs[1], vsz; B0, iterations=5)

# QSM multi-echo postaverage (inverse-variance-weighted averaging)
qsm_laplacian_average = qsm_average(args...; B0, iterations=5)
# QSM.jl
# qsm_laplacian_average = qsm_average(args...; B0, iterations=5, unwrapping=laplacianunwrap)
# qsm_romeo_average = qsm_average(args...; B0, iterations=5, unwrapping=romeo)

# QSM multi-echo phase combine
qsm_laplacian_combined = qsm_laplacian_combine(args...; B0, iterations=5)
qsm_romeo_B0_map = qsm_romeo_B0(args...; B0, iterations=5)
qsm_romeo_B0_map_automask = qsm_romeo_B0(phase, mag, nothing, TEs, vsz; B0, iterations=5)

# ---- assertions on the results above -------------------------------------
sz = size(phase)[1:3]

for (name, chi) in ("qsm_average" => qsm_laplacian_average,
                    "qsm_laplacian_combine" => qsm_laplacian_combined,
                    "qsm_romeo_B0" => qsm_romeo_B0_map,
                    "qsm_romeo_B0 automask" => qsm_romeo_B0_map_automask)
    @testset "$name" begin
        @test size(chi) == sz
        @test all(isfinite, chi)
        @test !all(iszero, chi)
        # A susceptibility map is in ppm and brain tissue sits within a few ppm of
        # water. Anything far outside that is a unit or scaling error, which is the
        # failure mode these calls are most likely to regress into.
        @test maximum(abs, chi) < 10
        # The reconstruction is driven by the mask it was given, so there must be
        # substantially less signal outside it than inside. (Not `iszero` outside:
        # TGV zeroes the region outside its bounding box, QSM.jl does not, and
        # this test file runs against whichever backend is loaded.)
        @test mean(abs, chi[.!mask]) < mean(abs, chi[mask])
    end
end

# The single-echo call takes one echo and must still produce a full-size map.
@test size(qsm_single) == sz
@test all(isfinite, qsm_single)

# The multi-echo routes are different estimators of the same physical quantity on
# the same data, so they have to agree in broad shape even at these low iteration
# counts. Correlate inside the mask; the background is zero by construction and
# would inflate any whole-volume number.
cor_in_mask(a, b) = cor(vec(a[mask]), vec(b[mask]))
@test cor_in_mask(qsm_laplacian_combined, qsm_romeo_B0_map) > 0.5
@test cor_in_mask(qsm_laplacian_average, qsm_laplacian_combined) > 0.5

# Deliberately not asserted: the correlation between the supplied-mask and
# automask runs. Passing no mask makes the routine build its own from the B0 map
# rather than from echo-1 phase, which is a genuinely different mask, and QSM is
# strongly mask-sensitive - measured r = 0.37 on the TGV backend and 0.15 on
# QSM.jl. Any bound tight enough to be meaningful would just be asserting which
# backend is loaded. The automask result is checked for the properties above
# instead.

end
