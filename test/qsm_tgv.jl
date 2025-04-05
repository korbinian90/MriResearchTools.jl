@testitem "Test QuantitativeSusceptibilityMappingTGV integration" begin

# using QSM
using QuantitativeSusceptibilityMappingTGV

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


# QSM multi-echo postaverage (inverse-variance-weighted averaging)
qsm_laplacian_average = qsm_average(args...; B0, iterations=5)
# QSM.jl
# qsm_laplacian_average = qsm_average(args...; B0, iterations=5, unwrapping=laplacianunwrap)
# qsm_romeo_average = qsm_average(args...; B0, iterations=5, unwrapping=romeo)

# QSM multi-echo phase combine
qsm_laplacian_combined = qsm_laplacian_combine(args...; B0, iterations=5)
qsm_romeo_B0_map = qsm_romeo_B0(args...; B0, iterations=5)
qsm_romeo_B0_map_automask = qsm_romeo_B0(phase, mag, nothing, TEs, vsz; B0, iterations=5)

end
