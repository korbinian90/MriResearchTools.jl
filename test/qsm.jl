@testset "Test QSM.jl integration" begin
using QSM
# using QSMExt
cd(@__DIR__)
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
TEs = 4:4:12

vsz = header(phase_nii).pixdim[2:4]
phase = Float32.(phase_nii)
mag = Float32.(mag_nii)
mask = qsm_mask_filled(phase[:,:,:,1])

args = (phase, mag, mask, TEs, vsz)

# QSM single-echo


# QSM multi-echo postaverage (inverse-variance-weighted averaging)
qsm_laplacian_average = qsm_average(args...; unwrapping=laplacianunwrap)
qsm_romeo_average = qsm_average(args...; unwrapping=romeo)

# QSM multi-echo phase combine
qsm_laplacian_combined = qsm_laplacian_combine(args...)
qsm_romeo_B0_map = qsm_romeo_B0(args...)

end
