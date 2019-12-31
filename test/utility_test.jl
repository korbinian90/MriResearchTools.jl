# Read and properly scale phase
fn_phase = "data/small/Phase.nii"
phase_nii = readphase(fn_phase)
@test maximum(phase_nii) ≈ π
@test minimum(phase_nii) ≈ -π

# Read and normalize mag
fn_mag = "data/small/Mag.nii"
mag_nii = readmag(fn_mag; normalize=true)
@test 1 ≤ maximum(mag_nii) ≤ 2
@test 0 ≤ minimum(mag_nii) ≤ 1

#TODO savenii
#TODO createniiforwriting

@test estimatequantile(1:1000, 0.8) ≈ 800 atol=1

hdr = similar(mag_nii.header)
@test hdr.scl_inter == 0
@test hdr.scl_slope == 1
@test hdr.dim == mag_nii.header.dim
