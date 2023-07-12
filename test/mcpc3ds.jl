@testset "mcpc3ds" begin
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
complex = mag_nii .* exp.(1im .* phase_nii)
TEs = 4:4:12

mcpc3ds(complex; TEs=TEs)
mcpc3ds(phase_nii, mag_nii; TEs=TEs)
mcpc3ds(phase_nii; TEs=TEs)

# MEEPI
phase_me = cat(float.(phase_nii), float.(phase_nii); dims=5)
mag_me = cat(float.(mag_nii), float.(mag_nii); dims=5)
mcpc3ds_meepi(phase_me, mag_me; TEs=TEs)
end
