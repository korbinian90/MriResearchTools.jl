@testset "mcpc3ds" begin
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
complex = mag_nii .* exp.(1im .* phase_nii)
TEs = 4:4:12

mcpc3ds(complex; TEs=TEs)

end
