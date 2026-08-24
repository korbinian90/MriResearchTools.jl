@testitem "mcpc3ds" begin
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
complex = mag_nii .* exp.(1im .* phase_nii)
TEs = 4:4:12

combined_complex = mcpc3ds(complex; TEs=TEs)
combined_phase, combined_mag = mcpc3ds(phase_nii, mag_nii; TEs=TEs)
combined_phaseonly = mcpc3ds(phase_nii; TEs=TEs)

# The channel dimension is combined away; the echo dimension survives.
@test size(combined_complex) == size(phase_nii)
@test size(combined_phase) == size(phase_nii)
@test size(combined_phaseonly) == size(phase_nii)

@test all(isfinite, combined_complex)
@test all(isfinite, combined_phase)
@test all(isfinite, combined_phaseonly)
@test !all(iszero, combined_phase)

# Phase output is an angle.
@test all(-pi .<= combined_phase .<= pi)
@test all(-pi .<= combined_phaseonly .<= pi)
@test all(combined_mag .>= 0)

# `combinewithPO` has a separate method for complex input and for PhaseMag input.
# They are two implementations of the same sum, so they have to agree - this is
# the assertion that stops one of them drifting from the other.
angdiff(a, b) = maximum(abs.(rem2pi.(a .- b, RoundNearest)))
@test angdiff(angle.(combined_complex), combined_phase) < 1e-4

# Phase offsets are what the method exists to estimate, so they must actually be
# written to the `po` output and must not come back empty.
po = zeros(Float32, (size(phase_nii)[1:3]..., 1))
mcpc3ds(complex; TEs=TEs, po)
@test all(isfinite, po)
@test !all(iszero, po)
@test all(-pi .<= po .<= pi)

# Removing an offset that was estimated from the data is not the identity.
@test angdiff(combined_phase, Float32.(phase_nii)) > 1e-3

# MEEPI
phase_me = cat(float.(phase_nii), float.(phase_nii); dims=5)
mag_me = cat(float.(mag_nii), float.(mag_nii); dims=5)
corrected_me = mcpc3ds_meepi(phase_me, mag_me; TEs=TEs)

@test size(corrected_me) == size(phase_me)
# The two timepoints are copies of each other and one shared set of phase offsets
# is removed from both, so the corrected timepoints must stay identical.
# `isequal` rather than `==` because of the NaNs below.
@test all(isequal.(corrected_me[:,:,:,:,1], corrected_me[:,:,:,:,2]))

# Known defect, recorded rather than hidden: the MEEPI path emits NaN where the
# ordinary path does not. On this dataset that is ~2.2% of voxels, and they are
# not confined to the background - the magnitude at NaN voxels reaches 58% of the
# image maximum. The assertion above for `mcpc3ds` shows the non-MEEPI output is
# NaN-free on the same input, so this is specific to `mcpc3ds_meepi`.
# When this is fixed, this test will fail and should be turned into a plain @test.
@test_broken all(isfinite, corrected_me)
end
