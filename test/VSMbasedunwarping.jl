@testitem "VSMbasedunwarping" begin
phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
mag = Float32.(readmag(magfile))

TEs=[4,8,12]

unwrapped = romeo(phase; mag, TEs)
B0 = calculateB0_unwrapped(unwrapped, mag, TEs)

rbw = 50_000
dim = 2

vsm = getVSM(B0, rbw, dim)

@test size(vsm) == size(B0)
@test all(isfinite, vsm)
# A voxel shift map is B0/rbw plus a forward threshold, so no field means no shift.
@test getVSM(zeros(size(B0)), rbw, dim) == zeros(size(B0))
unwarped = unwarp(vsm, mag, dim)

@test size(unwarped) == size(mag)
@test all(isfinite, unwarped)
@test !all(iszero, unwarped)

# The strong invariant: a zero shift map has to be the identity. Every sample
# lands back on its own grid point, so the interpolation is exact and no voxel
# falls outside the regrid range.
@test unwarp(zeros(Float32, size(B0)), mag, dim) ≈ mag
# ... in either readout direction.
@test unwarp(zeros(Float32, size(B0)), mag, 1) ≈ mag

# thresholdforward is what stops the map from folding over itself: after it, no
# forward difference along the readout may exceed the threshold it was given.
tmax = 5.0
v = getVSM(B0, rbw, 1, tmax)
@test maximum(diff(v; dims=1)) <= tmax + 1e-4

end
