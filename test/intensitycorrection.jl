@testitem "intensitycorrection" begin
using Statistics

filter_nan(x) = (y = copy(x); y[.!isfinite.(y)] .= 0; y)

mag_nii = readmag("data/small/Mag.nii")
mag = Float32.(mag_nii)[:,:,:,1]

# The NIfTI convenience form derives sigma from the voxel size in the header.
corrected_nii = makehomogeneous(mag_nii)
@test size(corrected_nii) == size(mag_nii)

# The array form takes sigma in voxels.
sigma = [5, 5, 5]
corrected = makehomogeneous(mag; sigma)

@test size(corrected) == size(mag)
@test all(isfinite, corrected)
@test !all(iszero, corrected)
@test all(corrected .>= 0) # a magnitude stays a magnitude

# The point of the correction: dividing out the bias field makes the signal
# inside the object more uniform. Measure that on the same voxels, before and
# after, with the spread relative to the mean.
mask = robustmask(mag)
cv(x) = std(x) / mean(x)
@test cv(corrected[mask]) < cv(mag[mask])

# makehomogeneous is exactly "divide by getsensitivity", so the bias field it
# uses has to reproduce it.
sens = getsensitivity(mag; sigma)
@test size(sens) == size(mag)
@test all(x -> !isfinite(x) || x >= 0, sens) # a sensitivity is never negative
@test corrected ≈ mag ./ sens

# A bias field must be smooth: its gradient is far smaller than the image's.
grad(x) = mean(abs, diff(x; dims=1))
@test grad(filter_nan(sens)) < grad(mag)

# The in-place form agrees with the copying one.
inplace = copy(mag)
makehomogeneous!(inplace; sigma)
@test inplace ≈ corrected
end
