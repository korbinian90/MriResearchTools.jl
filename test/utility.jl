@testitem "utility" begin
# sample
sample = MriResearchTools.sample
@test length(sample(1:10)) >= 10
@test 10 >= length(sample(1:10; n=3)) >= 3
@test isempty(sample([NaN]))
@test all(isfinite.(sample([1:10..., NaN])))
@test length(sample([1])) == 1
@test isempty(sample([]))

# estimatenoise
fn_mag = "data/small/Mag.nii"
mag_nii = readmag(fn_mag; rescale=true)
@test estimatenoise(mag_nii)[2] ≈ 0.03 atol=1e-2
R = rand(500, 500, 500)
R[:, 251:500, :] .= 10
μ, σ = estimatenoise(R)
@test μ ≈ 0.5 atol=1e-1
@test σ ≈ sqrt(1/12) atol=2e-2
R[1:10,:,:] .= NaN; R[:,1:10,:] .= NaN; R[:,:,1:10] .= NaN;
R[end-9:end,:,:] .= NaN; R[:,end-9:end,:] .= NaN; R[:,:,end-9:end] .= NaN
μ, σ = estimatenoise(R)
#@test μ ≈ 0.5 atol=1e-1
@test σ ≈ sqrt(1/12) atol=1e-2

# setindex!
mag_nii[1] = 1
mag_nii[1,1,1,1] = 2
mag_nii[CartesianIndex(1,2,3,1)] = 5

# close mmapped files
GC.gc()

@test estimatequantile(1:1000, 0.8) ≈ 800 atol=1

# to_dim
@test [1 2] == to_dim([1, 2], 2)
a = 50:75
@test reshape(a, 1, 1, 1, :) == to_dim(a, 4)
@test reshape([5], 1, 1, 1, 1) == to_dim(5, 4)

end
