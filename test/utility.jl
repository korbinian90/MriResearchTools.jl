# Read and properly scale phase
fn_phase = "data/small/Phase.nii"
phase_nii = readphase(fn_phase)
@test maximum(phase_nii) ≈ π atol=2e-3
@test minimum(phase_nii) ≈ -π atol=2e-3

# Read and normalize mag
fn_mag = "data/small/Mag.nii"
mag_nii = readmag(fn_mag; rescale=true)
@test 1 ≤ maximum(mag_nii) ≤ 2
@test 0 ≤ minimum(mag_nii) ≤ 1

# sample
sample = MriResearchTools.sample
@test length(sample(1:10)) >= 10
@test 10 >= length(sample(1:10; n=3)) >= 3
@test isempty(sample([NaN]))
@test all(isfinite.(sample([1:10..., NaN])))
@test length(sample([1])) == 1
@test isempty(sample([]))

# estimatenoise
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

# robust mask
mag = Float32.(readmag(fn_mag; rescale=true))
@test robustmask(mag) |> m -> count(.!m) / count(m) < 0.01
for i in 1:10
    mag[(end÷2):end,:,:,:] .= i .* 0.025 .* rand.()
    m = robustmask(mag)
    @test 1.05 < count(.!m) / count(m) < 1.2
end

# savenii
fn_temp = tempname()
savenii(mag, fn_temp)
mag2 = niread(fn_temp)
@test mag == mag2

# setindex!
mag_nii[1] = 1
mag_nii[1,1,1,1] = 2
mag_nii[CartesianIndex(1,2,3,1)] = 5

# close mmapped files
GC.gc()

@test estimatequantile(1:1000, 0.8) ≈ 800 atol=1

function header_test(hdr, hdr2)
    @test hdr.scl_inter == 0
    @test hdr.scl_slope == 1
    @test hdr.dim == hdr2.dim
end
# similar
header_test(similar(mag_nii.header), mag_nii.header)
# header
header_test(header(mag_nii), mag_nii.header)
header_test(header(phase_nii), phase_nii.header)

# to_dim
@test [1 2] == to_dim([1, 2], 2)
a = 50:75
@test reshape(a, 1, 1, 1, :) == to_dim(a, 4)
@test reshape([5], 1, 1, 1, 1) == to_dim(5, 4)
