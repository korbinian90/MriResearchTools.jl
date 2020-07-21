# Read and properly scale phase
fn_phase = "data/small/Phase.nii"
phase_nii = readphase(fn_phase)
@test maximum(phase_nii) ≈ π atol=2e-3
@test minimum(phase_nii) ≈ -π atol=2e-3

# Read and normalize mag
fn_mag = "data/small/Mag.nii"
mag_nii = readmag(fn_mag; normalize=true)
@test 1 ≤ maximum(mag_nii) ≤ 2
@test 0 ≤ minimum(mag_nii) ≤ 1

# estimatenoise
@test 0.02 < estimatenoise(mag_nii)[2] < 0.03
R = rand(200, 200, 200)
μ, σ = estimatenoise(R)
@test μ ≈ 0.5 atol=1e-1
@test σ ≈ sqrt(1/12) atol=1e-2
R[1:10,:,:] .= NaN; R[:,1:10,:] .= NaN; R[:,:,1:10] .= NaN;
R[end-9:end,:,:] .= NaN; R[:,end-9:end,:] .= NaN; R[:,:,end-9:end] .= NaN
μ, σ = estimatenoise(R)
@test μ ≈ 0.5 atol=1e-1
@test σ ≈ sqrt(1/12) atol=1e-2

# robust mask
mag = Float32.(readmag(fn_mag; normalize=true))
for _ in 1:10
    mag[(end÷2):end,:,:,:] .= 0.2 .+ 0.1 .* rand.()
    m = robustmask(mag)
    @test 1.1 < count(.!m) / count(m) < 1.2
end

# savenii
fn_temp = tempname()
savenii(mag, fn_temp)
mag2 = niread(fn_temp)
@test mag == mag2

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
