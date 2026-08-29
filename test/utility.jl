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

@testitem "getHIP element type" begin
# The accumulator follows the data: it is the largest transient in mcpc3ds, and
# the result is read only through abs() and angle().
mag = Float32.(4000 .* rand(8, 8, 4, 2, 8))
phase = Float32.(2pi .* rand(8, 8, 4, 2, 8) .- pi)
hip32 = getHIP(mag, phase)
@test eltype(hip32) === ComplexF32

# Float64 input is preserved: this follows the input, it does not force Float32.
@test eltype(getHIP(Float64.(mag), Float64.(phase))) === ComplexF64
# An integer magnitude must not produce a Complex{Int} accumulator that cis cannot
# be accumulated into.
@test eltype(getHIP(round.(UInt16, mag), phase)) === ComplexF32

# Agreement with the Float64 computation, on the two quantities anyone reads.
hip64 = getHIP(Float64.(mag), Float64.(phase))
@test maximum(abs.(abs.(hip64) .- abs.(hip32)) ./ abs.(hip64)) < 1e-5
@test maximum(abs.(mod.(angle.(hip64) .- angle.(hip32) .+ pi, 2pi) .- pi)) < 1e-5
end


@testitem "Aqua" begin
    using Aqua
    # piracies is off, and it is a real finding rather than a false positive:
    # copy, similar and setindex! are defined here on NIfTI's NIfTI1Header and
    # NIVolume (niftihandling.jl:79, :81, :204). They are long-standing
    # conveniences this package deliberately owns, so the check is disabled
    # rather than suppressed one by one - revisit if they ever move upstream.
    Aqua.test_all(MriResearchTools; piracies=false)
end
