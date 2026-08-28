@testitem "niftihandling" begin
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

fn_mag_gz = "data/small/Mag.nii.gz"
@test all(readmag(fn_mag_gz) .== readmag(fn_mag))

# Test int16 rescale
fn_int16 = "data/small/int16.nii"
int16_nii = readmag(fn_int16)
int16_nii[:]
int16p_nii = readphase(fn_int16)
int16p_nii[:]

# Test GE fix for crash
readphase(fn_phase; fix_ge=true)
readphase(fn_int16; fix_ge=true)

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

# savenii
fn_temp = tempname()
mag = Float32.(mag_nii)
savenii(mag, fn_temp)
mag2 = niread(fn_temp)
@test mag == mag2

dir_temp = tempdir()
savenii(mag, "name", dir_temp)
@test isfile(joinpath(dir_temp, "name.nii"))

dir_temp = tempdir()
savenii(mag, "name2.nii", dir_temp)
@test isfile(joinpath(dir_temp, "name2.nii"))

dir_temp = tempdir()
savenii(mag, "name3.nii.gz", dir_temp)
@test isfile(joinpath(dir_temp, "name3.nii.gz"))

@test filesize(joinpath(dir_temp, "name2.nii")) != filesize(joinpath(dir_temp, "name3.nii.gz")) > 0

rm.(joinpath.(dir_temp, ["name.nii", "name2.nii", "name3.nii.gz"]))

end

@testitem "savenii output type" begin
using NIfTI
# What savenii writes is decided by the eltype of the array it hands to niwrite -
# a header's own datatype is discarded - so this is the one place the output type
# is controlled, for every package that writes through it.
d = mktempdir()
ondisk(img; kw...) = niread(savenii(img, joinpath(d, "t.nii"); kw...))

# Float64 results are written as Float32: half the file, no less information.
v = ondisk(rand(Float64, 6, 6, 3))
@test eltype(v.raw) === Float32
@test v.header.bitpix == 32
# ... and the values are exactly the narrowed originals, not something rescaled.
a = rand(Float64, 6, 6, 3)
@test niread(savenii(a, joinpath(d, "vals.nii"))).raw == Float32.(a)

@test eltype(ondisk(rand(Float32, 6, 6, 3)).raw) === Float32   # already right, untouched
@test eltype(ondisk(rand(ComplexF64, 6, 6, 3)).raw) === ComplexF32

# Masks are promoted rather than written as UInt8, because NIfTI.jl 0.6.2 writes
# bitpix 32 for every integer width other than 32 bits - a UInt8 mask would get a
# header contradicting itself. The values must still be exactly 0 and 1.
m = rand(Bool, 6, 6, 3)
w = ondisk(m)
@test eltype(w.raw) === Float32
@test w.header.bitpix == 32
@test (w.raw .> 0.5) == m
@test all(x -> x == 0 || x == 1, w.raw)
@test eltype(ondisk(rand(UInt8, 6, 6, 3)).raw) === Float32

# Both escape hatches.
v = ondisk(rand(Float64, 6, 6, 3); datatype=nothing)
@test eltype(v.raw) === Float64 && v.header.bitpix == 64
v = ondisk(Float64.(rand(0:4, 6, 6, 3)); datatype=Int32)
@test eltype(v.raw) === Int32

# No copy when the array is already the target type.
x = rand(Float32, 4, 4, 2)
@test MriResearchTools.to_output_type(x, Float32) === x
end

