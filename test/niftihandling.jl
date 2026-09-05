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
# niread is exported by MriResearchTools, so this needs no NIfTI test dependency.
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

# Masks go out as UInt8, a quarter of the bytes. datatype and bitpix have to
# agree, which is what the NIfTI 0.6.3 lower bound is for: readers outside Julia
# trust bitpix, and before that fix this array was written as datatype 2 with
# bitpix 32. The values must still be exactly 0 and 1.
m = rand(Bool, 6, 6, 3)
w = ondisk(m)
@test eltype(w.raw) === UInt8
@test w.header.datatype == 2
@test w.header.bitpix == 8
@test (w.raw .> 0) == m
@test all(x -> x == 0 || x == 1, w.raw)
# Only Bool is a mask. A UInt8 array carries values, so it keeps the float path.
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


@testitem "static NIfTI" begin
    using MriResearchTools.NIfTI
    fn_phase = "data/small/Phase.nii"
    fn_mag = "data/small/Mag.nii"
    fn_int16 = "data/small/int16.nii"

    # the loaders return what readphase, readmag and niread return, as Float32 with five dimensions
    for (fn, read, load) in ((fn_phase, readphase, loadphase), (fn_mag, readmag, loadmag), (fn_phase, niread, loadnii),
                             (fn_int16, readphase, loadphase), (fn_int16, readmag, loadmag), ("data/small/Mag.nii.gz", readmag, loadmag))
        vol = read(fn)
        a, hdr = load(fn)
        @test a isa Array{Float32,5}
        @test size(a) == (size(vol)..., ntuple(_ -> 1, 5 - ndims(vol))...)
        @test all(Float32.(vol) .== reshape(a, size(vol)))
        @test hdr.dim == vol.header.dim
        @test hdr.scl_slope == 1 && hdr.scl_inter == 0 # the scaling is in the data
    end
    a, _ = loadmag(fn_int16; rescale=true)
    @test all(Float32.(readmag(fn_int16; rescale=true)) .== reshape(a, size(a)[1:4]))
    a, _ = loadphase(fn_int16; fix_ge=true)
    @test all(Float32.(readphase(fn_int16; fix_ge=true)) .== reshape(a, size(a)[1:4]))
    @test loadheader(fn_phase).dim == readphase(fn_phase).header.dim

    # the writer produces the bytes NIfTI.jl produces
    vol = readphase(fn_phase)
    hdr = header(vol)
    img = Float32.(vol)[:,:,:,1]
    tmp = mktempdir()
    for name in ("a.nii", "a.nii.gz")
        savenii(img, name, tmp, hdr)
        niwrite(joinpath(tmp, "b" * name[2:end]), NIVolume(hdr, img))
        @test read(joinpath(tmp, name)) == read(joinpath(tmp, "b" * name[2:end]))
        @test niread(joinpath(tmp, name)) == img
    end
    savenii(img .> 100, "mask", tmp, hdr)
    m, _ = loadnii(joinpath(tmp, "mask.nii"))
    @test eltype(niread(joinpath(tmp, "mask.nii")).raw) == UInt8
    @test all((m .!= 0)[:,:,:,1,1] .== (img .> 100))
    savenii(round.(img .* 100), "int16", tmp, hdr; datatype=Int16)
    @test eltype(niread(joinpath(tmp, "int16.nii")).raw) == Int16
    @test_throws ArgumentError savenii(img, "conflict", tmp, hdr; voxel_size=(1, 1, 1))

    # a big endian file reads the same
    bytes = read(fn_phase)
    be = copy(bytes)
    io = IOBuffer(); write(io, NIfTI.byteswap(copy(vol.header)))
    be[1:348] .= take!(io)
    be[353:end] .= reinterpret(UInt8, hton.(vec(collect(vol.raw))))
    write(joinpath(tmp, "be.nii"), be)
    a_be, hdr_be = loadphase(joinpath(tmp, "be.nii"))
    @test a_be == first(loadphase(fn_phase))
    @test hdr_be.dim == vol.header.dim
    @test all(Float32.(readphase(joinpath(tmp, "be.nii"))) .== reshape(a_be, size(vol)))
end
