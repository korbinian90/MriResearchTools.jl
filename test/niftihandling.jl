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
