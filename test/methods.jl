@testset "methods" begin
# homodyne
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
I = mag_nii .* exp.(1im .* phase_nii)

h1 = homodyne(mag_nii, phase_nii)
h2 = homodyne(Float32.(mag_nii), Float32.(phase_nii))
h3 = homodyne(I)
h4 = homodyne(mag_nii, phase_nii; σ=5)
h5 = homodyne(mag_nii, phase_nii; dims=1:3)

@test h1 == h2
@test h1 == h3
@test h1 != h4
@test h1 != h5
@test h4 != h5

# inplace test
I2 = copy(I)
@test I == I2
homodyne!(I2)
@test I != I2
@test h3 == I2

end
