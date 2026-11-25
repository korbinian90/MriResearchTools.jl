@testitem "methods" begin
# Data
fn_phase = "data/small/Phase.nii"
fn_mag = "data/small/Mag.nii"
phase_nii = readphase(fn_phase)
mag_nii = readmag(fn_mag)
I = mag_nii .* exp.(1im .* phase_nii)
TEs = 4:4:12

# homodyne
h1 = homodyne(mag_nii, phase_nii)
h2 = homodyne(Float32.(mag_nii), Float32.(phase_nii))
h3 = homodyne(I)
h4 = homodyne(mag_nii, phase_nii; sigma=[5,5])
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

# calculateB0
unwrapped = romeo(phase_nii; TEs=TEs)
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs)
@test all(isfinite.(B0))
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :phase_snr)
@test all(isfinite.(B0))
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :average)
@test all(isfinite.(B0))
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :TEs)
@test all(isfinite.(B0))
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :mag)
@test all(isfinite.(B0))
B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :simulated_mag)
@test all(isfinite.(B0))
snr = get_B0_snr(mag_nii, TEs)
@test all(isfinite.(snr))
m = "The phase weighting option 'magTEs' is not defined!"
@test_throws ErrorException(m) B0 = calculateB0_unwrapped(unwrapped, mag_nii, TEs, :magTEs)

# make border of image noise
phase = Float32.(phase_nii)
phase[1:3,:,:,:] .= 2π .* rand.() .- π
phase[(end-2):end,:,:,:] .= 2π .* rand.() .- π
phase[:,1:3,:,:] .= 2π .* rand.() .- π
phase[:,(end-2):end,:,:] .= 2π .* rand.() .- π
phase[:,:,1:3,:] .= 2π .* rand.() .- π
phase[:,:,(end-2):end,:] .= 2π .* rand.() .- π

# getvoxelquality
vq = romeovoxelquality(phase; TEs=TEs)
@test all(isfinite.(vq))
@test size(vq) == size(phase_nii)[1:3]

# mask_from_voxelquality
mask = mask_from_voxelquality(vq; threshold=0.5)
@test mask isa AbstractArray{<:Bool}
@test !all(mask .== true)
@test !all(mask .== false)

bm = brain_mask(mask)
@test bm isa AbstractArray{<:Bool}
@test !all(bm .== true)
@test !all(bm .== false)


# NumART2star
t2s_true = 30.0
decay(TEs) = reshape(exp.(-TEs ./ t2s_true), 1, 1, 1, :)
TEs_eq = 4:4:20
m = decay(TEs_eq)
old = (TEs_eq[end]-TEs_eq[1]) / 2(length(TEs_eq)-1) * (m[1]+m[end]+2sum(m[2:end-1])) / (m[1]-m[end])
@test NumART2star(m, TEs_eq)[1] ≈ old
TEs_uneq = [4, 5, 7, 12, 20]
m = decay(TEs_uneq)
trapz = sum(diff(TEs_uneq) .* (m[1:end-1] .+ m[2:end]) ./ 2)
@test NumART2star(m, TEs_uneq)[1] ≈ trapz / (m[1]-m[end])
@test abs(NumART2star(decay(4:0.1:20), 4:0.1:20)[1] - t2s_true) < 0.01

end
