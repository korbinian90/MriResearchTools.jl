@testset "Masking" begin
using ImageFiltering
# robust mask
fn_mag = "data/small/Mag.nii"
mag = Float32.(readmag(fn_mag; rescale=true))
@test robustmask(mag) |> m -> count(.!m) / count(m) < 0.01
for i in 1:8 # test with different noise levels
    mag[(end÷2):end,:,:,:] .= i .* 0.025 .* rand.()
    m = robustmask(mag)
    @test 0.9 < count(.!m) / count(m) < 1.1
end

# phase_based_mask
fn_phase = "data/small/Phase.nii"
phase = Float32.(readphase(fn_phase))
PB = phase_based_mask(phase)
phase_based_mask(phase; filter=false)
qm = phase_based_mask(phase; filter=false, threshold=nothing)
robustmask(qm)

# brain mask
brain_mask(PB)

# romeo mask
qm = romeovoxelquality(phase; mag, TEs=[4,8,12])
robustmask(qm)

end
