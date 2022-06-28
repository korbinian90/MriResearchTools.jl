@testset "Masking" begin
# robust mask
fn_mag = "data/small/Mag.nii"
mag = Float32.(readmag(fn_mag; rescale=true))
@test robustmask(mag) |> m -> count(.!m) / count(m) < 0.01
for i in 1:10 # test with different noise levels
    mag[(end÷2):end,:,:,:] .= i .* 0.025 .* rand.()
    m = robustmask(mag)
    @test 1.0 < count(.!m) / count(m) < 1.2
end

# phase_based_mask
fn_phase = "data/small/Phase.nii"
phase = Float32.(readphase(fn_phase))
PB = phase_based_mask(phase)
phase_based_mask(phase; filter=false)

# brain mask
brain_mask(PB)

# romeo mask
qm = romeovoxelquality(phase; mag, TEs=[4,8,12])
robustmask(qm)

end
