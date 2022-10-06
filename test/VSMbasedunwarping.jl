@testset "VSMbasedunwarping" begin
phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
mag = Float32.(readmag(magfile))

TEs=[4,8,12]

unwrapped = romeo(phase; mag, TEs)
B0 = calculateB0_unwrapped(unwrapped, mag, TEs)

rbw = 50_000
dim = 2
vsm = getVSM(B0, rbw, dim)

unwarp(vsm, mag, dim)

end
