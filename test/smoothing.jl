@testitem "Smoothing" begin

phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = readphase(phasefile)[:,:,:,1]
mag = readmag(magfile)[:,:,:,1]
    
sm1 = gaussiansmooth3d(mag)
sm2 = gaussiansmooth3d(mag, [10,10,10])
sm3 = gaussiansmooth3d(mag; padding=true)
sm3 = gaussiansmooth3d(mag, [2.2,2.2,2.2]; padding=true)
sm3 = gaussiansmooth3d(mag, (2.2,2.2,2.2); padding=true)
sm3 = gaussiansmooth3d(mag, (2,2,2); padding=true)
sm4 = gaussiansmooth3d(mag; nbox=1)
sm5 = gaussiansmooth3d(mag; dims=2:3)
sm6 = gaussiansmooth3d(mag; boxsizes=[[2,1], [5,3], [4,2]], nbox=2)
sm7 = gaussiansmooth3d(mag; mask=robustmask(mag))

@test sm1 != sm2
@test sm1 != sm3
@test sm1 != sm4
@test sm1 != sm5
@test sm1 != sm6
@test sm1 != sm7

ph1 = gaussiansmooth3d_phase(phase)
ph2 = gaussiansmooth3d_phase(phase; weight=mag)
@test ph1 != ph2

end
