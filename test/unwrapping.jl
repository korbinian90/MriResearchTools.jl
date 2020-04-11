phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
magni = readmag(magfile)

iswrap(x, y) = abs(rem2pi(x - y, RoundNearest)) < 1e-6

function test(f)
    unwrapped = f(phase; mag=magni)
    @test !all(unwrapped .== phase)
    @test all(iswrap.(unwrapped, phase))
    return unwrapped
end

test(romeo)
test(unwrap)
test(unwrap_individual)
@test !all(laplacianunwrap(phase) .== phase)
