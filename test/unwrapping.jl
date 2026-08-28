@testitem "unwrapping" begin
phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
magni = readmag(magfile)

iswrap(x, y) = abs(rem2pi(x - y, RoundNearest)) < 1e-6

function test(f)
    unwrapped = f(phase; mag=magni, TEs=[4,8,12])
    @test !all(unwrapped .== phase)
    @test all(iswrap.(unwrapped, phase))
    return unwrapped
end

test(romeo)
test(unwrap)
test(unwrap_individual)
@test !all(laplacianunwrap(phase) .== phase)

end

@testitem "laplacian unwrapping pqterm" begin
# pqterm's Int32 storage must not perturb the transforms: k is a difference of
# two nearly identical Laplacians, exactly zero wherever the phase is unwrapped,
# and the residual is amplified by length/(2pi)^N before it reaches the phase.
# So the assertion is bit-identical, not merely close.
pq = MriResearchTools.pqterm((6, 5, 4))
@test eltype(pq) === Int32
@test pq == [p^2 + q^2 + t^2 for p in 1:6, q in 1:5, t in 1:4]
@test eltype(MriResearchTools.pqterm((6,))) === Int32
@test eltype(MriResearchTools.pqterm((6, 5))) === Int32
@test eltype(MriResearchTools.pqterm((6, 5, 4, 3))) === Int32

# Reference k, with a plain Int64 pqterm rebuilt per transform.
ref_pq(sz) = [p^2 + q^2 + t^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3]]
ref_lap(x) = -(2pi)^ndims(x) / length(x) .* MriResearchTools.idct(ref_pq(size(x)) .* MriResearchTools.dct(x))
ref_ilap(x) = -length(x) / (2pi)^ndims(x) .* MriResearchTools.idct(MriResearchTools.dct(x) ./ ref_pq(size(x)))
ref_lapnw(p) = cos.(p) .* ref_lap(sin.(p)) .- sin.(p) .* ref_lap(cos.(p))
ref_k(p) = 1 / 2pi .* ref_ilap(ref_lapnw(p) - ref_lap(p))

phi = 2pi .* rand(12, 10, 6) .- pi
@test MriResearchTools.k(phi) == ref_k(phi)   # bit-identical, not merely close
end

