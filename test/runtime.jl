#fl = [MriResearchTools.gaussiansmooth3dpar!, MriResearchTools.gaussiansmooth3d!]
fl = [MriResearchTools.gaussiansmooth3d!]
a = rand(100,100,100)
w = rand(100,100,100)
m = ones(100,100,100)
a1 = gaussiansmooth3d(a; weight=copy(w), mask=m)
for f in fl
    a2 = copy(a)
    f(a2; weight=copy(w), mask=m)
    @show count((a1 .≈ a2) .== false)
    @show length(a1)
end

function testruntime(f, dtype=Float32)
    a = rand(dtype, 100,100,100)
    w = rand(dtype, 100,100,100)
    m = rand(Bool, 100,100,100)
    for sz in 2:25
        f(a, sz; mask=m)
    end
end

for f in fl
    testruntime(f)
end

for dtype in [Float32], f in fl
    println("$f: ")
    @time testruntime(f, dtype)
end