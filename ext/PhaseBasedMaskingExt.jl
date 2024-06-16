module PhaseBasedMaskingExt

using LocalFilters
using ImageFiltering
using OffsetArrays
using MriResearchTools

"""
    phase_based_mask(phase; filter=true, threshold=1.0)

Creates a mask from a phase image.
Morphological filtering is activated by default.
To return the mask before thresholding pass `threshold=nothing`.

# Examples
```julia-repl
julia> phase_mask = phase_based_mask(phase);
```

See also [`robustmask`](@ref), [`brain_mask`](@ref)

Original MATLAB algorithm:
    se=strel('sphere',6);
    L=del2(sign(wr));
    test=convn(abs(L),se.Neighborhood,'same');
    PB=mask{1}.*(test<500);
    PB=imclose(PB,se);
    mask{2}=round(imopen(PB,se));
"""
function MriResearchTools.phase_based_mask(phase; filter=true, threshold=1.0)
    strel = sphere(6, ndims(phase))
    laplacian = imfilter(sign.(phase), Kernel.Laplacian(1:ndims(phase), ndims(phase)))
    test = imfilter(abs.(laplacian), strel)
    if isnothing(threshold)
        return test * 500 * 6
    end
    PB = test .< (500 * 6 * threshold)
    if filter
        PB = LocalFilters.closing(PB, strel)
        PB = LocalFilters.opening(PB, strel)
    end
    return PB
end

function sphere(radius, dim=3)
    len = 2radius + 1
    arr = OffsetArray(falses(repeat([len], dim)...), repeat([-radius:radius], dim)...)
    for I in CartesianIndices(arr)
        arr[I] = sqrt(sum(Tuple(I).^2)) < radius
    end
    return arr
end

end