function phase_based_mask(phase)
    #laplacian = ∇²(sign.(phase))
    laplacian = imfilter(sign.(phase), Kernel.Laplacian())
    test = imfilter(abs.(laplacian), sphere(6, 3))
    PB = test < 500
    return PB

# se=strel('sphere',6);
# L=del2(sign(wr));
# test=convn(abs(L),se.Neighborhood,'same');
# PB=mask{1}.*(test<500);
# PB=imclose(PB,se);
# mask{2}=round(imopen(PB,se));
end

function sphere(radius, dim)
    len = 2radius + 1
    arr = OffsetArray(zeros(repeat([len], dim)...), repeat([-r:r], dim)...)
    for I in CartesianIndices(arr)
        arr[I] = sqrt(sum(Tuple(I).^2)) < radius
    end
    return arr
end

## Test
pth = raw"C:\MRI\Gisela data\v18"

phase = readphase(joinpath(pth, "aspire_phase.nii"))[:,:,:,1];
laplacian = imfilter(sign.(phase), Kernel.Laplacian(), NoPad(Pad(:replicate)))

laplacian = ∇²(sign.(phase))