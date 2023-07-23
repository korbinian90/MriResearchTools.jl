function rotate3D(image::AbstractArray{T,N}, angle_x, mode=Linear()) where {T,N}
    centered = OffsetArrays.centered(image)
    R = RotX(deg2rad(angle_x))
    rotated_points = [R * p for p in coordinates(image)]
    return resample_new_points(centered, rotated_points, mode)
end

function affine_transformation(image::AbstractArray{T,N}, affine, mode=Linear()) where {T,N}
    centered = OffsetArrays.centered(image)
    aff_map = AffineMap(affine[1:3,1:3], affine[1:3,4])
    new_points = map(aff_map, coordinates(image))
    return resample_new_points(centered, new_points, mode)
end

coordinates(image) = SVector.(Tuple.(CartesianIndices(image)))

function resample_new_points(array::AbstractArray{T,N}, new_points, mode) where {T,N}
    itp = interpolate(array, BSpline(mode))
    safe_apply(r) =
    if checkbounds(Bool, array, r...)
        itp(r...)
    else
        zero(T)
    end
    return T.(safe_apply.(new_points))
end
