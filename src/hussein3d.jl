const path = "$(@__DIR__)/hussein_3d.so"

hussein3d(wrapped::AbstractArray{Float64, 3}) = hussein3d!(similar(wrapped), wrapped)

function hussein3d!(unwrapped::AbstractArray{Float64, 3}, wrapped::AbstractArray{Float64, 3})
    maskUInt8 = zeros(UInt8, size(wrapped))
    ccall((:unwrap3D, path), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}, Cint, Cint, Cint, Cint, Cint, Cint, UInt8, UInt32), wrapped, unwrapped, maskUInt8, size(wrapped, 1), size(wrapped, 2), size(wrapped, 3), 0, 0, 0, 0, 0)
    unwrapped
end
