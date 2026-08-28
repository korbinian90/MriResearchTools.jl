using .FFTW
laplacianunwrap(ϕ) = laplacianunwrap!(copy(ϕ))
function laplacianunwrap!(ϕ::AbstractArray)
    FFTW.set_num_threads(Threads.nthreads())
    ϕ .+= 2π .* k(ϕ) # rounding k as suggested in the paper does not work
end

# Schofield and Zhu 2003, https://doi.org/10.1364/OL.28.001194
#
# The four transforms below run on arrays of the same size, so pqterm is built
# once here and passed down rather than rebuilt inside each - it is a full-volume
# array, 64 MB at 256^3.
k(ϕw) = (pq = pqterm(size(ϕw)); 1 / 2π .* ∇⁻²(∇²_nw(ϕw, pq) - ∇²(ϕw, pq), pq))  # (1)

∇²(x, pq=pqterm(size(x))) = -(2π)^ndims(x) / length(x) .* idct(pq .* dct(x))  # (2)

∇⁻²(x, pq=pqterm(size(x))) = -length(x) / (2π)^ndims(x) .* idct(dct(x) ./ pq)  # (3)

∇²_nw(ϕw, pq=pqterm(size(ϕw))) = cos.(ϕw) .* ∇²(sin.(ϕw), pq) .- sin.(ϕw) .* ∇²(cos.(ϕw), pq)  # (in text)


# Int32 holds sums of at most four squared indices exactly (4*26000^2 fits) and
# converts exactly to Float64 in the multiply and divide above, so the transforms
# are unaffected. The arithmetic itself has to stay Float64: (2) minus (3) in k is
# a difference of two nearly identical Laplacians, exactly zero wherever the phase
# is unwrapped, and the residual is amplified by length/(2π)^N before it reaches
# the phase.
pqterm(sz::NTuple{1}) = Int32[p^2 for p in 1:sz[1]]  # 1D case
pqterm(sz::NTuple{2}) = Int32[p^2 + q^2 for p in 1:sz[1], q in 1:sz[2]]  # 2D case
pqterm(sz::NTuple{3}) = Int32[p^2 + q^2 + t^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3]]  # 3D case
pqterm(sz::NTuple{4}) = Int32[p^2 + q^2 + t^2 + r^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3], r in 1:sz[4]]  # 4D case

"""
    laplacianunwrap(ϕ::AbstractArray)

Performs laplacian unwrapping on the input phase. (1D - 4D)
The phase has to be scaled to radians.
The implementation is close to the original publication: Schofield and Zhu 2003, https://doi.org/10.1364/OL.28.001194.
It is not the fastest implementation of laplacian unwrapping (doesn't use discrete laplacian).
"""
laplacianunwrap, laplacianunwrap!

# FFT variant — periodic-boundary discrete Laplacian, evaluated entirely in
# k-space. Matches the Bilgic / ICE Laplacian-unwrap implementation more
# closely than the default DCT-based version, which uses mirror boundaries.
# Works for 2D and 3D arrays. `z_weight` (3D only) rescales the through-slice
# coupling: 1 = isotropic, 0 = pure in-plane.
"""
    laplacianunwrap_fft(ϕ::AbstractArray, z_weight=1)

Periodic-boundary Laplacian phase unwrap (Schofield & Zhu, 2003) using a
discrete 5-point (2D) / 7-point (3D) Laplacian stencil evaluated by FFT.
Self-contained; no ImageFiltering dependency.
"""
function laplacianunwrap_fft(ϕ::AbstractArray, z_weight=1)
    FFTW.set_num_threads(min(4, Threads.nthreads()))
    del_op = _laplacian_kspace_kernel(size(ϕ), z_weight, eltype(ϕ))
    del_inv = 1 ./ del_op
    del_inv[.!isfinite.(del_inv)] .= 0
    cs, sn = cos.(ϕ), sin.(ϕ)
    ∇²(x) = real.(ifft(fft(x) .* del_op))
    del_phase = cs .* ∇²(sn) .- sn .* ∇²(cs)
    return real.(ifft(fft(del_phase) .* del_inv))
end

# Build the DFT of a centred discrete Laplacian stencil. For ND, the
# stencil has -2N at the centre and +1 at each of the 2N neighbours
# (with z taps scaled by z_weight in 3D).
function _laplacian_kspace_kernel(sz::NTuple{N,Int}, z_weight, T) where {N}
    kernel = zeros(float(T), sz)
    centre = ntuple(d -> (sz[d] ÷ 2) + 1, N)
    # in-plane (and z if z_weight==1) taps
    if N == 2
        kernel[centre...] = -4
        kernel[centre[1]-1, centre[2]  ] = 1
        kernel[centre[1]+1, centre[2]  ] = 1
        kernel[centre[1],   centre[2]-1] = 1
        kernel[centre[1],   centre[2]+1] = 1
    elseif N == 3
        zw = float(z_weight)
        kernel[centre...] = -(4 + 2 * zw)
        kernel[centre[1]-1, centre[2],   centre[3]  ] = 1
        kernel[centre[1]+1, centre[2],   centre[3]  ] = 1
        kernel[centre[1],   centre[2]-1, centre[3]  ] = 1
        kernel[centre[1],   centre[2]+1, centre[3]  ] = 1
        kernel[centre[1],   centre[2],   centre[3]-1] = zw
        kernel[centre[1],   centre[2],   centre[3]+1] = zw
    else
        error("laplacianunwrap_fft supports 2D or 3D arrays, got $(N)D")
    end
    return fft(ifftshift(kernel))
end
