using .FFTW
laplacianunwrap(ϕ) = laplacianunwrap!(copy(ϕ))
function laplacianunwrap!(ϕ::AbstractArray)
    FFTW.set_num_threads(Threads.nthreads())
    ϕ .+= 2π .* k(ϕ) # rounding k as suggested in the paper does not work
end

# Schofield and Zhu 2003, https://doi.org/10.1364/OL.28.001194
k(ϕw) = 1 / 2π .* ∇⁻²(∇²_nw(ϕw) - ∇²(ϕw))  # (1)

∇²(x) = -(2π)^ndims(x) / length(x) .* idct(pqterm(size(x)) .* dct(x))  # (2)

∇⁻²(x) = -length(x) / (2π)^ndims(x) .* idct(dct(x) ./ pqterm(size(x)))  # (3)

∇²_nw(ϕw) = cos.(ϕw) .* ∇²(sin.(ϕw)) .- sin.(ϕw) .* ∇²(cos.(ϕw))  # (in text)


pqterm(sz::NTuple{1}) = (1:sz[1]).^2  # 1D case
pqterm(sz::NTuple{2}) = [p^2 + q^2 for p in 1:sz[1], q in 1:sz[2]]  # 2D case
pqterm(sz::NTuple{3}) = [p^2 + q^2 + t^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3]]  # 3D case
pqterm(sz::NTuple{4}) = [p^2 + q^2 + t^2 + r^2 for p in 1:sz[1], q in 1:sz[2], t in 1:sz[3], r in 1:sz[4]]  # 4D case

"""
    laplacianunwrap(ϕ::AbstractArray)

Performs laplacian unwrapping on the input phase. (1D - 4D)
The phase has to be scaled to radians.
The implementation is close to the original publication: Schofield and Zhu 2003, https://doi.org/10.1364/OL.28.001194.
It is not the fastest implementation of laplacian unwrapping (doesn't use discrete laplacian).
"""
laplacianunwrap, laplacianunwrap!

# FFT variant
using OffsetArrays, ImageFiltering
function laplacianunwrap_fft(ϕ::AbstractArray, z_weight=1)
    FFTW.set_num_threads(min(4, Threads.nthreads()))

    kernel = float.(convert(AbstractArray, Kernel.Laplacian((true,true,true))))
    kernel[0,0,1] = kernel[0,0,-1] = z_weight
    kernel[0,0,0] += 2 * (1 - z_weight)

    ∇²(x) = imfilter(x, kernel)


    kernel_full = centered(zeros(size(ϕ)))
    kernel_full[CartesianIndices(kernel)] .= kernel
        
    del_op = fft(kernel_full)
    
    del_inv = 1 ./ del_op
    del_inv[.!isfinite.(del_inv)] .= 0
    
    del_phase = cos.(ϕ) .* ∇²(sin.(ϕ)) .- sin.(ϕ) .* ∇²(cos.(ϕ))
    
    unwrapped = real.(ifft( fft(del_phase) .* del_inv ))
    
    return unwrapped
end

function laplacianunwrap_mixed(ϕ::AbstractArray)
    FFTW.set_num_threads(min(4, Threads.nthreads()))

    kernel = Kernel.Laplacian((true,true,true))
    ∇²(x) = imfilter(x, kernel)
    
    del_phase = cos.(ϕ) .* ∇²(sin.(ϕ)) .- sin.(ϕ) .* ∇²(cos.(ϕ))
    
    unwrapped = ∇⁻²(del_phase)
    
    return unwrapped
end
