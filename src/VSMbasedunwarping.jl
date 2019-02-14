# VSMBasedUnwarping.jl
# Methods for unwarping geometric distortions in MRI images

using Interpolations

function unwarp(VSM, distorted, dim)
    if dim == 2
        distorted = switchdim(distorted)
        VSM = switchdim(VSM)
    end

    unwarped = unwarp(VSM, distorted)

    if dim == 2
        unwarped = switchdim(unwarped)
    end
    unwarped
end

switchdim(v) = permutedims(v, [2, 1, (3:ndims(v))...])

unwarp(VSM, distorted) = unwarp!(similar(distorted), VSM, distorted)

function unwarp!(unwarped, VSM, distorted)
    xi = 1:size(distorted, 1)
    for J in CartesianIndices(size(distorted)[4:end])
        for I in CartesianIndices(size(distorted)[2:3])
            xnew = xi .+ VSM[:,I]
            xregrid = (xnew[1] .<= xi .<= xnew[end]) # only use x values inside (no extrapolation)
            unwarped[.!xregrid,I,J] .= 0
            unwarped[xregrid,I,J] .= interpolate((xnew,), distorted[:,I,J], Gridded(Linear()))(xi[xregrid]) #TODO try better interpolation than linear
        end
    end
    unwarped
end


function getVSM(B0, rbw, threshold = 5.0)
    VSM = B0 ./ (2 * pi * rbw)
    VSM .= thresholdforward(VSM, -0.9, threshold)
end

threshold(VSM, tmin, tmax) = (thresholdforward(VSM, tmin, tmax) .+ thresholdforward(VSM[end:-1:1,:,:], -tmax, -tmin))[end:-1:1,:,:] ./ 2.0

function thresholdforward(VSM, tmin, tmax)
    deltaVSM = VSM[2:end,:,:] .- VSM[1:(end-1),:,:]
    VSMret = copy(VSM)
    nx = size(VSM, 1)
    for I in CartesianIndices(size(VSM)[2:3])
        for x in 1:(nx - 1)
            if deltaVSM[x,I] < tmin || deltaVSM[x,I] > tmax
                if deltaVSM[x,I] < tmin
                    diff = tmin - deltaVSM[x,I]
                else
                    diff = tmax - deltaVSM[x,I]
                end
                VSMret[x+1,I] += diff
                if x + 1 < nx
                    deltaVSM[x+1,I] -= diff
                end
            end
        end
    end
    VSMret
end
