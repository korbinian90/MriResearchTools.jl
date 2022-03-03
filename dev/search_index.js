var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MriResearchTools","category":"page"},{"location":"#MriResearchTools","page":"Home","title":"MriResearchTools","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MriResearchTools.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MriResearchTools]","category":"page"},{"location":"#MriResearchTools.NumART2star-Union{Tuple{T}, Tuple{AbstractArray{T, 4}, Any}} where T","page":"Home","title":"MriResearchTools.NumART2star","text":"NumART2star(image::AbstractArray{T,4}, TEs) where T\n\nPerforms T2* calculation on 4D-multi-echo magnitude data. https://doi.org/10.1002/mrm.10283\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.RSS-Tuple{Any}","page":"Home","title":"MriResearchTools.RSS","text":"RSS(mag; dims=ndims(mag))\n\nPerforms root-sum-of-squares combination along the last dimension of mag. The dimension can be specificed via the dims keyword argument.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.calculateB0_unwrapped-Tuple{Any, Any, Any}","page":"Home","title":"MriResearchTools.calculateB0_unwrapped","text":"calculateB0_unwrapped(unwrapped_phase, mag, TEs)\n\nCalculates B0 in [Hz] from unwrapped phase. The phase offsets have to be removed prior.\n\nSee also mcpc3ds and romeo\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.estimatenoise-Tuple{AbstractArray}","page":"Home","title":"MriResearchTools.estimatenoise","text":"estimatenoise(image::AbstractArray)\n\nEstimates the noise from the corners of the image. It assumes that at least one corner is without signal and only contains noise.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.estimatequantile-Tuple{Any, Any}","page":"Home","title":"MriResearchTools.estimatequantile","text":"estimatequantile(array, p)\n\nQuickly estimates the quantile p of a possibly large array by using a subset of the data.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.gaussiansmooth3d","page":"Home","title":"MriResearchTools.gaussiansmooth3d","text":"gaussiansmooth3d(image)\n\ngaussiansmooth3d(image, σ=[5,5,5];\n    mask=nothing,\n    nbox=ifelse(isnothing(mask), 3, 6), \n    weight=nothing, dims=1:min(ndims(image),3), \n    boxsizes=getboxsizes.(σ, nbox)\n    )\n\nPerforms Gaussian smoothing on image with σ as standard deviation of the Gaussian. By application of nbox times running average filters in each dimension. The length of σ and the length of the dims that are smoothed have to match. (Default 3)\n\nOptional arguments:\n\nmask: Smoothing can be performed using a mask to inter-/extrapolate missing values.\nnbox: Number of box applications. Default is 3 for normal smoothing and 6 for masked smoothing.\nweight: Apply weighted smoothing. Either weighted or masked smoothing can be porformed.\ndims: Specify which dims should be smoothed. Corresponds to manually looping of the other dimensions.\nboxizes: Manually specify the boxsizes, not using the provided σ. length(boxsizes)==length(dims) && length(boxsizes[1])==nbox\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.gaussiansmooth3d!","page":"Home","title":"MriResearchTools.gaussiansmooth3d!","text":"gaussiansmooth3d(image)\n\ngaussiansmooth3d(image, σ=[5,5,5];\n    mask=nothing,\n    nbox=ifelse(isnothing(mask), 3, 6), \n    weight=nothing, dims=1:min(ndims(image),3), \n    boxsizes=getboxsizes.(σ, nbox)\n    )\n\nPerforms Gaussian smoothing on image with σ as standard deviation of the Gaussian. By application of nbox times running average filters in each dimension. The length of σ and the length of the dims that are smoothed have to match. (Default 3)\n\nOptional arguments:\n\nmask: Smoothing can be performed using a mask to inter-/extrapolate missing values.\nnbox: Number of box applications. Default is 3 for normal smoothing and 6 for masked smoothing.\nweight: Apply weighted smoothing. Either weighted or masked smoothing can be porformed.\ndims: Specify which dims should be smoothed. Corresponds to manually looping of the other dimensions.\nboxizes: Manually specify the boxsizes, not using the provided σ. length(boxsizes)==length(dims) && length(boxsizes[1])==nbox\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.gaussiansmooth3d_phase","page":"Home","title":"MriResearchTools.gaussiansmooth3d_phase","text":"gaussiansmooth3d_phase(phase, σ=[5,5,5]; weight=1, kwargs...)\n\nSmoothes the phase via complex smoothing. A weighting image can be given. The same keyword arguments are supported as in gaussiansmooth3d:\n\ngaussiansmooth3d(image)\n\ngaussiansmooth3d(image, σ=[5,5,5];\n    mask=nothing,\n    nbox=ifelse(isnothing(mask), 3, 6), \n    weight=nothing, dims=1:min(ndims(image),3), \n    boxsizes=getboxsizes.(σ, nbox)\n    )\n\nPerforms Gaussian smoothing on image with σ as standard deviation of the Gaussian. By application of nbox times running average filters in each dimension. The length of σ and the length of the dims that are smoothed have to match. (Default 3)\n\nOptional arguments:\n\nmask: Smoothing can be performed using a mask to inter-/extrapolate missing values.\nnbox: Number of box applications. Default is 3 for normal smoothing and 6 for masked smoothing.\nweight: Apply weighted smoothing. Either weighted or masked smoothing can be porformed.\ndims: Specify which dims should be smoothed. Corresponds to manually looping of the other dimensions.\nboxizes: Manually specify the boxsizes, not using the provided σ. length(boxsizes)==length(dims) && length(boxsizes[1])==nbox\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.getHIP-Tuple{Any, Any}","page":"Home","title":"MriResearchTools.getHIP","text":"getHIP(mag, phase; echoes=[1,2])\n\ngetHIP(compl; echoes=[1,2])\n\nCalculates the Hermitian Inner Product between the specified echoes.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.getsensitivity","page":"Home","title":"MriResearchTools.getsensitivity","text":"getsensitivity(mag; σ, nbox=15)\n\ngetsensitivity(mag, pixdim; σ_mm=7, nbox=15)\n\ngetsensitivity(mag::NIVolume, datatype=eltype(mag); σ_mm=7, nbox=15)\n\nCalculates the bias field using the boxsegment approach. It assumes that there is a \"main tissue\" that is present in most areas of the object. Published in CLEAR-SWI.\n\nSee also makehomogeneous\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.header-Tuple{NIfTI.NIVolume}","page":"Home","title":"MriResearchTools.header","text":"header(v::NIVolume)\n\nReturns a copy of the header with the orientation information.\n\nExamples\n\njulia> vol = readmag(\"image.nii\")\njulia> hdr = header(vol)\njulia> savenii(vol .+ 10, \"vol10.nii\"; header=hdr)\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.homodyne","page":"Home","title":"MriResearchTools.homodyne","text":"homodyne(mag, phase)\n\nhomodyne(mag, phase; dims, σ)\n\nhomodyne(I; kwargs...)\n\nPerforms homodyne filtering via division of complex complex smoothing.\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.homodyne!","page":"Home","title":"MriResearchTools.homodyne!","text":"homodyne(mag, phase)\n\nhomodyne(mag, phase; dims, σ)\n\nhomodyne(I; kwargs...)\n\nPerforms homodyne filtering via division of complex complex smoothing.\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.makehomogeneous","page":"Home","title":"MriResearchTools.makehomogeneous","text":"makehomogeneous(mag::NIVolume; σ_mm=7, nbox=15)\n\nHomogeneity correction for NIVolume from NIfTI files.\n\nKeyword arguments:\n\nσ_mm: σ size for smoothing to obtain bias field. Takes NIfTI voxel size into account\nnbox: Number of boxes in each dimension for the box-segmentation step.\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.makehomogeneous-2","page":"Home","title":"MriResearchTools.makehomogeneous","text":"makehomogeneous(mag; σ, nbox=15)\n\nHomogeneity correction of 3D arrays. 4D volumes are corrected using the first 3D volume to obtain the bias field.\n\nKeyword arguments:\n\nσ: σ size in voxel for smoothing to obtain bias field. (mandatory)\nnbox: Number of boxes in each dimension for the box-segmentation step.\n\nLarger σ-values make the bias field smoother, but might not be able to catch the inhomogeneity. Smaller values can catch fast varying inhomogeneities but new inhomogeneities might be created. The stronger the bias field, the more boxes are required for segmentation. With too many boxes, it can happen that big darker structures are captured and appear overbrightened.\n\nCalculates the bias field using the boxsegment approach. It assumes that there is a \"main tissue\" that is present in most areas of the object. Published in CLEAR-SWI.\n\nSee also getsensitivity\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.makehomogeneous!","page":"Home","title":"MriResearchTools.makehomogeneous!","text":"makehomogeneous(mag; σ, nbox=15)\n\nHomogeneity correction of 3D arrays. 4D volumes are corrected using the first 3D volume to obtain the bias field.\n\nKeyword arguments:\n\nσ: σ size in voxel for smoothing to obtain bias field. (mandatory)\nnbox: Number of boxes in each dimension for the box-segmentation step.\n\nLarger σ-values make the bias field smoother, but might not be able to catch the inhomogeneity. Smaller values can catch fast varying inhomogeneities but new inhomogeneities might be created. The stronger the bias field, the more boxes are required for segmentation. With too many boxes, it can happen that big darker structures are captured and appear overbrightened.\n\nCalculates the bias field using the boxsegment approach. It assumes that there is a \"main tissue\" that is present in most areas of the object. Published in CLEAR-SWI.\n\nSee also getsensitivity\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.mask_from_voxelquality","page":"Home","title":"MriResearchTools.mask_from_voxelquality","text":"mask_from_voxelquality(qmap::AbstractArray, threshold=0.5)\n\nCreates a mask from a quality map. Another option is to use robustmask(qmap)\n\nExamples\n\njulia> qmap = romeovoxelquality(phase_3echo; TEs=[1,2,3]);\njulia> mask = mask_from_voxelquality(qmap);\n\nSee also romeovoxelquality, romeo, robustmask\n\n\n\n\n\n","category":"function"},{"location":"#MriResearchTools.mcpc3ds-Tuple{AbstractArray{<:Real}}","page":"Home","title":"MriResearchTools.mcpc3ds","text":"mcpc3ds(phase, mag; TEs, keyargs...)\n\nmcpc3ds(compl; TEs, keyargs...)\n\nmcpc3ds(phase; TEs, keyargs...)\n\nPerform MCPC-3D-S coil combination and phase offset removal on 4D (multi-echo) and 5D (multi-echo, uncombined) input.\n\nOptional Keyword Arguments\n\nechoes: only use the defined echoes. default: echoes=[1,2]\nσ: smoothing parameter for phase offsets. default: σ=[10,10,5]\nbipolar_correction: removes linear phase artefact. default: bipolar_correction=false\npo: phase offsets are stored in this array. Can be used to retrieve phase offsets or work with memory mapping.\n\nExamples\n\njulia> phase = readphase(\"phase5D.nii\")\njulia> mag = readmag(\"mag5D.nii\")\njulia> combined = mcpc3ds(phase, mag; TEs=[4,8,12])\n\nFor very large files that don't fit into memory, the uncombined data can be processed with memory mapped to disk:\n\njulia> phase = readphase(\"phase5D.nii\"; mmap=true)\njulia> mag = readmag(\"mag5D.nii\"; mmap=true)\njulia> po_size = (size(phase)[1:3]..., size(phase,5))\njulia> po = write_emptynii(po_size, \"po.nii\")\njulia> combined = mcpc3ds(phase, mag; TEs=[4,8,12], po)\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.r2s_from_t2s-Tuple{Any}","page":"Home","title":"MriResearchTools.r2s_from_t2s","text":"r2s_from_t2s(t2s) = 1000 ./ t2s\n\nConverts from T2* [ms] to R2* [1/s] values.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.readmag-Tuple{Any}","page":"Home","title":"MriResearchTools.readmag","text":"readmag(filename; rescale=false, keyargs...)\n\nReads the NIfTI magnitude with sanity checking and optional rescaling to [0;1].\n\nExamples\n\njulia> mag = readmag(\"mag.nii\")\n\nOptional keyargs are forwarded to niread:\n\nNo documentation found.\n\nNIfTI.niread is a Function.\n\n# 1 method for generic function \"niread\":\n[1] niread(file::AbstractString; mmap, mode) in NIfTI at /home/runner/.julia/packages/NIfTI/5JvZx/src/NIfTI.jl:326\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.readphase-Tuple{Any}","page":"Home","title":"MriResearchTools.readphase","text":"readphase(filename; rescale=true, keyargs...)\n\nReads the NIfTI phase with sanity checking and optional rescaling to [-π;π].\n\nExamples\n\njulia> phase = readphase(\"phase.nii\")\n\nOptional keyargs are forwarded to niread:\n\nNo documentation found.\n\nNIfTI.niread is a Function.\n\n# 1 method for generic function \"niread\":\n[1] niread(file::AbstractString; mmap, mode) in NIfTI at /home/runner/.julia/packages/NIfTI/5JvZx/src/NIfTI.jl:326\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.robustmask!-Tuple{Any}","page":"Home","title":"MriResearchTools.robustmask!","text":"robustmask!(image)\nrobustmask!(image; maskedvalue)\n\nCreates a mask and applies it inplace. It assumes that at least one corner is without signal and only contains noise.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.robustmask-Tuple{AbstractArray}","page":"Home","title":"MriResearchTools.robustmask","text":"robustmask(weight::AbstractArray)\n\nCreates a mask from a weights images. It assumes that at least one corner is without signal and only contains noise.\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.robustrescale-Tuple{Any, Any, Any}","page":"Home","title":"MriResearchTools.robustrescale","text":"robustrescale(array, newmin, newmax; threshold=false, mask=trues(size(array)), datatype=Float64)\n\nRescales the image to the the new range, disregarding outliers. Only values inside mask are used for estimating the rescaling option\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.romeovoxelquality-Tuple{Any}","page":"Home","title":"MriResearchTools.romeovoxelquality","text":"romeovoxelquality(phase::AbstractArray; keyargs...)\n\nCalculates a quality for each voxel. The voxel quality can be used to create a mask.\n\nExamples\n\njulia> qmap = romeovoxelquality(phase_3echo; TEs=[1,2,3]);\njulia> mask1 = mask_from_voxelquality(qmap);\njulia> mask2 = robustmask(qmap);\n\nTakes the same inputs as romeo/unwrap:\n\nunwrap(wrapped::AbstractArray; keyargs...)\n\nROMEO unwrapping for 3D and 4D data.\n\nOptional keyword arguments:\n\nTEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\nweights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\nmag: The magnitude is used to improve the unwrapping-path.\nmask: Unwrapping is only performed inside the mask.\nphase2: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\ncorrectglobal=false: If true corrects global n2π offsets.\nindividual=false: If true perform individual unwrapping of echos.   Type ?unwrap_individual for more information\ntemplate=2: echo that is spatially unwrapped (if individual is false)\nmaxseeds=1: higher values allow more seperate regions\nmerge_regions=false: spatially merge neighboring regions after unwrapping\ncorrect_regions=false: bring each regions median closest to 0 by adding n2π\nwrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\ntemporal_uncertain_unwrapping=false: uses spatial unwrapping on voxels that   have high uncertainty values after temporal unwrapping\n\nExamples\n\njulia> using MriResearchTools\njulia> phase = readphase(\"phase_3echo.nii\")\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\njulia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase))\n\nSee also mask_from_voxelquality, romeo, robustmask\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.savenii-Tuple{AbstractArray, Any}","page":"Home","title":"MriResearchTools.savenii","text":"savenii(image::AbstractArray, filepath; header=nothing)\n\nsavenii(image::AbstractArray, name, writedir, header=nothing)\n\nWarning: MRIcro can only open images with types Int32, Int64, Float32, Float64\n\nExamples\n\njulia> savenii(ones(64,64,5), \"image.nii\")\n\njulia> savenii(ones(64,64,5), \"image2\", \"folder\")\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.to_dim-Tuple{Real, Int64}","page":"Home","title":"MriResearchTools.to_dim","text":"to_dim(V::AbstractVector, dim::Int)\n\nto_dim(a::Real, dim::Int)\n\nConverts a vector or number to a higher dimension.\n\nExamples\n\njulia> to_dim(5, 3)\n1×1×1 Array{Int64, 3}:\n[:, :, 1] =\n 5\njulia> to_dim([1,2], 2)\n 1×2 Matrix{Int64}:\n  1  2\n\n\n\n\n\n","category":"method"},{"location":"#MriResearchTools.write_emptynii-Tuple{Any, Any}","page":"Home","title":"MriResearchTools.write_emptynii","text":"write_emptynii(size, path; datatype=Float32, header=NIVolume(zeros(datatype, 1)).header)\n\nWrites an empty NIfTI image to disk that can be used for memory-mapped access.\n\nExamples\n\njulia> vol = write_emptynii((64,64,64), \"empty.nii\")\njulia> vol.raw[:,:,1] .= ones(64,64) # synchronizes mmapped file on disk\n\nWarning: MRIcro can only open images with types Int32, Int64, Float32, Float64\n\n\n\n\n\n","category":"method"}]
}
