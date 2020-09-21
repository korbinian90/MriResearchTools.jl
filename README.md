# MriResearchTools

[![Build Status](https://travis-ci.com/korbinian90/MriResearchTools.jl.svg?branch=master)](https://travis-ci.com/korbinian90/MriResearchTools.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/MriResearchTools.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/MriResearchTools-jl)
[![Codecov](https://codecov.io/gh/korbinian90/MriResearchTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/MriResearchTools.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/MriResearchTools.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/MriResearchTools.jl?branch=master)

## Prerequisites
A Julia installation v1.x is required.

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

## Installing
Open the Julia REPL and type

```julia
julia> ] # enter julia package manager
(v1.5) pkg> add MriResearchTools
(v1.5) pkg> # type backspace to get back to the julia REPL
julia>
```

## Quick Start
Open multi-echo 4D NIfTI phase and magnitude files and perform ROMEO phase unwrapping.

```julia
using MriResearchTools
# input images
TEs = [4,8,12]
nifti_folder = joinpath("test", "data", "small")
magfile = joinpath(nifti_folder, "Mag.nii") # Path to the magnitude image in nifti format, must be .nii or .hdr
phasefile = joinpath(nifti_folder, "Phase.nii") # Path to the phase image
# load images
mag = readmag(magfile)
phase = readphase(phasefile)
# unwrap
unwrapped = romeo(phase; mag=mag, TEs=TEs)
# save unwrapped image
outputfolder = "outputFolder"
mkpath(outputfolder)
savenii(unwrapped, "unwrapped", outputfolder, header(phase))
```

## Included Functionality

[ROMEO](https://github.com/korbinian90/ROMEO.jl) 3D/4D Phase Unwrapping  
`romeo` `unwrap` `unwrap_individual`

Laplacian unwrapping  
`laplacianunwrap`

MCPC3Ds multi-echo coil combination (at the moment only monopolar)  
`mcpc3ds`

Reading, writing and other functions for NIfTI files (adapted from JuliaIO/NIfTI)  
`readphase` `readmag` `niread` `savenii` `header` `write_emptynii`

Magnitude homogeneity correction ([example](https://github.com/korbinian90/Magnitude-Intensity-Correction/blob/master/Intensity%20Correction.ipynb))  
`makehomogeneous`

Simple robust masking (threshold)  
`robustmask`

Combine multiple echoes  
`RSS`

Unwarping of B0 dependent shifts  
`getVSM` `thresholdforward` `unwarp`

Fast gaussian smoothing for real and complex data  
`gaussiansmooth3d`
  - standard
  - weighted
  - with missing values

Other functions
`robustrescale` `getHIP` `getsensitivity` `getscaledimage` `estimatequantile` `estimatenoise`

## Publications
### ROMEO
Dymerska, B., Eckstein, K., Bachrata, B., Siow, B., Trattnig, S., Shmueli, K., Robinson, S.D., 2020. Phase Unwrapping with a Rapid Opensource Minimum Spanning TreE AlgOrithm (ROMEO). bioRxiv 2020.07.24.214551. https://doi.org/10.1101/2020.07.24.214551

### MCPC3Ds 
Eckstein, K., Dymerska, B., Bachrata, B., Bogner, W., Poljanc, K., Trattnig, S., Robinson, S.D., 2018. Computationally Efficient Combination of Multi-channel Phase Data From Multi-echo Acquisitions (ASPIRE). Magnetic Resonance in Medicine 79, 2996–3006. https://doi.org/10.1002/mrm.26963

### Homogeneity Correction
Eckstein, K., Trattnig, S., Robinson, S.D., 2019. A Simple Homogeneity Correction for Neuroimaging at 7T, in: Proceedings of the 27th Annual Meeting ISMRM. Presented at the ISMRM, Montréal, Québec, Canada. https://index.mirasmart.com/ISMRM2019/PDFfiles/2716.html


## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/MriResearchTools.jl/blob/master/LICENSE) for details


## TODO
more Tests, Examples and Documentation
