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
(v1.3) pkg> add MriResearchTools
(v1.3) pkg> # type backspace to get back to the julia REPL
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

ROMEO 3D/4D Phase Unwrapping\
`romeo` `unwrap` `unwrap_individual`

Reading, writing and other functions for NIfTI files (adapted from JuliaIO/NIfTI)\
`readphase` `readmag` `niread` `savenii` `header`

Magnitude homogeneity correction ([example](https://github.com/korbinian90/Magnitude-Intensity-Correction/blob/master/Intensity%20Correction.ipynb))\
`makehomogeneous`

Simple robust masking (threshold)\
`robustmask`

Combine multiple echoes\
`RSS`

Laplacian unwrapping\
`laplacianunwrap`

Unwarping of B0 dependent shifts\
`getVSM` `thresholdforward` `unwarp`

Fast gaussian smoothing\
`gaussiansmooth3d`
  - standard
  - weighted
  - with missing values

Other functions
`robustrescale` `getHIP` `getsensitivity` `getscaledimage` `estimatequantile` `estimatenoise`

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/MriResearchTools.jl/blob/master/LICENSE) for details


## TODO
Tests, Examples and Documentation
