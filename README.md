# MriResearchTools

[![Build Status](https://travis-ci.com/korbinian90/MriResearchTools.jl.svg?branch=master)](https://travis-ci.com/korbinian90/MriResearchTools.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/MriResearchTools.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/MriResearchTools-jl)
[![Codecov](https://codecov.io/gh/korbinian90/MriResearchTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/MriResearchTools.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/MriResearchTools.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/MriResearchTools.jl?branch=master)

### Prerequisites
A Julia installation v1.x is required.

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing
Open the Julia REPL and type

```julia
julia> ] # enter julia package manager
(v1.3) pkg> add MriResearchTools
(v1.3) pkg> # type backspace to get back to the julia REPL
julia>
```

### Included Functionality

ROMEO 3D/4D Phase Unwrapping

Reading and writing NIfTI files (adapted from JuliaIO/NIfTI)

Magnitude homogeneity correction ([example](https://github.com/korbinian90/Magnitude-Intensity-Correction/blob/master/Intensity%20Correction.ipynb))

Simple robust masking (threshold)

Combine multiple echoes

Laplacian unwrapping

Unwarping of B0 dependent shifts

Fast gaussian smoothing
- standard
- weighted
- with missing values

TODO: Tests and Examples
