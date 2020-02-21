# MriResearchTools

[![Build Status](https://travis-ci.com/korbinian90/MRI.jl.svg?branch=master)](https://travis-ci.com/korbinian90/MRI.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/korbinian90/MRI.jl?svg=true)](https://ci.appveyor.com/project/korbinian90/MRI-jl)
[![Codecov](https://codecov.io/gh/korbinian90/MRI.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/korbinian90/MRI.jl)
[![Coveralls](https://coveralls.io/repos/github/korbinian90/MRI.jl/badge.svg?branch=master)](https://coveralls.io/github/korbinian90/MRI.jl?branch=master)

### Prerequisites
A Julia installation v1.1 or higher

Magnitude and Phase images in NIfTI fileformat (4D images with echoes in the 4th dimension)

### Installing
Currently, MRI depends on a not yet registered version of NIfTI.jl

Open the REPL in Julia and type

```julia
import Pkg;
Pkg.add(PackageSpec(name="NIfTI", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/korbinian90/MRI.jl"))
```

### Included Functionality

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
