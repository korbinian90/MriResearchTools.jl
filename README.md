# MriResearchTools - Version matching the CLEAR SWI publication

## Prerequisites
A Julia installation v1.x is required.

## Quick Start
### SNR Plots
Generate the SNR graph plot and duty cycle table appearing in the CLEAR SWI publication:

```julia
julia> ] add MriResearchTools#clearswi_publication
julia> using MriResearchTools
julia> snr_graph()
julia> dutycycle_table()
```
The first line downloads and installs the package (with dependencies). To leave the julia package manager (shell prompt `(v1.3) pkg>`) press `backspace`.

### Files
`scanner.jl` and `tissue.jl` contain constants

`snr.jl` contains the SNR theory and simulation

`snr_plots.jl` contains the settings and functions to create the graph and table

## License
This project is licensed under the MIT License - see the [LICENSE](https://github.com/korbinian90/MriResearchTools.jl/blob/master/LICENSE) for details
