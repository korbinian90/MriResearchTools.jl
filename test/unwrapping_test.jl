phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
magni = readmag(magfile)

@test sizeof(unwrap(phase; mag=magni)) == sizeof(phase)
@test sizeof(unwrap_individual(phase; mag=magni)) == sizeof(phase)
