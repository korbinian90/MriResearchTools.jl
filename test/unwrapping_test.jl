phasefile = joinpath("data", "small", "Phase.nii")
magfile = joinpath("data", "small", "Mag.nii")
phase = Float32.(readphase(phasefile))
magni = readmag(magfile)

@test sizeof(unwrap(phase; mag=magni)) == sizeof(phase)
@test sizeof(unwrap_individual(phase; mag=magni)) == sizeof(phase)

#=
outpath = raw"F:\MRI\debug\romeo\simon20200130"
phasefile = joinpath(outpath, "combined_phase_not_nan.nii")
phase = Float32.(niread(phasefile))

savenii(unwrap(phase), "unwrapped.nii", outpath)


phasefile = raw"F:\MRI\scanner_nifti\Paper\SWI_paper_7T_volunteers\19930503JSPC_201907041530\nifti\4\reform\Image.nii"
magfile = raw"F:\MRI\scanner_nifti\Paper\SWI_paper_7T_volunteers\19930503JSPC_201907041530\nifti\3\reform\Image.nii"
phase = Float32.(readphase(phasefile))
magni = readmag(magfile)
mag = Float32.(magni)

unwrapped = unwrap_individual(phase; mag=mag)
hdr = magni.header
hdr.scl_inter = 0
hdr.scl_slope = 1
savenii(unwrapped, raw"F:\MRI\Analysis\Volunteer7T\phase\romeo\unwrapped.nii"; header=hdr)

unwrapped[:,:,:,1] .+= 50
unwrapped[:,:,:,1] .-= 50

TEs = reshape(1:size(phase,4),1,1,1,:)
B0 = 1000 * sum(unwrapped .* mag; dims=4)
B0 ./= sum(mag .* TEs; dims=4)

savenii(B0, raw"F:\MRI\Analysis\Volunteer7T\phase\romeo\B0_50.nii"; header=hdr)

savenii(mag[:,:,:,1], raw"F:\MRI\Analysis\Volunteer7T\phase\romeo\mag_1.nii"; header=hdr)

M1 = copy(mag)
M1[:,:,:,1] .= 5000

MC = sum(M1 .* mag .* TEs; dims=4)
MC ./= sum(M1 .* TEs; dims=4)
savenii(MC, raw"F:\MRI\Analysis\Volunteer7T\phase\romeo\MC_5000_T.nii"; header=hdr)
=#
