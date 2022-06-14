module DRRs

# Utilities
export CT, read_dicom
include("io.jl")

export Camera, Detector, make_plane, get_rays, trace
include("camera.jl")

# Sampling
export make_coordinate_matrix, make_inverse_coordinate_matrix, trilinear
include("sampling/trilinear.jl")

export sample
include("sampling/sampler.jl")

export DRR
include("drr.jl")

end