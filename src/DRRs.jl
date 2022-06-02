module DRRs

# Utilities
export read_dicom
include("io.jl")

export Camera, Detector, make_plane, get_rays, trace
include("camera.jl")

# Sampling
export make_coordinate_matrix, make_inverse_coordinate_matrix, interpolate
include("sampling/trilinear.jl")

end
