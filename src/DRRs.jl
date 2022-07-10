module DRRs

# Utilities
export CT, read_dicom
include("io.jl")

export Camera, Detector, Ray, make_xrays, trace, length
include("camera.jl")

export get_basis
include("utils.jl")

# Sampling
export make_coordinate_matrix, make_inverse_coordinate_matrix, trilinear
include("sampling/trilinear.jl")

export sample
include("sampling/sampler.jl")

export siddon
include("sampling/siddon.jl")

export sampler, siddon
include("drr.jl")

# Visualization
export plot, plot_image
include("viz.jl")

end