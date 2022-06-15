module DRRs

# Utilities
export CT, read_dicom
include("io.jl")

export Camera, Detector, make_xrays, trace
include("camera.jl")

export get_basis
include("utils.jl")

# Sampling
export make_coordinate_matrix, make_inverse_coordinate_matrix, trilinear
include("sampling/trilinear.jl")

export sample
include("sampling/sampler.jl")

export DRR
include("drr.jl")

# Visualization
export plot, plot_image
include("viz.jl")

end